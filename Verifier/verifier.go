package verifier

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"math/big"
	"os"
	"time"

	Parameters "vSIS-Signature/System"
	measure "vSIS-Signature/measure"
	prof "vSIS-Signature/prof"
	vsishash "vSIS-Signature/vSIS-HASH"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// PublicKey mirrors public_key/public_key.json
type PublicKey struct {
	A    [][]uint64 `json:"A"`
	Base uint64     `json:"base"`
	K    int        `json:"k"`
}

// SignatureData mirrors Signature/Signature.json, now including the target syndrome
type SignatureData struct {
	Message   []uint64   `json:"message"`   // In Coeffs
	X0        []uint64   `json:"x0"`        // In Coeffs
	X1        []uint64   `json:"x1"`        // In Coeffs
	Target    []uint64   `json:"target"`    // t (syndrome) in NTT domain
	Signature [][]uint64 `json:"signature"` // each row already NTT
}

func Verify() bool {
	defer prof.Track(time.Now(), "Verify")
	// 1) load system params
	pp, err := loadParams("Parameters/Parameters.json")
	if err != nil {
		log.Fatalf("loading params: %v", err)
	}

	// 2) instantiate R_q
	ringQ, err := ring.NewRing(pp.N, []uint64{pp.Q})
	if err != nil {
		log.Fatalf("ring.NewRing: %v", err)
	}
	if measure.Enabled {
		q := new(big.Int).SetUint64(pp.Q)
		bytesF := measure.BytesField(q)
		bq := q.BitLen()
		fmt.Printf("[measure] Params: q=%d bq=%d BytesF=%d phi=%d\n", pp.Q, bq, bytesF, pp.N)
	}

	// 3) load public key
	pk, err := loadPublicKey("public_key/public_key.json")
	if err != nil {
		log.Fatalf("loading public key: %v", err)
	}

	// 4) load B‐matrix
	const prec = 256
	Bcyclo, err := loadBMatrix("Parameters/Bmatrix.json", ringQ)
	if err != nil {
		log.Fatalf("loading Bmatrix: %v", err)
	}

	// 5) load signature bundle (including target)
	sig, err := loadSignature("Signature/Signature.json")
	if err != nil {
		log.Fatalf("loading signature: %v", err)
	}
	if measure.Enabled {
		q := new(big.Int).SetUint64(pp.Q)
		bytesF := measure.BytesField(q)
		measure.Global.Add("verify/message", int64(len(sig.Message)*bytesF))
		measure.Global.Add("verify/x0", int64(len(sig.X0)*bytesF))
		measure.Global.Add("verify/x1", int64(len(sig.X1)*bytesF))
		measure.Global.Add("verify/target", int64(len(sig.Target)*bytesF))
		bytesR := measure.BytesRing(pp.N, q)
		measure.Global.Add("verify/signature", int64(len(sig.Signature))*int64(bytesR))
		fmt.Printf("[measure] Loaded signature: |u|=%s |x0|=%s |x1|=%s |s|=%s\n",
			measure.Human(int64(len(sig.Message)*bytesF)),
			measure.Human(int64(len(sig.X0)*bytesF)),
			measure.Human(int64(len(sig.X1)*bytesF)),
			measure.Human(int64(len(sig.Signature))*int64(bytesR)))
	}

	// b₀ check: lengths
	if len(sig.Message) != pp.N || len(sig.X0) != pp.N || len(sig.X1) != pp.N || len(sig.Target) != pp.N {
		log.Println("❌ b₀-check failed: incorrect lengths in signature data")
		return false
	}

	// 6) rebuild mEval, x0Eval, x1Eval
	mEval := ringQ.NewPoly()
	copy(mEval.Coeffs[0], sig.Message)
	x0Eval := ringQ.NewPoly()
	copy(x0Eval.Coeffs[0], sig.X0)
	x1Eval := ringQ.NewPoly()
	copy(x1Eval.Coeffs[0], sig.X1)

	// 7) recompute t = BBSHash
    tCmp, err := vsishash.ComputeBBSHash(ringQ, Bcyclo, mEval, x0Eval, x1Eval)
    if err != nil {
        log.Fatalf("ComputeBBSHash: %v", err)
    }

	// 8) load target syndrome from signature
	targetPoly := ringQ.NewPoly()
	copy(targetPoly.Coeffs[0], sig.Target)

	// compare recomputed vs target
	if !ringQ.Equal(tCmp, targetPoly) {
		log.Println("❌ target-check failed: recomputed hash ≠ target stored in signature")
		return false
	}
	log.Println("✅ target-check passed: recomputed hash matches target")

	// 9) build Aeval (NTT-domain)
	Aeval := make([]*ring.Poly, len(pk.A))
	for i, coeffs := range pk.A {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs)
		Aeval[i] = p
	}

	// 10) build SEval
	if len(sig.Signature) != len(Aeval) {
		log.Fatalf("signature length %d, expected %d", len(sig.Signature), len(Aeval))
	}
	SEval := make([]*ring.Poly, len(sig.Signature))
	for i, coeffs := range sig.Signature {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs)
		SEval[i] = p
	}

	// 11) check A·s == target
	acc := ringQ.NewPoly()
	tmp := ringQ.NewPoly()
	for i, si := range SEval {
		ringQ.MulCoeffs(Aeval[i], si, tmp)
		ringQ.Add(acc, tmp, acc)
	}
	if !ringQ.Equal(acc, targetPoly) {
		log.Println("❌ b₁-check failed: A·s ≠ target")
		return false
	}
	// ---------------------------------------------------------------------------
	//  12. Norm accounting
	//     – ℓ₂ on each polynomial row  s_i, using
	//       dist(c) = min{|c - (−q/2)|, |c - 0|, |c - (+q/2)|}
	//     – ∞-over-2  :  max_i √Σ_j dist(c_{i,j})²
	//     – 2-over-2  :  √Σ_i Σ_j dist(c_{i,j})²
	// ---------------------------------------------------------------------------
	var (
		globalInfOverL2 float64 // max_i ‖s_i‖₂
		sumRowL2Sq      float64 // Σ_i  ‖s_i‖₂²
	)

	qval := int64(pp.Q)
	halfQ := qval / 2
	centers := []int64{-halfQ, 0, halfQ}

	abs64 := func(x int64) int64 {
		if x < 0 {
			return -x
		}
		return x
	}

	for i, p := range SEval {
		// 1) bring p back to coefficient domain
		coefPoly := ringQ.NewPoly()
		ringQ.InvNTT(p, coefPoly)

		var (
			rowL2Sq float64
			rowInf  int64
		)

		// 2) for each coefficient, compute its minimal “distance”
		for _, cu := range coefPoly.Coeffs[0] {
			// first center into (−q/2, +q/2]
			var v int64
			if cu > uint64(halfQ) {
				v = int64(cu) - qval
			} else {
				v = int64(cu)
			}

			// now compute dist = min_j |v - centers[j]|
			dist := abs64(v - centers[0])
			for _, c0 := range centers[1:] {
				if d := abs64(v - c0); d < dist {
					dist = d
				}
			}

			// accumulate into ℓ₂² and ∞
			rowL2Sq += float64(dist * dist)
			if dist > rowInf {
				rowInf = dist
			}
		}

		// row‐level norm
		rowL2 := math.Sqrt(rowL2Sq)
		log.Printf("row %2d : ‖s_i‖₂ (using min-distance) = %.4f", i, rowL2)

		// accumulate global norms
		if rowL2 > globalInfOverL2 {
			globalInfOverL2 = rowL2
		}
		sumRowL2Sq += rowL2Sq
	}

	// ---------------------------------------------------------------------------
	// Global figures
	// ---------------------------------------------------------------------------
	globalL2OverL2 := math.Sqrt(sumRowL2Sq) // √Σ_i ‖s_i‖₂²
	log.Printf("GLOBAL : max_i‖s_i‖₂ = %.4f    √Σ‖s_i‖₂² = %.4f",
		globalInfOverL2, globalL2OverL2)

	// ---------------------------------------------------------------------------
	// Global bound check (using q on the ∞-over-2 value)
	// ---------------------------------------------------------------------------
	if globalInfOverL2 > float64(pp.Q) {
		log.Printf("❌ global bound failed: max_i‖s_i‖₂ = %.4f  >  q = %d",
			globalInfOverL2, pp.Q)
	} else {
		log.Printf("✅ norm checks passed: max_i‖s_i‖₂ = %.4f,   √Σ‖s_i‖₂² = %.4f",
			globalInfOverL2, globalL2OverL2)
	}

	// ---------------------------------------------------------------------------
	// Extract centered coefficients into coefSigs (unchanged)
	// ---------------------------------------------------------------------------
	coefSigs := make([][]int64, len(SEval))
	for i, p := range SEval {
		coefPoly := ringQ.NewPoly()
		ringQ.InvNTT(p, coefPoly)

		row := make([]int64, pp.N)
		for j, cu := range coefPoly.Coeffs[0] {
			if cu > uint64(halfQ) {
				row[j] = int64(cu) - qval
			} else {
				row[j] = int64(cu)
			}
		}
		coefSigs[i] = row
	}
	// prepare JSON blob
	out := struct {
		Timestamp string    `json:"timestamp"` // YYYYMMDD_HHMMSS
		Coeffs    [][]int64 `json:"coeffs"`
	}{
		Timestamp: time.Now().Format("20060102_150405"),
		Coeffs:    coefSigs,
	}

	// ensure directory
	_ = os.MkdirAll("Signature_Reading", 0755)
	fname := fmt.Sprintf("Signature_Reading/coeffsig_%s.json", out.Timestamp)
	f, err := os.Create(fname)
	if err != nil {
		log.Printf("⚠️ could not write coeff‐sig JSON: %v", err)
	} else {
		defer f.Close()
		enc := json.NewEncoder(f)
		enc.SetIndent("", "  ")
		if err := enc.Encode(out); err != nil {
			log.Printf("⚠️ error encoding coeff‐sig JSON: %v", err)
		} else {
			log.Printf("✔ saved coefficient‐domain signature to %s", fname)
		}
	}

	if globalInfOverL2 > float64(pp.Q) {
		log.Printf("❌ global bound failed: max_i‖s_i‖₂ = %.4f  >  q = %d",
			globalInfOverL2, pp.Q)
		return false
	}

	log.Printf("✅ norm checks passed: max_i‖s_i‖₂ = %.4f,   √Σ‖s_i‖₂² = %.4f",
		globalInfOverL2, globalL2OverL2)

	log.Println("✅ signature valid")
	return true
}

func loadParams(path string) (*Parameters.SystemParams, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var p Parameters.SystemParams
	if err := json.Unmarshal(data, &p); err != nil {
		return nil, err
	}
	return &p, nil
}

func loadPublicKey(path string) (*PublicKey, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var pk PublicKey
	if err := json.Unmarshal(data, &pk); err != nil {
		return nil, err
	}
	return &pk, nil
}

func loadSignature(path string) (*SignatureData, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var sd SignatureData
	if err := json.Unmarshal(data, &sd); err != nil {
		return nil, err
	}
	return &sd, nil
}

// tiny helper
func abs64(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}

// ----------------------------------------------------------------------------
// loadBMatrix – now returns the four B-polynomials **already in NTT form**
// ----------------------------------------------------------------------------
func loadBMatrix(path string, ringQ *ring.Ring) ([]*ring.Poly, error) {

	// --- read JSON -----------------------------------------------------------
	type bjson struct {
		B [][]uint64 `json:"B"`
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var bj bjson
	if err = json.Unmarshal(raw, &bj); err != nil {
		return nil, err
	}
	if len(bj.B) != 4 {
		return nil, fmt.Errorf("expected 4 polys in B, got %d", len(bj.B))
	}

	// --- re-hydrate and lift to NTT ------------------------------------------
	B := make([]*ring.Poly, 4)
	for i := 0; i < 4; i++ {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], bj.B[i]) // coefficients → poly (coeff. domain)
		ringQ.NTT(p, p)            // single lift to evaluation / NTT domain
		B[i] = p
	}

	return B, nil
}
