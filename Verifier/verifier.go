package verifier

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"

	Parameters "vSIS-Signature/System"
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
	if err != nil {
		log.Fatalf("ToPolyNTT: %v", err)
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
	//     – ℓ₂ on each polynomial row  s_i
	//     – ∞-over-2  :  max_i ‖s_i‖₂
	//     – 2-over-2  :  √Σ_i ‖s_i‖₂²
	//
	// ---------------------------------------------------------------------------
	var (
		globalInfOverL2 float64 // max_i ‖s_i‖₂        (∞ over 2)
		sumRowL2Sq      float64 // Σ_i  ‖s_i‖₂²  → √   (2 over 2)
	)

	abs64 := func(x int64) int64 {
		if x < 0 {
			return -x
		}
		return x
	}

	for i, p := range SEval {
		coef := ringQ.NewPoly()
		ringQ.InvNTT(p, coef) // back to coefficients

		var (
			rowL2Sq float64
			rowInf  int64
		)

		for _, c := range coef.Coeffs[0] {
			// centre coefficient in (−q/2 , q/2]
			var v int64
			if c > pp.Q/2 {
				v = int64(c) - int64(pp.Q)
			} else {
				v = int64(c)
			}

			rowL2Sq += float64(v * v)
			if a := abs64(v); a > rowInf {
				rowInf = a
			}
		}

		//-----------------------------------------------------------------------
		//  Row-level coefficient bound :  ‖s_i‖_∞ ≤ q
		//-----------------------------------------------------------------------
		if rowInf > int64(pp.Q) {
			log.Printf("❌ row %d : ‖s_i‖_∞ = %d  >  q = %d", i, rowInf, pp.Q)
			return false
		}

		rowL2 := math.Sqrt(rowL2Sq)
		log.Printf("row %2d : ‖s_i‖₂ = %.4f   ‖s_i‖_∞ = %d", i, rowL2, rowInf)

		//-----------------------------------------------------------------------
		//  Accumulate global norms
		//-----------------------------------------------------------------------
		if rowL2 > globalInfOverL2 {
			globalInfOverL2 = rowL2
		}
		sumRowL2Sq += rowL2Sq
	}

	// ---------------------------------------------------------------------------
	//
	//	Global figures
	//
	// ---------------------------------------------------------------------------
	globalL2OverL2 := math.Sqrt(sumRowL2Sq) // √Σ_i ‖s_i‖₂²

	log.Printf("GLOBAL :  max_i‖s_i‖₂ = %.4f    √Σ‖s_i‖₂² = %.4f",
		globalInfOverL2, globalL2OverL2)

	// ---------------------------------------------------------------------------
	//
	//	Global bound  (still using bound = q  on the  ∞-over-2  value)
	//
	// ---------------------------------------------------------------------------
	if globalInfOverL2 > float64(pp.Q)*float64(pp.Q) {
		log.Printf("❌ global bound failed: max_i‖s_i‖₂ = %.4f  >  q**2 = %d",
			globalInfOverL2, pp.Q*pp.Q)
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
