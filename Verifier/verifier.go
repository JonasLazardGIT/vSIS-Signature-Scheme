package verifier

import (
	"encoding/json"
	"fmt"
	"log"
	"os"

	ps "vSIS-Signature/Preimage_Sampler"
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

	// 4) load B‐matrix as CyclotomicFieldElem[]
	const prec = 256
	Bcyclo, err := loadBMatrix("Parameters/Bmatrix.json", ringQ, prec)
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
	tEval, err := vsishash.ComputeBBSHash(ringQ, Bcyclo, mEval, x0Eval, x1Eval, prec)
	if err != nil {
		log.Fatalf("ComputeBBSHash: %v", err)
	}
	tCmp, err := vsishash.ToPolyNTT(tEval, ringQ)
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

	// // 12) norm bound
	// var sumSq float64
	// for i, p := range SEval {
	// 	coef := ringQ.NewPoly()
	// 	ringQ.InvNTT(p, coef)
	// 	for j, c := range coef.Coeffs[0] {
	// 		var v int64
	// 		if c > pp.Q/2 {
	// 			v = int64(c) - int64(pp.Q)
	// 		} else {
	// 			v = int64(c)
	// 		}
	// 		sumSq += float64(v * v)
	// 		if sumSq > pp.Bound*pp.Bound {
	// 			log.Printf("❌ norm-check failed: partial norm²=%.2f > bound²=%.2f at s[%d][%d]", sumSq, pp.Bound*pp.Bound, i, j)
	// 			return false
	// 		}
	// 	}
	// }
	// if math.Sqrt(sumSq) > pp.Bound {
	// 	log.Printf("❌ norm-check failed: ‖s‖=%.2f > %.2f", math.Sqrt(sumSq), pp.Bound)
	// 	return false
	// }

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

func loadBMatrix(path string, ringQ *ring.Ring, prec uint) ([]*ps.CyclotomicFieldElem, error) {
	type bjson struct {
		B [][]uint64 `json:"B"`
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var bj bjson
	if err := json.Unmarshal(raw, &bj); err != nil {
		return nil, err
	}
	if len(bj.B) != 4 {
		return nil, fmt.Errorf("expected 4 polys in B, got %d", len(bj.B))
	}
	Bcyclo := make([]*ps.CyclotomicFieldElem, 4)
	for i := 0; i < 4; i++ {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], bj.B[i])
		Bcyclo[i] = ps.ConvertFromPolyBig(ringQ, p, prec)
	}
	return Bcyclo, nil
}
