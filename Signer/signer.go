package signer

import (
	"encoding/json"
	"fmt"
	"log"
	"os"

	ps "vSIS-Signature/Preimage_Sampler"
	Parameters "vSIS-Signature/System"
	vsishash "vSIS-Signature/vSIS-HASH"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// PublicKey holds the public part of the trapdoor.
type PublicKey struct {
	A    [][]uint64 `json:"A"`
	Base uint64     `json:"base"`
	K    int        `json:"k"`
}

// PrivateKey holds the secret part of the trapdoor.
type PrivateKey struct {
	R0   [][]uint64 `json:"R0"`
	R1   [][]uint64 `json:"R1"`
	Base uint64     `json:"base"`
	K    int        `json:"k"`
}

type SignatureData struct {
	Message   []uint64   `json:"message"`
	X0        []uint64   `json:"x0"`
	X1        []uint64   `json:"x1"`
	Target    []uint64   `json:"target"`
	Signature [][]uint64 `json:"signature"`
}

func Sign() {
	// 1) load public parameters
	params, err := loadParams("Parameters/Parameters.json")
	if err != nil {
		log.Fatalf("loadParams: %v", err)
	}
	log.Printf("Params: N=%d Q=%d Base=%d K=%d", params.N, params.Q, params.Base, params.K)

	// 2) instantiate ring R_q
	ringQ, err := ring.NewRing(params.N, []uint64{params.Q})
	if err != nil {
		log.Fatalf("ring.NewRing: %v", err)
	}

	// 3) keygen → trapdoor
	trap := ps.TrapGen(ringQ, params.Base, params.SigmaT)
	log.Println("Trapdoor generated.")

	// 4) write public and private keys
	if err := savePublicKey("public_key/public_key.json", ringQ, &trap); err != nil {
		log.Fatalf("savePublicKey: %v", err)
	}
	if err := savePrivateKey("private_key/private_key.json", ringQ, &trap); err != nil {
		log.Fatalf("savePrivateKey: %v", err)
	}

	// 5) load B‐matrix from JSON and convert to CyclotomicFieldElem
	const prec = 256
	Bcyclo, err := loadBMatrix("Parameters/Bmatrix.json", ringQ, prec)
	if err != nil {
		log.Fatalf("loadBMatrix: %v", err)
	}

	// 6) sample a random message m ∈ R_q
	prng, _ := utils.NewPRNG()
	uni := ring.NewUniformSampler(prng, ringQ)
	mPoly := ringQ.NewPoly()
	uni.Read(mPoly)
	mCoeffs := append([]uint64(nil), mPoly.Coeffs[0]...)
	// 7) sample randomness x0, x1 ∈ R_q
	x0Poly := ringQ.NewPoly()
	x1Poly := ringQ.NewPoly()
	uni.Read(x0Poly)
	uni.Read(x1Poly)
	x0Coeffs := append([]uint64(nil), x0Poly.Coeffs[0]...)
	x1Coeffs := append([]uint64(nil), x1Poly.Coeffs[0]...)

	// 8) compute BBS‐syndrome in Eval domain
	tEval, err := vsishash.ComputeBBSHash(ringQ, Bcyclo, mPoly, x0Poly, x1Poly, prec)
	if err != nil {
		log.Fatalf("ComputeBBSHash: %v", err)
	}
	tNTT, err := vsishash.ToPolyNTT(tEval, ringQ)
	if err != nil {
		log.Fatalf("ToPolyNTT: %v", err)
	}
	// extract target syndrome coefficients
	tNTTCoeffs := append([]uint64(nil), tNTT.Coeffs[0]...)
	fmt.Printf("Computed tNTT.Coeffs[0][:8]=%v\n", tNTTCoeffs[:8])
	trap.K = params.K
	// 9) Gaussian pre‐image sampling in Eval domain
	sEval := ps.GaussSamp(
		ringQ,
		trap.A, trap.R[0], trap.R[1],
		tNTT,
		params.Sigma, params.Bound,
		trap.Base, trap.K,
	)

	// 10) pull out s in NTT domain directly
	sNTT := make([][]uint64, len(sEval))
	for i, p := range sEval {
		sNTT[i] = append([]uint64(nil), p.Coeffs[0]...)
	}

	// 11) save everything in NTT domain, including target
	if err := saveSignature(
		"Signature/Signature.json",
		mCoeffs, x0Coeffs, x1Coeffs,
		tNTTCoeffs,
		sNTT,
	); err != nil {
		log.Fatalf("saveSignature: %v", err)
	}
	log.Println("✔ Signing complete")
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
	for i := range bj.B {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], bj.B[i])
		Bcyclo[i] = ps.ConvertFromPolyBig(ringQ, p, prec)
	}
	return Bcyclo, nil
}

func savePublicKey(path string, ringQ *ring.Ring, trap *ps.Trapdoor) error {
	_ = os.MkdirAll("public_key", 0755)
	pub := PublicKey{Base: trap.Base, K: trap.K, A: make([][]uint64, len(trap.A))}
	for i, a := range trap.A {
		pub.A[i] = append([]uint64(nil), a.Coeffs[0]...)
	}
	out, _ := json.MarshalIndent(pub, "", "  ")
	return os.WriteFile(path, out, 0644)
}

func savePrivateKey(path string, ringQ *ring.Ring, trap *ps.Trapdoor) error {
	_ = os.MkdirAll("private_key", 0755)
	priv := PrivateKey{
		Base: trap.Base,
		K:    trap.K,
		R0:   make([][]uint64, len(trap.R[0])),
		R1:   make([][]uint64, len(trap.R[1])),
	}
	for i, r0 := range trap.R[0] {
		priv.R0[i] = append([]uint64(nil), r0.Coeffs[0]...)
	}
	for i, r1 := range trap.R[1] {
		priv.R1[i] = append([]uint64(nil), r1.Coeffs[0]...)
	}
	out, _ := json.MarshalIndent(priv, "", "  ")
	return os.WriteFile(path, out, 0644)
}

func saveSignature(path string, message, x0, x1, target []uint64, signature [][]uint64) error {
	_ = os.MkdirAll("Signature", 0755)
	sig := SignatureData{
		Message:   message,
		X0:        x0,
		X1:        x1,
		Target:    target,
		Signature: signature,
	}
	data, _ := json.MarshalIndent(sig, "", "  ")
	return os.WriteFile(path, data, 0644)
}
