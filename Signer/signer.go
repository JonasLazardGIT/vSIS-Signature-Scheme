package signer

import (
	"encoding/json"
	"fmt"
	"log"
	"os"
	"time"

	ps "vSIS-Signature/Preimage_Sampler"
	Parameters "vSIS-Signature/System"
	prof "vSIS-Signature/prof"
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
	defer prof.Track(time.Now(), "Sign")
	// 1) load public parameters
	params, err := LoadParams("Parameters/Parameters.json")
	if err != nil {
		log.Fatalf("loadParams: %v", err)
	}
	log.Printf("Params: N=%d Q=%d Base=%d K=%d", params.N, params.Q, params.Base, params.K)

	// 2) instantiate ring R_q
	ringQ, err := ring.NewRing(params.N, []uint64{params.Q})
	if err != nil {
		log.Fatalf("ring.NewRing: %v", err)
	}

	var trap ps.Trapdoor
	if fileExists("public_key/public_key.json") && fileExists("private_key/private_key.json") {
		log.Println("Loading existing keypair...")
		pk, err := LoadPublicKey("public_key/public_key.json")
		if err != nil {
			log.Fatalf("loadPublicKey: %v", err)
		}
		sk, err := loadPrivateKey("private_key/private_key.json")
		if err != nil {
			log.Fatalf("loadPrivateKey: %v", err)
		}
		trap = assembleTrapdoor(ringQ, pk, sk)
	} else {
		log.Println("Generating new keypair...")
		trap = ps.TrapGen(ringQ, params.Base, params.SigmaT)
		if err := savePublicKey("public_key/public_key.json", ringQ, &trap); err != nil {
			log.Fatalf("savePublicKey: %v", err)
		}
		if err := savePrivateKey("private_key/private_key.json", ringQ, &trap); err != nil {
			log.Fatalf("savePrivateKey: %v", err)
		}
	}

	// 5) load B‐matrix from JSON and convert to CyclotomicFieldElem
	Bcyclo, err := loadBMatrix("Parameters/Bmatrix.json", ringQ)
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
	tNTT, err := vsishash.ComputeBBSHash(ringQ, Bcyclo, mPoly, x0Poly, x1Poly)
	if err != nil {
		log.Fatalf("ComputeBBSHash: %v", err)
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

func LoadParams(path string) (*Parameters.SystemParams, error) {
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

func savePublicKey(path string, ringQ *ring.Ring, trap *ps.Trapdoor) error {
	_ = os.MkdirAll("public_key", 0755)
	pub := PublicKey{Base: trap.Base, K: trap.K, A: make([][]uint64, len(trap.A))}
	for i, a := range trap.A {
		pub.A[i] = append([]uint64(nil), a.Coeffs[0]...)
	}
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	if _, err := f.WriteString("{\n  \"A\": [\n"); err != nil {
		return err
	}
	for i, poly := range pub.A {
		line, _ := json.Marshal(poly)
		if i < len(pub.A)-1 {
			fmt.Fprintf(f, "    %s,\n", line)
		} else {
			fmt.Fprintf(f, "    %s\n", line)
		}
	}
	if _, err := f.WriteString("  ],\n"); err != nil {
		return err
	}
	fmt.Fprintf(f, "  \"base\": %d,\n", pub.Base)
	fmt.Fprintf(f, "  \"k\": %d\n", pub.K)
	_, err = f.WriteString("}\n")
	return err
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
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	if _, err := f.WriteString("{\n  \"R0\": [\n"); err != nil {
		return err
	}
	for i, poly := range priv.R0 {
		line, _ := json.Marshal(poly)
		if i < len(priv.R0)-1 {
			fmt.Fprintf(f, "    %s,\n", line)
		} else {
			fmt.Fprintf(f, "    %s\n", line)
		}
	}
	if _, err := f.WriteString("  ],\n  \"R1\": [\n"); err != nil {
		return err
	}
	for i, poly := range priv.R1 {
		line, _ := json.Marshal(poly)
		if i < len(priv.R1)-1 {
			fmt.Fprintf(f, "    %s,\n", line)
		} else {
			fmt.Fprintf(f, "    %s\n", line)
		}
	}
	if _, err := f.WriteString("  ],\n"); err != nil {
		return err
	}
	fmt.Fprintf(f, "  \"base\": %d,\n", priv.Base)
	fmt.Fprintf(f, "  \"k\": %d\n", priv.K)
	_, err = f.WriteString("}\n")
	return err
}

func saveSignature(path string, message, x0, x1, target []uint64, signature [][]uint64) error {
	_ = os.MkdirAll("Signature", 0755)
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	msg, _ := json.Marshal(message)
	x0b, _ := json.Marshal(x0)
	x1b, _ := json.Marshal(x1)
	tgt, _ := json.Marshal(target)
	if _, err := fmt.Fprintf(f, "{\n  \"message\": %s,\n  \"x0\": %s,\n  \"x1\": %s,\n  \"target\": %s,\n  \"signature\": [\n", msg, x0b, x1b, tgt); err != nil {
		return err
	}
	for i, poly := range signature {
		line, _ := json.Marshal(poly)
		if i < len(signature)-1 {
			fmt.Fprintf(f, "    %s,\n", line)
		} else {
			fmt.Fprintf(f, "    %s\n", line)
		}
	}
	_, err = f.WriteString("  ]\n}\n")
	return err
}

func LoadPublicKey(path string) (*PublicKey, error) {
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

func loadPrivateKey(path string) (*PrivateKey, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var sk PrivateKey
	if err := json.Unmarshal(data, &sk); err != nil {
		return nil, err
	}
	return &sk, nil
}

func assembleTrapdoor(ringQ *ring.Ring, pk *PublicKey, sk *PrivateKey) ps.Trapdoor {
	trap := ps.Trapdoor{Base: pk.Base, K: pk.K, Rows: 1, Cols: 2}
	trap.A = make([]*ring.Poly, len(pk.A))
	for i, a := range pk.A {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], a)
		trap.A[i] = p
	}
	trap.R[0] = make([]*ring.Poly, len(sk.R0))
	for i, r0 := range sk.R0 {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], r0)
		trap.R[0][i] = p
	}
	trap.R[1] = make([]*ring.Poly, len(sk.R1))
	for i, r1 := range sk.R1 {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], r1)
		trap.R[1][i] = p
	}
	return trap
}

func fileExists(path string) bool {
	if _, err := os.Stat(path); err == nil {
		return true
	}
	return false
}
