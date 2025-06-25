// ProverZKP/gate_builder.go
package proverzkp

import (
	"encoding/json"
	"fmt"
	"log"
	"os"
	qgate "vSIS-Signature/Quadratic_Gate"
	signer "vSIS-Signature/Signer"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// -----------------------------------------------------------------------------
// Tiny helpers for loading what the Signer already produced
// -----------------------------------------------------------------------------

type signatureJSON struct {
	Message   []uint64   `json:"message"`
	X0        []uint64   `json:"x0"`
	X1        []uint64   `json:"x1"`
	Target    []uint64   `json:"target"`
	Signature [][]uint64 `json:"signature"`
}

func loadSignature(path string) (*signatureJSON, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var s signatureJSON
	return &s, json.Unmarshal(raw, &s)
}

func loadBmatrixCoeffs(path string) ([][]uint64, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var tmp struct {
		B [][]uint64 `json:"B"`
	}
	return tmp.B, json.Unmarshal(raw, &tmp)
}

// --------------------------------------------------------------------------
//
//	BuildGateFromDisk  –  corrected version
//
// --------------------------------------------------------------------------
func BuildGateFromDisk() (qgate.GatePublic, qgate.GatePrivate, error) {

	// ‣ 0. parameters ----------------------------------------------------------
	par, err := signer.LoadParams("Parameters/Parameters.json")
	if err != nil {
		return qgate.GatePublic{}, qgate.GatePrivate{}, err
	}

	ringQ, err := ring.NewRing(par.N, []uint64{par.Q})
	if err != nil {
		return qgate.GatePublic{}, qgate.GatePrivate{}, err
	}

	// convenience: explicit in-place lift
	toNTT := func(p *ring.Poly) { ringQ.NTT(p, p) }

	// ‣ 1. matrix A  (already stored in NTT, **do not** touch) -----------------
	pk, err := signer.LoadPublicKey("public_key/public_key.json")
	if err != nil {
		return qgate.GatePublic{}, qgate.GatePrivate{}, err
	}

	A := make([][]*ring.Poly, 1)
	A[0] = make([]*ring.Poly, len(pk.A))
	for i, coeffs := range pk.A {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs) // pk.A is already NTT
		A[0][i] = p               // no toNTT
	}

	// ‣ 2. B-matrix columns  (stored in coefficient domain → lift) -------------
	Bcoeffs, err := loadBmatrixCoeffs("Parameters/Bmatrix.json")
	if err != nil {
		return qgate.GatePublic{}, qgate.GatePrivate{}, err
	}
	if len(Bcoeffs) != 4 {
		return qgate.GatePublic{}, qgate.GatePrivate{}, fmt.Errorf("expected 4 polys in B")
	}

	B0Const := make([]*ring.Poly, 1)
	B0Msg := [][]*ring.Poly{make([]*ring.Poly, 1)}
	B0Rnd := [][]*ring.Poly{make([]*ring.Poly, 1)}
	b1 := make([]*ring.Poly, 1)

	makePolyNTT := func(raw []uint64) *ring.Poly {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], raw)
		toNTT(p) // lift to NTT **once**
		return p
	}

	B0Const[0] = makePolyNTT(Bcoeffs[0])  //   B₀,const
	B0Msg[0][0] = makePolyNTT(Bcoeffs[1]) //   B₀,msg
	B0Rnd[0][0] = makePolyNTT(Bcoeffs[2]) //   B₀,rnd
	b1[0] = makePolyNTT(Bcoeffs[3])       //   b₁

	// ‣ 3. ρ  (compression vector) --------------------------------------------
	prng, _ := utils.NewPRNG()
	rho := make([]uint64, 1)
	rbuf := make([]byte, 8)
	prng.Read(rbuf)
	rho[0] = uint64(rbuf[0]) % par.Q // small is fine

	// ‣ 4. signature blobs -----------------------------------------------------
	sig, err := loadSignature("Signature/Signature.json")
	if err != nil {
		return qgate.GatePublic{}, qgate.GatePrivate{}, err
	}

	// message  u   and masks  x₀, x₁  – stored in coeff domain → lift
	m := makePolyNTT(sig.Message)
	x0 := makePolyNTT(sig.X0)
	x1 := makePolyNTT(sig.X1)

	// signature vector s  – already NTT (output of GaussSamp) ------------------
	s := make([]*ring.Poly, len(sig.Signature))
	for i, coeffs := range sig.Signature {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs) // keep as is (NTT)
		s[i] = p
	}

	// ‣ 5. build quadratic gate -----------------------------------------------
	pub, priv, err := qgate.BuildGate(
		ringQ,
		A, b1,
		B0Const, B0Msg, B0Rnd,
		rho,
		/*private*/ s, x1, []*ring.Poly{m}, []*ring.Poly{x0})

	if err != nil {
		log.Fatal(err)
	}

	if err := qgate.VerifyGate(ringQ, pub, priv); err != nil {
		log.Fatalf("gate verification: %v", err)
	}
	// fmt.Println("✓   xᵀR₂x + r₁ᵀx + r₀ = 0  holds")
	// fmt.Print("public key: ")
	// fmt.Printf("R0 : %v\n", pub.R0.Coeffs)
	return pub, priv, nil
}
