package decs

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func testParams(ringQ *ring.Ring, eta int, ell int) Params {
	return Params{Degree: int(ringQ.N - 1), Eta: eta, NonceBytes: 16}
}

func TestDECS_CommitEval_Accepts(t *testing.T) {
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, err := ring.NewRing(N, moduli)
	if err != nil {
		t.Fatal(err)
	}
	eta, ell := 2, 64
	params := testParams(ringQ, eta, ell)

	r := 5
	Ps := make([]*ring.Poly, r)
	prng, _ := utils.NewPRNG()
	us := ring.NewUniformSampler(prng, ringQ)
	for j := 0; j < r; j++ {
		Ps[j] = ringQ.NewPoly()
		us.Read(Ps[j])
		for i := params.Degree + 1; i < int(N); i++ {
			Ps[j].Coeffs[0][i] = 0
		}
	}
	prover := NewProverWithParams(ringQ, Ps, params)
	root, err := prover.CommitInit()
	if err != nil {
		t.Fatal(err)
	}
	verifier := NewVerifierWithParams(ringQ, r, params)
	Gamma := verifier.DeriveGamma(root)
	R := prover.CommitStep2(Gamma)
	if !verifier.VerifyCommit(root, R, Gamma) {
		t.Fatal("VerifyCommit failed (should accept)")
	}

	E := make([]int, ell)
	for i := range E {
		x, _ := rand.Int(rand.Reader, big.NewInt(int64(N)))
		E[i] = int(x.Uint64())
	}
	open := prover.EvalOpen(E)
	if !verifier.VerifyEvalAt(root, Gamma, R, open, E) {
		t.Fatal("VerifyEvalAt failed (should accept)")
	}
}

func TestDECS_Rejects_HighDegree(t *testing.T) {
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, _ := ring.NewRing(N, moduli)
       eta := 2
       params := Params{Degree: int(ringQ.N - 2), Eta: eta, NonceBytes: 16}

	r := 3
	Ps := make([]*ring.Poly, r)
	for j := 0; j < r; j++ {
		Ps[j] = ringQ.NewPoly()
	}
	Ps[0].Coeffs[0][params.Degree+1] = 1

	prover := NewProverWithParams(ringQ, Ps, params)
	root, err := prover.CommitInit()
	if err != nil {
		t.Fatal(err)
	}
	verifier := NewVerifierWithParams(ringQ, r, params)
	Gamma := verifier.DeriveGamma(root)
	R := prover.CommitStep2(Gamma)

	if verifier.VerifyCommit(root, R, Gamma) {
		t.Fatal("VerifyCommit accepted R with degree > d")
	}
}
