package decs

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func TestDECSProtocol(t *testing.T) {
	// 1) Setup ring
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, err := ring.NewRing(N, moduli)
	if err != nil {
		t.Fatal(err)
	}

	// 2) sample random P_j of deg ≤ Degree
	r := 5
	Ps := make([]*ring.Poly, r)
	prng, _ := utils.NewPRNG()
	us := ring.NewUniformSampler(prng, ringQ)
	for j := 0; j < r; j++ {
		Ps[j] = ringQ.NewPoly()
		us.Read(Ps[j])
		for i := Degree + 1; i < N; i++ {
			Ps[j].Coeffs[0][i] = 0
		}
	}

	// 3) prover.CommitInit → root
	prover := NewProver(ringQ, Ps)
	root, err := prover.CommitInit()
	if err != nil {
		t.Fatal(err)
	}

	// 4) verifier → Gamma
	verifier := NewVerifier(ringQ, r, Eta)
	Gamma := verifier.DeriveGamma(root)

	// 5) prover → R
	R := prover.CommitStep2(Gamma)

	// 6) verifier.VerifyCommit
	if !verifier.VerifyCommit(root, R, Gamma) {
		t.Fatal("VerifyCommit failed")
	}

	// 7) pick random E of size ℓ
	ℓ := 64
	E := make([]int, ℓ)
	for i := 0; i < ℓ; i++ {
		idxBig, _ := rand.Int(rand.Reader, big.NewInt(int64(N)))
		E[i] = int(idxBig.Int64())
	}

	// 8) prover opens
	open := prover.EvalOpen(E)

	// 9) verifier.VerifyEval
	if !verifier.VerifyEval(root, Gamma, R, open) {
		t.Error("VerifyEval failed")
	}
}
