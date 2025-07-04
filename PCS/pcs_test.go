package pcs_test

import (
	"math/rand"
	"testing"
	"time"

	pcs "vSIS-Signature/PCS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

func TestPCSEndToEnd(t *testing.T) {
	// 1) Set up the NTT ring (degree N=1024, a prime ≡1 mod 2N)
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, err := ring.NewRing(N, moduli)
	if err != nil {
		t.Fatalf("failed to create ring: %v", err)
	}

	// 2) PCS parameters: DECS repetition η, block size μ, mask‐rows ℓ′
	const (
		Eta      = 4
		Mu       = 3
		EllPrime = 5
	)
	// Degrees of the two test polynomials
	Dj := []int{7, 11}
	npcs := len(Dj)

	// 3) Random polynomials P_j of degree Dj[j]
	rand.Seed(time.Now().UnixNano())
	P := make([][]uint64, npcs)
	for j := 0; j < npcs; j++ {
		P[j] = make([]uint64, Dj[j]+1)
		for i := range P[j] {
			P[j][i] = rand.Uint64() % ringQ.Modulus[0]
		}
	}

	// 4) Prover’s first commit step
	com, pk, err := pcs.CommitInitPolynomial(ringQ, P, Dj, Mu, EllPrime)
	if err != nil {
		t.Fatalf("CommitInitPolynomial error: %v", err)
	}

	// 5) Verifier setup & CommitStep1
	vk := pcs.NewVerifierPolynomial(ringQ, Dj, Mu, EllPrime, Eta)
	vk.CommitStep1Polynomial(com)

	// 6) Prover finishes commit by sending R_k(X)
	Rpolys, err := pcs.CommitFinishPolynomial(pk)
	if err != nil {
		t.Fatalf("CommitFinishPolynomial error: %v", err)
	}
	vk.CommitStep2Polynomial(Rpolys)

	// 7) Prepare a random linear‐map C[k][j]
	m := 3
	nrows := Mu + EllPrime
	C := make([][]uint64, m)
	for k := 0; k < m; k++ {
		C[k] = make([]uint64, nrows)
		for j := 0; j < nrows; j++ {
			C[k][j] = rand.Uint64() % ringQ.Modulus[0]
		}
	}

	// 8) Prover masks those rows: bar[k][i] = Σ_j C[k][j]*mask_j[i]
	bar, err := pcs.EvalInitPolynomial(pk, C)
	if err != nil {
		t.Fatalf("EvalInitPolynomial error: %v", err)
	}

	// 9) Verifier chooses ℓ′ random subset E
	E := vk.ChooseEPolynomial()

	// 10) Prover opens via DECS at indices E
	open, err := pcs.EvalFinishPolynomial(pk, E)
	if err != nil {
		t.Fatalf("EvalFinishPolynomial error: %v", err)
	}

	// 11) Verifier checks Merkle, degree, and final polynomial relation
	if ok := vk.EvalStep2Polynomial(bar, E, open, C); !ok {
		t.Fatalf("EvalStep2Polynomial failed")
	}
}
