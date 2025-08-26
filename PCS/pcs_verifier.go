// Package pcs implements the SmallWood polynomial‐commitment wrapper
// on top of LVCS (and DECS).  This file contains the verifier‐side routines.
package pcs

import (
	lvcs "vSIS-Signature/LVCS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// VerifierState holds all the verifier’s state for a PCS instance.
type VerifierState struct {
	RingQ    *ring.Ring // NTT context
	LvcsVer  *lvcs.VerifierState
	Dj       []int // degrees of each P_j
	Mu       int   // block size μ
	EllPrime int   // mask‐row count ℓ′
	Nu       []int // column counts ν_j
	Nrows    int   // total row‐count = Mu+EllPrime
}

// NewVerifierPolynomial instantiates a PCS verifier for polynomials of degrees Dj,
// using stripe‐size Mu, mask size EllPrime, and DECS repetition Eta.
func NewVerifierPolynomial(
	ringQ *ring.Ring,
	Dj []int,
	Mu, EllPrime, Eta int,
) *VerifierState {
	npcs := len(Dj)
	Nu := make([]int, npcs)
	sumNu := 0
	for j := 0; j < npcs; j++ {
		d := Dj[j] + 1 - EllPrime
		if d <= 0 {
			Nu[j] = 1
		} else {
			Nu[j] = (d + Mu - 1) / Mu
		}
		sumNu += Nu[j]
	}
	nrows := Mu + EllPrime
	// Initialize the LVCS verifier with the same row‐count, DECS η, and ncols
	lvcsVer := lvcs.NewVerifier(ringQ, nrows, Eta, sumNu)
	return &VerifierState{
		RingQ:    ringQ,
		LvcsVer:  lvcsVer,
		Dj:       Dj,
		Mu:       Mu,
		EllPrime: EllPrime,
		Nu:       Nu,
		Nrows:    nrows,
	}
}

// CommitStep1Polynomial ingests the PCS commitment root and derives the LVCS challenge Gamma.
func (vs *VerifierState) CommitStep1Polynomial(com Commitment) {
	vs.LvcsVer.CommitStep1(com.Root)
}

// CommitStep2Polynomial ingests the prover’s R_k(X) polynomials.
func (vs *VerifierState) CommitStep2Polynomial(Rpolys []*ring.Poly) {
	vs.LvcsVer.CommitStep2(Rpolys)
}

// ChooseEPolynomial picks a random ℓ′-subset on the masked tail for opening.
func (vs *VerifierState) ChooseEPolynomial() []int {
	sumNu := 0
	for _, nu := range vs.Nu {
		sumNu += nu
	}
	return vs.LvcsVer.ChooseE(vs.EllPrime, sumNu)
}

// EvalStep2Polynomial completes the PCS open protocol:
// 1) runs LVCS.EvalStep2 to check Merkle, degree, and per‐row linear constraints,
// 2) returns true iff all checks pass.
func (vs *VerifierState) EvalStep2Polynomial(
	bar [][]uint64, // masked targets bar[k][i]
	E []int, // challenge indices
	open *lvcs.Opening, // LVCS opening
	C [][]uint64, // coefficient matrix C[k][j]
) bool {
	// 1) LVCS‐level verification
	if !vs.LvcsVer.EvalStep2(bar, E, open.DECSOpen, C) {
		return false
	}
	// 2) (No extra work needed: LVCS.EvalStep2 already enforces
	//    Q_k(e) = Σ_j C[k][j]*P_j(e) on the opened subset.)
	return true
}
