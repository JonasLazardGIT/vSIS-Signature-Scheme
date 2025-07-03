package lvcs

import (
	"math/rand"
	"time"

	decs "vSIS-Signature/DECS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// VerifierState holds verifier‐side LVCS state.
type VerifierState struct {
	RingQ  *ring.Ring
	r, eta int

	Root  [32]byte
	Gamma [][]uint64
	R     []*ring.Poly
}

// NewVerifier constructs the LVCS verifier.
func NewVerifier(ringQ *ring.Ring, r, eta int) *VerifierState {
	return &VerifierState{RingQ: ringQ, r: r, eta: eta}
}

// CommitStep1: record the root, sample Γ via the exported method.
func (v *VerifierState) CommitStep1(root [32]byte) {
	v.Root = root
	decv := decs.NewVerifier(v.RingQ, v.r, v.eta)
	v.Gamma = decv.DeriveGamma(root) // <- fixed: use DeriveGamma on Verifier, not decs.DeriveGamma
}

// CommitStep2: record the prover’s R_k polys
func (v *VerifierState) CommitStep2(R []*ring.Poly) {
	v.R = R
}

// ChooseE picks a uniformly random size‐ℓ subset of [0..N).
func (v *VerifierState) ChooseE(ell int) []int {
	N := v.RingQ.N
	E := make([]int, ell)
	rnd := rand.New(rand.NewSource(time.Now().UnixNano()))
	for i := range E {
		E[i] = rnd.Intn(N)
	}
	return E
}

// EvalStep2 does the final LVCS.Eval verification:
//
//	(1) delegate Merkle+DECS.Eval checks to decs.Verifier.VerifyEval
//	(2) check Qk(e) = Σ_j C[k][j]·Pj(e) for each k, t
func (v *VerifierState) EvalStep2(
	bar [][]uint64, // prover’s ¯v_k
	E []int, // challenge set
	open *decs.DECSOpening,
	C [][]uint64, // coefficient matrix
) bool {
	// 1) Merkle + mask‐relation check
	decv := decs.NewVerifier(v.RingQ, v.r, v.eta)
	if !decv.VerifyEval(v.Root, v.Gamma, v.R, open) {
		return false
	}

	// 2) check for each k,t: Q_k(e_t) ?= Σ_j C[k][j]·P_j(e_t)
	mod := v.RingQ.Modulus[0]
	for k := range bar {
		for t := range open.Pvals {
			acc := uint64(0)
			for j := 0; j < v.r; j++ {
				acc = (acc + C[k][j]*open.Pvals[t][j]) % mod
			}
			if acc != bar[k][t] {
				return false
			}
		}
	}
	return true
}
