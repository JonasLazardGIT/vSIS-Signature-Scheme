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

// CommitStep1 – §4.1 steps 1–3:
// Commit all those polynomials via DECS and record the commitment root.
func (v *VerifierState) CommitStep1(root [32]byte) {
	v.Root = root
	decv := decs.NewVerifier(v.RingQ, v.r, v.eta)
	v.Gamma = decv.DeriveGamma(root) // <- fixed: use DeriveGamma on Verifier, not decs.DeriveGamma
}

// CommitStep2 stores the prover's R_k polynomials.
func (v *VerifierState) CommitStep2(R []*ring.Poly) {
	v.R = R
}

// ChooseE – §4.1 step 3:
// Later, open masked linear combinations on a small random subset EE of size ℓ.
func (v *VerifierState) ChooseE(ell int) []int {
	N := v.RingQ.N
	E := make([]int, ell)
	rnd := rand.New(rand.NewSource(time.Now().UnixNano()))
	for i := range E {
		E[i] = rnd.Intn(N)
	}
	return E
}

// EvalStep2 – §4.1 step 4:
// Verify Merkle + low-degree + linear checks, leaking only the ℓ masked positions.
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

	// 2) check only the masked positions
	mod := v.RingQ.Modulus[0]
	ncols := v.RingQ.N - len(bar[0])
	for t, idx := range open.Indices {
		if idx < ncols {
			continue
		}
		maskedPos := idx - ncols
		for k := range bar {
			acc := uint64(0)
			for j := 0; j < v.r; j++ {
				acc = (acc + C[k][j]*open.Pvals[t][j]) % mod
			}
			if acc != bar[k][maskedPos] {
				return false
			}
		}
	}
	return true
}
