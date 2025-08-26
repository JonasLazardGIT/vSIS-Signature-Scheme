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
	ncols  int

	Root  [32]byte
	Gamma [][]uint64
	R     []*ring.Poly
}

// NewVerifier constructs the LVCS verifier.
func NewVerifier(ringQ *ring.Ring, r, eta, ncols int) *VerifierState {
	return &VerifierState{RingQ: ringQ, r: r, eta: eta, ncols: ncols}
}

// CommitStep1 – §4.1 steps 1–3:
// Commit all those polynomials via DECS and record the commitment root.
func (v *VerifierState) CommitStep1(root [32]byte) [][]uint64 {
	v.Root = root
	decv := decs.NewVerifier(v.RingQ, v.r, v.eta)
	v.Gamma = decv.DeriveGamma(root)
	return v.Gamma
}

// CommitStep2 stores the prover's R_k polynomials.
func (v *VerifierState) CommitStep2(R []*ring.Poly) bool {
	v.R = R
	for _, p := range R {
		// quick degree bound using coeff-form length
		if deg(p) > decs.Degree {
			return false
		}
	}
	return true
}

// tiny local helper
func deg(p *ring.Poly) int {
	c := p.Coeffs[0]
	for i := len(c) - 1; i >= 0; i-- {
		if c[i] != 0 {
			return i
		}
	}
	return 0
}

// ChooseE – choose ℓ distinct indices on the MASKED TAIL [ncols, ncols+ℓ).
// Pass ncols := N - ℓ so E ⊆ Ω′ (the blinded coordinates), per §4.1.
func (v *VerifierState) ChooseE(ell, ncols int) []int {
	N := v.RingQ.N
	if ell <= 0 || ncols < 0 || ncols+ell > N {
		return nil
	}
	E := make([]int, 0, ell)
	rnd := rand.New(rand.NewSource(time.Now().UnixNano()))
	used := make(map[int]struct{}, ell)
	for len(E) < ell {
		idx := ncols + rnd.Intn(ell)
		if _, ok := used[idx]; ok {
			continue
		}
		used[idx] = struct{}{}
		E = append(E, idx)
	}
	return E
}

// *simplified* verifier: only Merkle + degree test (bar / C not needed
// by the current one-shot simulation).
func (v *VerifierState) EvalStep2Slim(open *decs.DECSOpening) bool {
	decv := decs.NewVerifier(v.RingQ, v.r, v.eta)
	return decv.VerifyEval(v.Root, v.Gamma, v.R, open)
}

// EvalStep2 – §4.1 step 4:
// Verify Merkle + low-degree + linear checks, leaking only the ℓ masked positions.
func (v *VerifierState) EvalStep2(
	bar [][]uint64, // prover’s ¯v_k
	E []int, // challenge set
	open *decs.DECSOpening,
	C [][]uint64, // coefficient matrix
) bool {
	// enforce tail-only and exact cardinality
	ncols := v.ncols
	if len(E) != len(bar[0]) {
		return false
	}
	ell := len(bar[0])
	for _, idx := range E {
		if idx < ncols || idx >= ncols+ell {
			return false
		}
	}

	// 0) Bind the openings to the challenge set E (set equality)
	if !equalSets(open.Indices, E) {
		return false
	}

	// 1) Merkle + masked-relation check (now bound to E)
	decv := decs.NewVerifier(v.RingQ, v.r, v.eta)
	if !decv.VerifyEvalAt(v.Root, v.Gamma, v.R, open, E) {
		return false
	}

	// 2) check only the masked positions
	mod := v.RingQ.Modulus[0]
	for t, idx := range open.Indices {
		if idx < ncols || idx >= ncols+ell {
			return false
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

// equalSets checks multisets equality of int slices.
func equalSets(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	seen := make(map[int]int, len(a))
	for _, x := range a {
		seen[x]++
	}
	for _, y := range b {
		if seen[y] == 0 {
			return false
		}
		seen[y]--
	}
	return true
}
