package lvcs

import (
	"fmt"
	"math/rand"
	"os"
	"time"

	decs "vSIS-Signature/DECS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

var debugLVCS = os.Getenv("DEBUG_LVCS") != ""

// VerifierState holds verifier‐side LVCS state.
type VerifierState struct {
	RingQ  *ring.Ring
	r      int
	params decs.Params

	Root  [32]byte
	Gamma [][]uint64
	R     []*ring.Poly
}

// NewVerifierWithParams constructs the LVCS verifier with explicit DECS params.
func NewVerifierWithParams(ringQ *ring.Ring, r int, params decs.Params) *VerifierState {
	return &VerifierState{RingQ: ringQ, r: r, params: params}
}

// NewVerifier is a backwards-compatible helper taking η and ncols. ncols is
// ignored (EvalStep2 derives it from bar), but kept for API compatibility.
func NewVerifier(ringQ *ring.Ring, r, eta, ncols int) *VerifierState {
	params := decs.DefaultParams
	params.Eta = eta
	if params.Degree >= int(ringQ.N) {
		params.Degree = int(ringQ.N) - 1
	}
	return NewVerifierWithParams(ringQ, r, params)
}

// CommitStep1 – §4.1 steps 1–3:
// Commit all those polynomials via DECS and record the commitment root.
func (v *VerifierState) CommitStep1(root [32]byte) [][]uint64 {
	v.Root = root
	decv := decs.NewVerifierWithParams(v.RingQ, v.r, v.params)
	v.Gamma = decv.DeriveGamma(root)
	return v.Gamma
}

// CommitStep2 stores the prover's R_k polynomials.
func (v *VerifierState) CommitStep2(R []*ring.Poly) bool {
	v.R = R
	for _, p := range R {
		// quick degree bound using coeff-form length
		if deg(p) > v.params.Degree {
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
	decv := decs.NewVerifierWithParams(v.RingQ, v.r, v.params)
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
	if len(bar) == 0 || len(bar[0]) == 0 {
		if debugLVCS {
			fmt.Println("[LVCS] FAIL: empty bar")
		}
		return false
	}
	ell := len(bar[0])
	ncols := int(v.RingQ.N) - ell
	if len(E) != ell {
		if debugLVCS {
			fmt.Printf("[LVCS] FAIL: |E|=%d != ell=%d\n", len(E), ell)
		}
		return false
	}
	for _, idx := range E {
		if idx < ncols || idx >= ncols+ell {
			if debugLVCS {
				fmt.Printf("[LVCS] FAIL: E contains head index %d (ncols=%d ell=%d)\n", idx, ncols, ell)
			}
			return false
		}
	}

	// 0) Bind the openings to the challenge set E (set equality)
	if !equalSets(open.Indices, E) {
		if debugLVCS {
			fmt.Println("[LVCS] FAIL: open.Indices != E")
		}
		return false
	}

	// 1) Merkle + masked-relation check (now bound to E)
	decv := decs.NewVerifierWithParams(v.RingQ, v.r, v.params)
	if !decv.VerifyEvalAt(v.Root, v.Gamma, v.R, open, E) {
		if debugLVCS {
			fmt.Println("[LVCS] FAIL: DECS.VerifyEvalAt")
		}
		return false
	}

	// 2) check only the masked positions
	mod := v.RingQ.Modulus[0]
	for t, idx := range open.Indices {
		maskedPos := idx - ncols
		for k := range bar {
			acc := uint64(0)
			for j := 0; j < v.r; j++ {
				acc = (acc + C[k][j]*open.Pvals[t][j]) % mod
			}
			if acc != bar[k][maskedPos] {
				if debugLVCS {
					fmt.Printf("[LVCS] FAIL linear at t=%d idx=%d (pos=%d) k=%d: acc=%d bar=%d\n", t, idx, maskedPos, k, acc, bar[k][maskedPos])
					fmt.Printf("       C[k]=%v\n", C[k])
					fmt.Printf("       Pvals[t]=%v\n", open.Pvals[t])
				}
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
