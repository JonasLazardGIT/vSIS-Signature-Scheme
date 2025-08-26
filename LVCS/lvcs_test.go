// lvcs_test.go
package lvcs

import (
	"crypto/rand"
	"math/big"
	"testing"

	decs "vSIS-Signature/DECS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

func TestLVCSCommitAndEval(t *testing.T) {
	// ───────────────────────────────────────────────────────────
	// 1) SETUP
	// Choose a small NTT ring. In real use you’d pick a secure
	// parameter set; here we use toy values so the test runs quickly.
	// N must be a power of two and q ≡ 1 mod 2N.
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, err := ring.NewRing(N, moduli)
	if err != nil {
		t.Fatalf("ring.NewRing: %v", err)
	}

	// Protocol parameters
	nrows := 4 // number of row-vectors r_j
	ell := 8   // masking‐size ℓ

	// Generate random row-vectors r_j ∈ F_q^N
	rows := make([][]uint64, nrows)
	q0 := ringQ.Modulus[0]
	for j := 0; j < nrows; j++ {
		rows[j] = make([]uint64, ringQ.N-ell)
		for i := 0; i < ringQ.N-ell; i++ {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			rows[j][i] = x.Uint64()
		}
	}

	// ───────────────────────────────────────────────────────────
	// 2) LVCS.Commit
	// Prover: init
	root, proverKey, err := CommitInit(ringQ, rows, ell)
	if err != nil {
		t.Fatalf("CommitInit: %v", err)
	}

	// Verifier: record commitment & sample Γ
	ver := NewVerifier(ringQ, nrows, decs.Eta, ringQ.N-ell)
	ver.CommitStep1(root)

	// Prover: finish (recv Γ, send R)
	R := CommitFinish(proverKey, ver.Gamma)

	// Verifier: receive R
	ver.CommitStep2(R)

	// ───────────────────────────────────────────────────────────
	// 3) LVCS.Eval
	// Both sides agree on a coefficient matrix C of size m×nrows
	m := 3
	C := make([][]uint64, m)
	for k := 0; k < m; k++ {
		C[k] = make([]uint64, nrows)
		for j := 0; j < nrows; j++ {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			C[k][j] = x.Uint64()
		}
	}

	// Prover: step1 → send bar = EvalInit(...)
	bar := EvalInit(ringQ, proverKey, C)

	// Verifier: choose E ⊆ Ω′ (masked tail)
	ncols := ringQ.N - ell
	E := ver.ChooseE(ell, ncols)

	// Prover: step3+4 → send opening = EvalFinish(...)
	opening := EvalFinish(proverKey, E)

	// Verifier: final check
	ok := ver.EvalStep2(bar, E, opening.DECSOpen, C)
	if !ok {
		t.Fatal("LVCS.EvalStep2 failed: proof rejected")
	}

	t.Log("LVCS commit+eval succeeded")
}

// Negative tests for LVCS EvalStep2
func TestLVCSRejectsBadOpenings(t *testing.T) {
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, err := ring.NewRing(N, moduli)
	if err != nil {
		t.Fatalf("ring.NewRing: %v", err)
	}
	nrows := 3
	ell := 4
	q0 := ringQ.Modulus[0]
	rows := make([][]uint64, nrows)
	for j := 0; j < nrows; j++ {
		rows[j] = make([]uint64, ringQ.N-ell)
		for i := range rows[j] {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			rows[j][i] = x.Uint64()
		}
	}
	root, proverKey, err := CommitInit(ringQ, rows, ell)
	if err != nil {
		t.Fatalf("CommitInit: %v", err)
	}
	ver := NewVerifier(ringQ, nrows, decs.Eta, ringQ.N-ell)
	ver.CommitStep1(root)
	R := CommitFinish(proverKey, ver.Gamma)
	ver.CommitStep2(R)
	m := 2
	C := make([][]uint64, m)
	for k := 0; k < m; k++ {
		C[k] = make([]uint64, nrows)
		for j := 0; j < nrows; j++ {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			C[k][j] = x.Uint64()
		}
	}
	bar := EvalInit(ringQ, proverKey, C)
	ncols := ringQ.N - ell

	// Case 1: prover opens at E but verifier expects Ebad
	E := ver.ChooseE(ell, ncols)
	opening := EvalFinish(proverKey, E)
	Ebad := make([]int, len(E))
	copy(Ebad, E)
	if Ebad[0] == ncols {
		Ebad[0]++
	} else {
		Ebad[0] = ncols
	}
	if ver.EvalStep2(bar, Ebad, opening.DECSOpen, C) {
		t.Fatal("EvalStep2 accepted opening for mismatched E")
	}

	// Case 2: verifier’s E includes a head index
	Ehead := make([]int, len(E))
	copy(Ehead, E)
	Ehead[0] = 0 // head position
	openHead := EvalFinish(proverKey, Ehead)
	if ver.EvalStep2(bar, Ehead, openHead.DECSOpen, C) {
		t.Fatal("EvalStep2 accepted opening with head index")
	}
}
