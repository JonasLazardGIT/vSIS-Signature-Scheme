// lvcs_test.go
package lvcs

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	decs "vSIS-Signature/DECS"
)

func TestLVCSCommitAndEval(t *testing.T) {
	// ───────────────────────────────────────────────────────────
	// 1) SETUP
	// Choose a small NTT ring. In real use you’d pick a secure
	// parameter set; here we use toy values so the test runs quickly.
	// N must be a power of two and q ≡ 1 mod 2N.
	N := 1 << 11
	// 0x1000000000001 is prime and ≡1 mod 32
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
	ver := NewVerifier(ringQ, nrows, decs.Eta)
	ver.CommitStep1(root)

	// Prover: finish (recv Γ, send R)
	R, err := CommitFinish(proverKey, ver.Gamma)
	if err != nil {
		t.Fatalf("CommitFinish: %v", err)
	}

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

	// Verifier: choose E uniformly of size ℓ
	E := ver.ChooseE(ell)

	// Prover: step3+4 → send opening = EvalFinish(...)
	opening := EvalFinish(proverKey, E)

	// Verifier: final check
	ok := ver.EvalStep2(bar, E, opening, C)
	if !ok {
		t.Fatal("LVCS.EvalStep2 failed: proof rejected")
	}

	t.Log("LVCS commit+eval succeeded")
}
