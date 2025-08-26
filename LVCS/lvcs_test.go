package lvcs

import (
	"crypto/rand"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	decs "vSIS-Signature/DECS"
)

func lvcsParams(ringQ *ring.Ring, eta, ell int) decs.Params {
	return decs.Params{Degree: int(ringQ.N - 1), Eta: eta, NonceBytes: 16}
}

func TestLVCS_CommitAndEval_Accepts(t *testing.T) {
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, err := ring.NewRing(N, moduli)
	if err != nil {
		t.Fatalf("ring.NewRing: %v", err)
	}

	nrows := 4
	ell := 8
	ncols := int(ringQ.N) - ell
	params := lvcsParams(ringQ, 2, ell)

	rows := make([][]uint64, nrows)
	q0 := ringQ.Modulus[0]
	for j := 0; j < nrows; j++ {
		rows[j] = make([]uint64, ncols)
		for i := 0; i < ncols; i++ {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			rows[j][i] = x.Uint64()
		}
	}

	root, proverKey, err := CommitInitWithParams(ringQ, rows, ell, params)
	if err != nil {
		t.Fatalf("CommitInit: %v", err)
	}

	ver := NewVerifierWithParams(ringQ, nrows, params, ncols)
	ver.CommitStep1(root)

	R := CommitFinish(proverKey, ver.Gamma)
	ver.CommitStep2(R)

	m := 3
	C := make([][]uint64, m)
	for k := 0; k < m; k++ {
		C[k] = make([]uint64, nrows)
		for j := 0; j < nrows; j++ {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			C[k][j] = x.Uint64()
		}
	}

	bar := EvalInit(ringQ, proverKey, C)

	E := ver.ChooseE(ell, ncols)

	opening := EvalFinish(proverKey, E)

	if !ver.EvalStep2(bar, E, opening.DECSOpen, C) {
		t.Fatal("LVCS EvalStep2 rejected a valid proof")
	}
}

func TestLVCS_Rejects_MismatchedE_AndHeadIndex(t *testing.T) {
	N := 1 << 11
	moduli := []uint64{(1<<32 - (1 << 20) + 1)}
	ringQ, _ := ring.NewRing(N, moduli)

	nrows := 3
	ell := 4
	ncols := int(ringQ.N) - ell
	params := lvcsParams(ringQ, 2, ell)

	q0 := ringQ.Modulus[0]
	rows := make([][]uint64, nrows)
	for j := 0; j < nrows; j++ {
		rows[j] = make([]uint64, ncols)
		for i := 0; i < ncols; i++ {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			rows[j][i] = x.Uint64()
		}
	}
	root, proverKey, _ := CommitInitWithParams(ringQ, rows, ell, params)
	ver := NewVerifierWithParams(ringQ, nrows, params, ncols)
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

	E := ver.ChooseE(ell, ncols)
	opening := EvalFinish(proverKey, E)
	Ebad := append([]int(nil), E...)
	if Ebad[0] == ncols {
		Ebad[0]++
	} else {
		Ebad[0] = ncols
	}
	if ver.EvalStep2(bar, Ebad, opening.DECSOpen, C) {
		t.Fatal("accepted opening for mismatched E")
	}

	Ehead := append([]int(nil), E...)
	Ehead[0] = 0
	openHead := EvalFinish(proverKey, Ehead)
	if ver.EvalStep2(bar, Ehead, openHead.DECSOpen, C) {
		t.Fatal("accepted opening containing a head index")
	}
}
