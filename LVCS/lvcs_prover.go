package lvcs

import (
	"crypto/rand"
	"math/big"

	decs "vSIS-Signature/DECS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// ProverKey holds everything the prover needs between Commit and Eval.
type ProverKey struct {
	RingQ      *ring.Ring   // so we can grab q later without touching unexported decs.Prover.ringQ
	DecsProver *decs.Prover // underlying DECS prover

	RowData  [][]uint64 // original r_j
	MaskData [][]uint64 // random ̄r_j
}

// CommitInit – §4.1 steps 1–2:
// Lift each row vector to a degree-(N+ℓ−1) polynomial by appending ℓ random masks
// and commit all those polynomials via DECS.
func CommitInit(
	ringQ *ring.Ring,
	rows [][]uint64, // r_j in F_q^N
	ell int, // ℓ
) (
	root [32]byte,
	prover *ProverKey,
	err error,
) {
	nrows := len(rows)
	q0 := ringQ.Modulus[0]

	// 1a) sample masks ̄r_j ∈ F_q^ℓ
	masks := make([][]uint64, nrows)
	for j := range rows {
		masks[j] = make([]uint64, ell)
		for i := 0; i < ell; i++ {
			x, _ := rand.Int(rand.Reader, big.NewInt(int64(q0)))
			masks[j][i] = uint64(x.Int64())
		}
	}

	// 1b) interpolate each (r_j, mask_j) into P_j(X)
	polys := make([]*ring.Poly, nrows)
	for j := range rows {
		ncols := len(rows[j])
		if polys[j], err = interpolateRow(ringQ, rows[j], masks[j], ncols, ell); err != nil {
			return
		}
	}

	// 2) DECS.CommitInit
	dprover := decs.NewProver(ringQ, polys)
	if root, err = dprover.CommitInit(); err != nil {
		return
	}

	prover = &ProverKey{
		RingQ:      ringQ,
		DecsProver: dprover,
		RowData:    rows,
		MaskData:   masks,
	}
	return
}

// CommitFinish – §4.1 step 3:
// Later, open masked linear combinations on a small random subset EE of size ℓ.
func CommitFinish(
	prover *ProverKey,
	Gamma [][]uint64,
) []*ring.Poly {
	// nothing exported in decs.Prover needs ringQ here
	return prover.DecsProver.CommitStep2(Gamma)
}

// EvalInit – §4.1 step 1:
// Compute each \bar{v}_k = Σ_j C[k][j] · mask_j.
func EvalInit(
	ringQ *ring.Ring,
	prover *ProverKey,
	C [][]uint64, // C[k][j]
) [][]uint64 {
	nrows := len(prover.RowData)
	m := len(C)
	ell := len(prover.MaskData[0])
	q0 := ringQ.Modulus[0]

	bar := make([][]uint64, m)
	for k := 0; k < m; k++ {
		bar[k] = make([]uint64, ell)
		for j := 0; j < nrows; j++ {
			cij := C[k][j]
			row := prover.MaskData[j]
			for i := 0; i < ell; i++ {
				bar[k][i] = (bar[k][i] + cij*row[i]) % q0
			}
		}
	}
	return bar
}

// EvalFinish – §4.1 steps 3–4:
// Open the masked positions via DECS.EvalOpen.
func EvalFinish(
	prover *ProverKey,
	E []int,
) *decs.DECSOpening {
	return prover.DecsProver.EvalOpen(E)
}
