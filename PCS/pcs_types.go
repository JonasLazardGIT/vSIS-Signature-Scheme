// Package pcs implements the SmallWood Polynomial Commitment Scheme
// on top of the LVCS and DECS layers (see Fig. 3 in the paper).
package pcs

import (
	lvcs "vSIS-Signature/LVCS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// ProverKey holds all the prover’s secret state between Commit and Eval.
type ProverKey struct {
	RingQ    *ring.Ring      // NTT ring context
	LvcsKey  *lvcs.ProverKey // LVCS prover state from CommitInit
	Dj       []int           // degrees of each P_j(X)
	Mu       int             // row‐block parameter μ
	EllPrime int             // zero‐knowledge mask count ℓ′
	Nu       []int           // number of columns ν_j per P_j
}

// Commitment is the public PCS commitment (just a root).
type Commitment struct {
	Root [32]byte
}
