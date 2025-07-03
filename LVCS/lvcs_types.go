package lvcs

import decs "vSIS-Signature/DECS"

// LVCSProverKey holds the prover’s secret state.
type LVCSProverKey struct {
	// keyDECS is the opening key/state returned by DECS.Commit.
	KeyDECS decs.Prover
	// RowData[j] = original vector r_j ∈ F_q^{ncols}.
	RowData [][]uint64
	// MaskData[j] = random mask vector  ̄r_j ∈ F_q^ℓ.
	MaskData [][]uint64
}

// LVCSOpening is what the prover sends in LVCS.Eval.
type LVCSOpening struct {
	// BarTargets[k] =  ̄v_k = ∑_j c_{k,j} · ̄r_j
	BarTargets [][]uint64
	// DECSOpen contains P|E, nonces, and Merkle paths.
	DECSOpen *decs.DECSOpening
	// SubsetE is the set E ⊆ {0,…,ncols-1} the verifier challenged.
	SubsetE []int
}
