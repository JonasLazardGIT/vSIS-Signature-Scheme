package decs

// DECSOpening holds the data sent by the prover in DECS.Eval.
type DECSOpening struct {
	Indices []int      // the subset E ⊆ [0..N)
	Pvals   [][]uint64 // P_j(e) for each e∈E, j∈[0..r)
	Mvals   [][]uint64 // M_k(e) for each e∈E, k∈[0..η)
	Paths   [][][]byte // Merkle auth paths: Paths[t][level] = sibling hash
	Nonces  [][]byte   // ρ_e for each e∈E
}

// Protocol constants – must match prover & verifier.
const (
	Eta        = 2    // number of mask polys η
	Degree     = 4095 // max degree d
	NonceBytes = 16   // size of each ρ
)
