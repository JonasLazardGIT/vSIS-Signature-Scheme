package decs

// DECSOpening holds the data sent by the prover in DECS.Eval.
type DECSOpening struct {
	Indices []int      // the subset E ⊆ [0..N)
	Pvals   [][]uint64 // P_j(e) for each e∈E, j∈[0..r)
	Mvals   [][]uint64 // M_k(e) for each e∈E, k∈[0..η)
	Paths   [][][]byte // Merkle auth paths: Paths[t][level] = sibling hash
	Nonces  [][]byte   // ρ_e for each e∈E
}

// Params bundles the protocol parameters for DECS.
type Params struct {
	Degree     int // max degree d ≤ N-1
	Eta        int // number of mask polynomials η
	NonceBytes int // size of each nonce ρ_e in bytes
}

// DefaultParams provides legacy parameters for callers that do not
// explicitly configure DECS. It preserves the previous behaviour
// (Degree=4095, Eta=2, NonceBytes=16).
var DefaultParams = Params{Degree: 4095, Eta: 2, NonceBytes: 16}
