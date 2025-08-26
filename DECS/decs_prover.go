package decs

import (
	"crypto/rand"
	"crypto/sha256"
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Prover encapsulates the prover state for DECS.
type Prover struct {
	ringQ  *ring.Ring
	P      []*ring.Poly // r input polys (coeff form)
	M      []*ring.Poly // η mask polys (coeff form)
	nonces [][]byte     // length N nonces
	mt     *MerkleTree
	Pvals  []*ring.Poly // NTT(P)
	Mvals  []*ring.Poly // NTT(M)
	root   [32]byte
	R      []*ring.Poly // η output polys in coeff form
	params Params
}

// NewProverWithParams returns a new DECS prover for polynomials P and the
// provided protocol parameters. It panics if params.Degree is not in [0, N).
func NewProverWithParams(ringQ *ring.Ring, P []*ring.Poly, params Params) *Prover {
	if params.Degree < 0 || params.Degree >= int(ringQ.N) {
		panic("decs: invalid degree parameter")
	}
	return &Prover{ringQ: ringQ, P: P, params: params}
}

// NewProver returns a new DECS prover with DefaultParams. It is provided for
// backwards compatibility.
func NewProver(ringQ *ring.Ring, P []*ring.Poly) *Prover {
	params := DefaultParams
	if params.Degree >= int(ringQ.N) {
		params.Degree = int(ringQ.N) - 1
	}
	return NewProverWithParams(ringQ, P, params)
}

// CommitInit does DECS.Commit step 1: sample M, nonces; build Merkle tree; NTT(P,M).
func (pr *Prover) CommitInit() ([32]byte, error) {
	r := len(pr.P)
	N := pr.ringQ.N

	// sampler
	prng, err := utils.NewPRNG()
	if err != nil {
		return [32]byte{}, err
	}
	us := ring.NewUniformSampler(prng, pr.ringQ)

	// 1a) sample η mask polys
	pr.M = make([]*ring.Poly, pr.params.Eta)
	for k := 0; k < pr.params.Eta; k++ {
		pr.M[k] = pr.ringQ.NewPoly()
		us.Read(pr.M[k])
		for i := pr.params.Degree + 1; i < N; i++ {
			pr.M[k].Coeffs[0][i] = 0
		}
	}

	// 1b) NTT-transform P and M
	pr.Pvals = make([]*ring.Poly, r)
	for j := range pr.P {
		pr.Pvals[j] = pr.ringQ.NewPoly()
		pr.ringQ.NTT(pr.P[j], pr.Pvals[j])
	}
	pr.Mvals = make([]*ring.Poly, pr.params.Eta)
	for k := 0; k < pr.params.Eta; k++ {
		pr.Mvals[k] = pr.ringQ.NewPoly()
		pr.ringQ.NTT(pr.M[k], pr.Mvals[k])
	}

	// 1c) build leaves
	leaves := make([][]byte, N)
	pr.nonces = make([][]byte, N)
	for i := 0; i < N; i++ {
		buf := make([]byte, 4*(r+pr.params.Eta)+4+pr.params.NonceBytes)
		off := 0
		for j := 0; j < r; j++ {
			binary.LittleEndian.PutUint32(buf[off:], uint32(pr.Pvals[j].Coeffs[0][i]))
			off += 4
		}
		for k := 0; k < pr.params.Eta; k++ {
			binary.LittleEndian.PutUint32(buf[off:], uint32(pr.Mvals[k].Coeffs[0][i]))
			off += 4
		}
		binary.LittleEndian.PutUint32(buf[off:], uint32(i))
		off += 4
		rho := make([]byte, pr.params.NonceBytes)
		rand.Read(rho)
		copy(buf[off:], rho)
		pr.nonces[i] = rho

		// store the raw buffer; BuildMerkleTree will hash it
		leaves[i] = append([]byte(nil), buf...)
	}

	// 1d) Merkle tree
	pr.mt = BuildMerkleTree(leaves)
	pr.root = pr.mt.Root()

	return pr.root, nil
}

// CommitStep2 does DECS.Commit steps 2+3: given Γ, compute R_k = M_k + Σ_j Γ[k][j]*P_j.
func (pr *Prover) CommitStep2(Gamma [][]uint64) []*ring.Poly {
	r := len(pr.P)
	pr.R = make([]*ring.Poly, pr.params.Eta)
	tmp := pr.ringQ.NewPoly()
	tmp2 := pr.ringQ.NewPoly()

	for k := 0; k < pr.params.Eta; k++ {
		// inv-NTT(M_k) → tmp
		pr.ringQ.InvNTT(pr.Mvals[k], tmp)
		pr.R[k] = tmp.CopyNew()
		for j := 0; j < r; j++ {
			pr.ringQ.InvNTT(pr.Pvals[j], tmp)
			pr.ringQ.MulScalar(tmp, Gamma[k][j], tmp2) // tmp2 = tmp * Γ[k][j]
			pr.ringQ.Add(pr.R[k], tmp2, pr.R[k])       // R[k] += tmp2
		}
		// keep R[k] in coefficient form; verifier will NTT as needed
	}
	return pr.R
}

// EvalOpen does DECS.Eval step 1: given E, returns Pvals,Mvals,Paths,Nonces.
func (pr *Prover) EvalOpen(E []int) *DECSOpening {
	r := len(pr.P)
	open := &DECSOpening{
		Indices: E,
		Pvals:   make([][]uint64, len(E)),
		Mvals:   make([][]uint64, len(E)),
		Paths:   make([][][]byte, len(E)),
		Nonces:  make([][]byte, len(E)),
	}
	for t, idx := range E {
		open.Pvals[t] = make([]uint64, r)
		for j := 0; j < r; j++ {
			open.Pvals[t][j] = pr.Pvals[j].Coeffs[0][idx]
		}
		open.Mvals[t] = make([]uint64, pr.params.Eta)
		for k := 0; k < pr.params.Eta; k++ {
			open.Mvals[t][k] = pr.Mvals[k].Coeffs[0][idx]
		}
		open.Paths[t] = pr.mt.Path(idx)
		open.Nonces[t] = pr.nonces[idx]
	}
	return open
}

// DeriveGamma expands root→η×r matrix Γ with entries uniform in [0,q).
// Uses SHA256(root || ctr) as a PRF and 64-bit rejection sampling for exact uniformity.
func DeriveGamma(root [32]byte, eta, r int, q uint64) [][]uint64 {
	out := make([][]uint64, eta)
	var ctr uint64
	for k := 0; k < eta; k++ {
		out[k] = make([]uint64, r)
		for j := 0; j < r; j++ {
			for {
				var buf [40]byte
				copy(buf[:32], root[:])
				binary.LittleEndian.PutUint64(buf[32:], ctr)
				h := sha256.Sum256(buf[:])
				x := binary.LittleEndian.Uint64(h[:8])
				ctr++
				limit := (^uint64(0) / q) * q
				if x < limit {
					out[k][j] = x % q
					break
				}
			}
		}
	}
	return out
}
