package decs

import (
	"encoding/binary"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Verifier holds DECS verification parameters.
type Verifier struct {
	ringQ  *ring.Ring
	r, eta int
}

// NewVerifier constructs a DECS verifier for r×η.
func NewVerifier(ringQ *ring.Ring, r, eta int) *Verifier {
	return &Verifier{ringQ: ringQ, r: r, eta: eta}
}

// DeriveGamma runs step 2 (commit → Γ).
func (v *Verifier) DeriveGamma(root [32]byte) [][]uint64 {
	q := v.ringQ.Modulus[0]
	return DeriveGamma(root, v.eta, v.r, q)
}

// VerifyCommit checks deg R_k <= Degree (DECS §3 Step 3).
func (v *Verifier) VerifyCommit(root [32]byte, R []*ring.Poly, Gamma [][]uint64) bool {
	for _, p := range R {
		coeffs := p.Coeffs[0] // coeff domain
		for i := Degree + 1; i < len(coeffs); i++ {
			if coeffs[i] != 0 {
				return false // degree too large
			}
		}
	}
	return true
}

// VerifyEval runs DECS.Eval checks: Merkle paths + masked relation.
func (v *Verifier) VerifyEval(
	root [32]byte, Gamma [][]uint64, R []*ring.Poly,
	open *DECSOpening,
) bool {
	// pre-NTT R
	Re := make([]*ring.Poly, v.eta)
	for k := 0; k < v.eta; k++ {
		Re[k] = v.ringQ.NewPoly()
		v.ringQ.NTT(R[k], Re[k])
	}

	mod := v.ringQ.Modulus[0]

	for t, idx := range open.Indices {
		// rebuild leaf
		buf := make([]byte, 4*(v.r+v.eta)+4+NonceBytes)
		off := 0
		for j := 0; j < v.r; j++ {
			binary.LittleEndian.PutUint32(buf[off:], uint32(open.Pvals[t][j]))
			off += 4
		}
		for k := 0; k < v.eta; k++ {
			binary.LittleEndian.PutUint32(buf[off:], uint32(open.Mvals[t][k]))
			off += 4
		}
		binary.LittleEndian.PutUint32(buf[off:], uint32(idx))
		off += 4
		copy(buf[off:], open.Nonces[t])

		if !VerifyPath(buf, open.Paths[t], root, idx) {
			return false
		}

		for k := 0; k < v.eta; k++ {
			lhs := Re[k].Coeffs[0][idx]
			rhs := open.Mvals[t][k]
			for j := 0; j < v.r; j++ {
				mul := (open.Pvals[t][j] * Gamma[k][j]) % mod
				rhs = (rhs + mul) % mod
			}
			if lhs != rhs {
				return false
			}
		}
	}
	return true
}

// VerifyEvalAt enforces that the prover opened exactly the challenged set E,
// then runs the standard DECS checks.
func (v *Verifier) VerifyEvalAt(
	root [32]byte, Gamma [][]uint64, R []*ring.Poly,
	open *DECSOpening, E []int,
) bool {
	if len(open.Indices) != len(E) {
		return false
	}
	seen := make(map[int]int, len(E))
	for _, x := range E {
		seen[x]++
	}
	for _, y := range open.Indices {
		if seen[y] == 0 {
			return false
		}
		seen[y]--
	}
	return v.VerifyEval(root, Gamma, R, open)
}
