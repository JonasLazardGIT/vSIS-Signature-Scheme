package decs

import (
	"encoding/binary"
	"fmt"
	"os"

	"github.com/tuneinsight/lattigo/v4/ring"
)

var debugDECS = os.Getenv("DEBUG_DECS") != ""

// Verifier holds DECS verification parameters.
type Verifier struct {
	ringQ  *ring.Ring
	r      int
	params Params
}

// NewVerifierWithParams constructs a DECS verifier for r×η with the provided
// parameters. It panics if params.Degree is not in [0,N).
func NewVerifierWithParams(ringQ *ring.Ring, r int, params Params) *Verifier {
	if params.Degree < 0 || params.Degree >= int(ringQ.N) {
		panic("decs: invalid degree parameter")
	}
	return &Verifier{ringQ: ringQ, r: r, params: params}
}

// NewVerifier constructs a verifier with DefaultParams and the supplied η.
// It is retained for backwards compatibility.
func NewVerifier(ringQ *ring.Ring, r, eta int) *Verifier {
	params := DefaultParams
	params.Eta = eta
	if params.Degree >= int(ringQ.N) {
		params.Degree = int(ringQ.N) - 1
	}
	return NewVerifierWithParams(ringQ, r, params)
}

// DeriveGamma runs step 2 (commit → Γ).
func (v *Verifier) DeriveGamma(root [32]byte) [][]uint64 {
	q := v.ringQ.Modulus[0]
	return DeriveGamma(root, v.params.Eta, v.r, q)
}

// VerifyCommit checks deg R_k <= Degree (DECS §3 Step 3).
func (v *Verifier) VerifyCommit(root [32]byte, R []*ring.Poly, Gamma [][]uint64) bool {
	for _, p := range R {
		coeffs := p.Coeffs[0] // coeff domain
		for i := v.params.Degree + 1; i < len(coeffs); i++ {
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
	// NTT(R)
	Re := make([]*ring.Poly, v.params.Eta)
	for k := 0; k < v.params.Eta; k++ {
		Re[k] = v.ringQ.NewPoly()
		v.ringQ.NTT(R[k], Re[k])
	}

	mod := v.ringQ.Modulus[0]

    for t, idx := range open.Indices {
        // rebuild leaf (uint64 packing for values)
        buf := make([]byte, 8*(v.r+v.params.Eta)+4+v.params.NonceBytes)
        off := 0
        for j := 0; j < v.r; j++ {
            binary.LittleEndian.PutUint64(buf[off:], open.Pvals[t][j])
            off += 8
        }
        for k := 0; k < v.params.Eta; k++ {
            binary.LittleEndian.PutUint64(buf[off:], open.Mvals[t][k])
            off += 8
        }
        binary.LittleEndian.PutUint32(buf[off:], uint32(idx))
        off += 4
        copy(buf[off:], open.Nonces[t])

		if !VerifyPath(buf, open.Paths[t], root, idx) {
			if debugDECS {
				fmt.Printf("[DECS] FAIL Merkle at t=%d idx=%d (path mismatch)\n", t, idx)
			}
			return false
		}

		for k := 0; k < v.params.Eta; k++ {
			lhs := Re[k].Coeffs[0][idx]
			rhs := open.Mvals[t][k]
			for j := 0; j < v.r; j++ {
				mul := (open.Pvals[t][j] * Gamma[k][j]) % mod
				rhs = (rhs + mul) % mod
			}
			if lhs != rhs {
				if debugDECS {
					fmt.Printf("[DECS] FAIL masked-relation at t=%d idx=%d k=%d\n", t, idx, k)
					fmt.Printf("       lhs=NTT(R[%d])[idx]=%d, rhs=Mvals+ΣΓ·P=%d\n", k, lhs, rhs)
					fmt.Printf("       Pvals[t]=%v\n", open.Pvals[t])
					fmt.Printf("       Gamma[k]=%v\n", Gamma[k])
				}
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
		if debugDECS {
			fmt.Printf("[DECS] FAIL set-size: |open.Indices|=%d != |E|=%d\n", len(open.Indices), len(E))
		}
		return false
	}
	seen := make(map[int]int, len(E))
	for _, x := range E {
		seen[x]++
	}
	for _, y := range open.Indices {
		if seen[y] == 0 {
			if debugDECS {
				fmt.Printf("[DECS] FAIL set-bind: extra index %d in opening\n", y)
			}
			return false
		}
		seen[y]--
	}
	return v.VerifyEval(root, Gamma, R, open)
}
