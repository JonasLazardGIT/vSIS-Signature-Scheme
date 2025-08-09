package PIOP

import "math/big"

// BoundSpec holds parameters for the integer L2 bound proof.
type BoundSpec struct {
	Q     uint64   // modulus q
	R     uint64   // radix R = 2^W (must satisfy R < q)
	W     int      // bits per limb
	LS    int      // number of limbs for beta^2 in base R
	Beta2 []uint64 // beta^2 limbs (little-endian) in base R
}

// NewBoundSpec returns the radix decomposition of beta^2 in base R=2^w.
// R must satisfy R < q.
func NewBoundSpec(q, beta uint64, w int) BoundSpec {
	R := uint64(1) << uint(w)
	bb := new(big.Int).SetUint64(beta)
	bb.Mul(bb, bb) // beta^2

	BR := new(big.Int).SetUint64(R)
	zero := new(big.Int)
	var limbs []uint64
	tmp := new(big.Int).Set(bb)
	for tmp.Cmp(zero) > 0 {
		rem := new(big.Int)
		tmp.DivMod(tmp, BR, rem)
		limbs = append(limbs, rem.Uint64())
	}
	if len(limbs) == 0 {
		limbs = []uint64{0}
	}
	return BoundSpec{Q: q, R: R, W: w, LS: len(limbs), Beta2: limbs}
}
