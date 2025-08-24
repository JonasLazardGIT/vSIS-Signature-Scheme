package PIOP

import "math/big"

// BoundSpec holds parameters for the integer L2 bound proof.
type BoundSpec struct {
	Q     uint64   // modulus q
	R     uint64   // radix R = 2^W (must satisfy R < q)
	W     int      // bits per limb
	LS    int      // number of limbs for beta^2 in base R
	Beta2 []uint64 // beta^2 limbs (little-endian) in base R
	Ell   int      // blinding points per row
}

// NewBoundSpec returns the radix decomposition of beta^2 in base R=2^w.
// R must satisfy R < q and beta^2 < q.
func NewBoundSpec(q, beta uint64, w, ell int) BoundSpec {
	if (q & 1) == 0 {
		panic("q must be odd")
	}
	if ell < 1 {
		panic("ell must be ≥ 1")
	}
	R := uint64(1) << uint(w)
	if R >= q {
		panic("R >= q")
	}
	bb := new(big.Int).SetUint64(beta)
	bb.Mul(bb, bb) // beta^2
	if bb.Cmp(new(big.Int).SetUint64(q)) >= 0 {
		panic("beta^2 >= q")
	}

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
	return BoundSpec{Q: q, R: R, W: w, LS: len(limbs), Beta2: limbs, Ell: ell}
}
