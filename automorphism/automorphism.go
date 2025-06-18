package automorphism

import "github.com/tuneinsight/lattigo/v4/ring"

// ApplyNTT applies σ_ρ (ρ odd) to an NTT-domain polynomial.
// Both input and output stay in plain NTT form.
func ApplyNTT(r *ring.Ring, pIn *ring.Poly, rho int) *ring.Poly {
	if rho&1 == 0 {
		panic("rho must be odd")
	}

	d := r.N
	twoD := 2 * d
	q0 := r.Modulus[0]

	coeff := r.NewPoly()
	r.InvNTT(pIn, coeff)

	out := r.NewPoly()
	for j := 0; j < d; j++ {
		idx := (j * rho) % twoD
		if idx < d {
			out.Coeffs[0][j] = coeff.Coeffs[0][idx]
		} else {
			v := coeff.Coeffs[0][idx-d]
			out.Coeffs[0][j] = (q0 - v) % q0
		}
	}

	r.NTT(out, out)
	return out
}
