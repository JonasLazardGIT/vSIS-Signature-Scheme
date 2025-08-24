package PIOP

import (
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
	prof "vSIS-Signature/prof"
)

type DecompCols struct {
	D   []*ring.Poly   // [LS] digits
	T   []*ring.Poly   // [LS+1] remainders
	Bit [][]*ring.Poly // [LS][W] digit bits
}

// appendDecompositionColumns allocates witness columns for digits, remainders and bits.
func appendDecompositionColumns(r *ring.Ring, LS, W int) DecompCols {
	defer prof.Track(time.Now(), "appendDecompositionColumns")
	mk := func() *ring.Poly { p := r.NewPoly(); r.NTT(p, p); return p }
	cols := DecompCols{
		D:   make([]*ring.Poly, LS),
		T:   make([]*ring.Poly, LS+1),
		Bit: make([][]*ring.Poly, LS),
	}
	for i := 0; i < LS; i++ {
		cols.D[i] = mk()
		cols.Bit[i] = make([]*ring.Poly, W)
		for u := 0; u < W; u++ {
			cols.Bit[i][u] = mk()
		}
	}
	for i := 0; i <= LS; i++ {
		cols.T[i] = mk()
	}
	return cols
}

// buildFparIntegerDecomp emits parallel rows for bitness, digit formation and remainder chain.
func buildFparIntegerDecomp(r *ring.Ring, Sqs *ring.Poly, spec BoundSpec, cols DecompCols) (Fpar []*ring.Poly) {
	defer prof.Track(time.Now(), "buildFparIntegerDecomp")
	q := r.Modulus[0]
	LS, W, R := spec.LS, spec.W, spec.R

	// Bitness and digit formation
	for l := 0; l < LS; l++ {
		for u := 0; u < W; u++ {
			Fpar = append(Fpar, bitnessPoly(r, cols.Bit[l][u]))
		}
		sum := r.NewPoly()
		tmp := r.NewPoly()
		for u := 0; u < W; u++ {
			scalePolyNTT(r, cols.Bit[l][u], (1<<uint(u))%q, tmp)
			r.Add(sum, tmp, sum)
		}
		p := r.NewPoly()
		r.Sub(cols.D[l], sum, p)
		Fpar = append(Fpar, p)
	}

	// Link T0 to Sqs: T0 - Sqs = 0
	p0 := r.NewPoly()
	r.Sub(cols.T[0], Sqs, p0)
	Fpar = append(Fpar, p0)

	// Remainder chain: T_l - D_l - R*T_{l+1} = 0 for l = 0..LS-1
	for l := 0; l < LS; l++ {
		tmp := r.NewPoly()
		scalePolyNTT(r, cols.T[l+1], R%q, tmp)
		r.Add(tmp, cols.D[l], tmp)
		p := r.NewPoly()
		r.Sub(cols.T[l], tmp, p)
		Fpar = append(Fpar, p)
	}

	// Final remainder T_LS should be zero
	Fpar = append(Fpar, cols.T[LS].CopyNew())
	return
}

// bitnessPoly returns polynomial enforcing bitness: b^2 - b, with Hadamard square in coeff domain.
func bitnessPoly(r *ring.Ring, b *ring.Poly) *ring.Poly {
	q := r.Modulus[0]
	coeff := r.NewPoly()
	r.InvNTT(b, coeff)
	for i := 0; i < r.N; i++ {
		v := coeff.Coeffs[0][i] % q
		t := (v * v) % q
		coeff.Coeffs[0][i] = (t + q - v) % q
	}
	r.NTT(coeff, coeff)
	return coeff
}
