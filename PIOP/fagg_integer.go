package PIOP

import (
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
	prof "vSIS-Signature/prof"
)

type GlobCarry struct {
	C     []*ring.Poly
	CBits [][]*ring.Poly
}

// appendGlobalCarrys allocates global carry polys and bit columns.
func appendGlobalCarrys(r *ring.Ring, LS, Wc int) GlobCarry {
	defer prof.Track(time.Now(), "appendGlobalCarrys")
	mk := func() *ring.Poly { p := r.NewPoly(); r.NTT(p, p); return p }
	g := GlobCarry{
		C:     make([]*ring.Poly, LS+1),
		CBits: make([][]*ring.Poly, LS+1),
	}
	for l := 0; l <= LS; l++ {
		g.C[l] = mk()
		g.CBits[l] = make([]*ring.Poly, Wc)
		for u := 0; u < Wc; u++ {
			g.CBits[l][u] = mk()
		}
	}
	return g
}

// buildFparGlobCarryBits ties each carry to its bit decomposition.
func buildFparGlobCarryBits(r *ring.Ring, S0inv uint64, Wc int, g GlobCarry) (Fpar []*ring.Poly) {
	defer prof.Track(time.Now(), "buildFparGlobCarryBits")
	q := r.Modulus[0]
	for l := 0; l < len(g.C); l++ {
		for u := 0; u < Wc; u++ {
			Fpar = append(Fpar, bitnessPoly(r, g.CBits[l][u]))
		}
		sum := r.NewPoly()
		tmp := r.NewPoly()
		for u := 0; u < Wc; u++ {
			scale := (S0inv * (1 << uint(u))) % q
			scalePolyNTT(r, g.CBits[l][u], scale, tmp)
			r.Add(sum, tmp, sum)
		}
		p := r.NewPoly()
		r.Sub(g.C[l], sum, p)
		Fpar = append(Fpar, p)
	}
	// Boundary carries must be zero
	Fpar = append(Fpar, g.C[0].CopyNew())
	Fpar = append(Fpar, g.C[len(g.C)-1].CopyNew())
	return
}

// buildFaggIntegerSum constructs aggregated limb-sum rows without slack.
func buildFaggIntegerSum(r *ring.Ring, spec BoundSpec, S0, S0inv uint64, cols DecompCols, g GlobCarry, omega []uint64) (Fagg []*ring.Poly) {
	defer prof.Track(time.Now(), "buildFaggIntegerSum")
	LS, R, q := spec.LS, spec.R, r.Modulus[0]
	scratch := r.NewPoly()
	for l := 0; l < LS; l++ {
		sumOmegaD := sumEvals(r, cols.D[l], omega, scratch)
		constD := makeConstRow(r, modMul(sumOmegaD%q, S0inv%q, q))

		p := r.NewPoly()
		r.Add(p, constD, p)
		r.Add(p, g.C[l], p)

		Bscaled := modMul(spec.Beta2[l]%q, S0inv%q, q)
		constB := makeConstRow(r, Bscaled)
		r.Sub(p, constB, p)

		tmp := r.NewPoly()
		scalePolyNTT(r, g.C[l+1], R%q, tmp)
		r.Sub(p, tmp, p)
		Fagg = append(Fagg, p)
	}
	return
}

// buildFparGlobSlackBits ties slack digits to their bit decompositions.
func buildFparGlobSlackBits(r *ring.Ring, S0inv uint64, W int, s GlobSlack) (Fpar []*ring.Poly) {
	defer prof.Track(time.Now(), "buildFparGlobSlackBits")
	q := r.Modulus[0]
	for l := 0; l < len(s.D); l++ {
		for u := 0; u < W; u++ {
			Fpar = append(Fpar, bitnessPoly(r, s.DBits[l][u]))
		}
		sum := r.NewPoly()
		tmp := r.NewPoly()
		for u := 0; u < W; u++ {
			scale := (S0inv * (1 << uint(u))) % q
			scalePolyNTT(r, s.DBits[l][u], scale, tmp)
			r.Add(sum, tmp, sum)
		}
		p := r.NewPoly()
		r.Sub(s.D[l], sum, p)
		Fpar = append(Fpar, p)
	}
	return
}

// buildFaggIntegerSumDelta includes slack digits in the aggregated rows.
func buildFaggIntegerSumDelta(r *ring.Ring, spec BoundSpec, S0, S0inv uint64, cols DecompCols, g GlobCarry, s GlobSlack, omega []uint64) (Fagg []*ring.Poly) {
	defer prof.Track(time.Now(), "buildFaggIntegerSumDelta")
	LS, R, q := spec.LS, spec.R, r.Modulus[0]
	scratch := r.NewPoly()
	for l := 0; l < LS; l++ {
		sumOmegaD := sumEvals(r, cols.D[l], omega, scratch)
		constD := makeConstRow(r, modMul(sumOmegaD%q, S0inv%q, q))

		p := r.NewPoly()
		r.Add(p, constD, p)
		r.Add(p, g.C[l], p)
		r.Add(p, s.D[l], p)

		Bscaled := modMul(spec.Beta2[l]%q, S0inv%q, q)
		constB := makeConstRow(r, Bscaled)
		r.Sub(p, constB, p)

		tmp := r.NewPoly()
		scalePolyNTT(r, g.C[l+1], R%q, tmp)
		r.Sub(p, tmp, p)
		Fagg = append(Fagg, p)
	}
	return
}
