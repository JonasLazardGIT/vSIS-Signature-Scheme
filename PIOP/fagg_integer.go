package PIOP

import "github.com/tuneinsight/lattigo/v4/ring"

type GlobCarry struct {
	C     []*ring.Poly
	CBits [][]*ring.Poly
}

// appendGlobalCarrys allocates global carry polys and bit columns.
func appendGlobalCarrys(r *ring.Ring, LS, Wc int) GlobCarry {
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
	return
}

// buildFaggIntegerSum constructs aggregated limb-sum rows without slack.
func buildFaggIntegerSum(r *ring.Ring, spec BoundSpec, _ uint64, S0inv uint64, cols DecompCols, g GlobCarry) (Fagg []*ring.Poly) {
	LS, R, q := spec.LS, spec.R, r.Modulus[0]
	for l := 0; l < LS; l++ {
		p := r.NewPoly()
		r.Add(p, cols.D[l], p)
		r.Add(p, g.C[l], p)

		Bscaled := (spec.Beta2[l] % q) * (S0inv % q) % q
		constB := constPolyNTT(r, Bscaled)
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
func buildFaggIntegerSumDelta(r *ring.Ring, spec BoundSpec, _ uint64, S0inv uint64, cols DecompCols, g GlobCarry, s GlobSlack) (Fagg []*ring.Poly) {
	LS, R, q := spec.LS, spec.R, r.Modulus[0]
	for l := 0; l < LS; l++ {
		coeff := r.NewPoly()
		r.InvNTT(cols.D[l], coeff)
		sum := uint64(0)
		for j := 0; j < r.N; j++ {
			sum = (sum + coeff.Coeffs[0][j]) % q
		}
		constD := constPolyNTT(r, (sum*S0inv)%q)

		p := r.NewPoly()
		r.Add(p, constD, p)
		r.Add(p, g.C[l], p)
		r.Add(p, s.D[l], p)

		Bscaled := (spec.Beta2[l] % q) * (S0inv % q) % q
		constB := constPolyNTT(r, Bscaled)
		r.Sub(p, constB, p)

		tmp := r.NewPoly()
		scalePolyNTT(r, g.C[l+1], R%q, tmp)
		r.Sub(p, tmp, p)
		Fagg = append(Fagg, p)
	}
	return
}
