package Preimage_Sampler

import (
	"fmt"
	"log"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// ---------------------------------------------------------------------
// Sample2zField ⇐ ZSampleSigma2x2
// ---------------------------------------------------------------------
func Sample2zField(
	a, b, d *CyclotomicFieldElem, // in EVAL domain
	c0, c1 *CyclotomicFieldElem, // in COEFF domain
	n int,
	modulus uint64,
	prec uint,
) (q0, q1 *CyclotomicFieldElem) {

	// 1) draw q1 ← SampleFZBig(d, c1)
	dCoeff := FloatToCoeffNegacyclic(d, prec) // d is in EVAL, convert to COEFF
	// fmt.Printf("Sample2zField: d = %s\n", dCoeff.Coeffs[0].Real.Text('g', 10))
	// fmt.Printf("Sample2zField: c1 = %s\n", c1.Coeffs[0].Real.Text('g', 10))
	q1 = SampleFZBig(dCoeff, c1, modulus, prec)

	// 2) delta = q1 – c1  (both in COEFF), then lift to EVAL
	delta := FieldSubBig(q1, c1)

	delta = FloatToEvalNegacyclic(delta, prec) // delta is now in EVAL

	// 3) invD = conj(d) / |d|²   (remains in EVAL)
	invD, norms := FieldInverseDiagWithNorm(d)
	invD = FieldScalarDiv(invD, norms)
	invD.Domain = Eval

	//! Verify that invD is truly the inverse of d (i.e. d * invD ≈ 1 in Eval)
	prod := FieldMulBig(d, invD) // pointwise product in Eval
	prod.Domain = Eval

	// Build a “1 + 0i” constant for comparison
	oneReal := new(big.Float).SetPrec(prec).SetFloat64(1.0)
	zero := new(big.Float).SetPrec(prec).SetFloat64(0.0)
	one := &BigComplex{Real: oneReal, Imag: zero}

	// Tolerance for big.Float comparison (adjust as needed)
	tol := new(big.Float).SetPrec(prec).SetFloat64(1e-30)

	for i := 0; i < n; i++ {
		// compute difference prod.Coeffs[i] − (1 + 0i)
		diffReal := new(big.Float).Sub(prod.Coeffs[i].Real, one.Real)
		diffImag := new(big.Float).Sub(prod.Coeffs[i].Imag, one.Imag)

		// take absolute values
		absReal := new(big.Float).Abs(diffReal)
		absImag := new(big.Float).Abs(diffImag)

		if absReal.Cmp(tol) > 0 || absImag.Cmp(tol) > 0 {
			log.Panicf("invD check failed at index %d: d*invD = %s + %si", i,
				prod.Coeffs[i].Real.Text('g', 10),
				prod.Coeffs[i].Imag.Text('g', 10))
		}
	}
	//! End of invD sanity check
	tmp := FieldMulBig(invD, delta)
	tmp.Domain = Eval
	scaledDelta := FieldMulBig(b, tmp)

	// scaledDelta.Coeffs is still Eval‐domain; fix the flag:
	scaledDelta.Domain = Eval
	scaledDelta = FloatToCoeffNegacyclic(scaledDelta, prec)

	// 5) c0′ = c0 + scaledDelta   (all COEFF)
	c0p := FieldAddBig(c0, scaledDelta) //! ADD OR SUB ?

	// 6) cond = b ⋅ d⁻¹ ⋅ bᵀ  (compute in EVAL, then to COEFF)
	//   - form true Hermitian transpose of b in EVAL
	Bhermit := HermitianTransposeFieldElem(b) // now Domain=Eval

	// - multiply out in EVAL
	tmpCond := FieldMulBig(invD, Bhermit)
	tmpCond.Domain = Eval

	condEval := FieldMulBig(b, tmpCond)
	condEval.Domain = Eval
	condEval = FloatToCoeffNegacyclic(condEval, prec) // now Domain=Coeff

	// 7) aPr = a – cond   (both in COEFF)
	aPr := a.Copy()
	aPr = FloatToCoeffNegacyclic(aPr, prec) // a is in EVAL, convert to COEFF
	aPr = FieldSubBig(aPr, condEval)

	// 8) q0 ← SampleFZBig(aPr, c0′)
	fmt.Printf("Sample2zField: aPr = %s\n", aPr.Coeffs[0].Real.Text('g', 10))
	fmt.Printf("Sample2zField: c0′ = %s\n", c0p.Coeffs[0].Real.Text('g', 10))

	q0 = SampleFZBig(aPr, c0p, modulus, prec)
	return
}

func SampleFZBig(
	f *CyclotomicFieldElem, // in COEFF domain
	c *CyclotomicFieldElem, // in COEFF domain
	modulus uint64, // modulus of the ring Q
	prec uint,
) *CyclotomicFieldElem {
	m := f.N
	// fmt.Printf("SampleFZBig: m = %d\n", m)

	if m == 1 {
		varv, _ := f.Coeffs[0].Real.Float64()
		// fmt.Printf("SampleFZBig: varv = %f\n", varv)
		if varv < 0 {
			fmt.Printf("SampleFZBig: varv = %f < 0\n PANIC\n", varv)
		}
		std := math.Sqrt(varv)
		mean, _ := c.Coeffs[0].Real.Float64()
		dg := NewDiscreteGaussian(std)
		zInt := dg.Draw(mean) // fmt.Printf("SampleFZBig: zInt = %d\n", zInt)
		out := NewFieldElemBig(1, prec)
		out.Coeffs[0] = NewBigComplex(float64(zInt), 0, prec)
		return out
	}

	// 1) split even/odd parts
	f0 := f.Copy().ExtractEven()
	f1 := f.Copy().ExtractOdd()
	half := m / 2
	f0 = FloatToEvalNegacyclic(f0, prec)
	f1 = FloatToEvalNegacyclic(f1, prec)

	// 2) permute centers
	c0 := c.Copy().ExtractEven()
	c1 := c.Copy().ExtractOdd()

	// 3) two-block sample returns *two* halves
	q0, q1 := Sample2zField(
		f0, f1, f0,
		NewFieldElemBig(m/2, prec).SetCoeffs(c0),
		NewFieldElemBig(m/2, prec).SetCoeffs(c1),
		half,
		modulus, prec,
	)

	// 4) pack them into one length-m element
	out := NewFieldElemBig(m, prec)
	for i := 0; i < m/2; i++ {
		out.Coeffs[i] = q0.Coeffs[i]

		out.Coeffs[i+m/2] = q1.Coeffs[i]
	}

	// 5) inverse-permute to interleave even/odd
	InversePermuteFieldElem(out)

	return out
}

// ---------------------------------------------------------------------
// Permute / InversePermute  (if you ever need them)
// ---------------------------------------------------------------------
func Permute(p *Matrix[int64]) *Matrix[int64] {
	N := p.Rows
	out := NewMatrix[int64](N, 1)
	even, odd := 0, N/2
	for i := 0; i < N; i++ {
		if i%2 == 0 {
			out.Set(even, 0, p.At(i, 0))
			even++
		} else {
			out.Set(odd, 0, p.At(i, 0))
			odd++
		}
	}
	return out
}

// SamplePz implements Algorithm 4 “perturbation generation”
// – ringQ: the Rₙ ring (NTT domain)
// – s, alpha: real parameters
// – Ttilde = [r̂; ê] in EVAL (NTT) domain, each of length k
// – expectedLength = 2 + k
func SamplePz(
	ringQ *ring.Ring,
	s, alpha float64,
	Ttilde [2][]*ring.Poly,
	expectedLength int,
	prec uint,
) []*ring.Poly {
	n := ringQ.N
	k := len(Ttilde[0])

	// 1) Compute z = (1/α² − 1/s²)^(−1) with full `prec`
	one := new(big.Float).SetPrec(prec).SetFloat64(1.0)
	s2 := new(big.Float).SetPrec(prec).SetFloat64(s * s)
	α2 := new(big.Float).SetPrec(prec).SetFloat64(alpha * alpha)
	invS2 := new(big.Float).SetPrec(prec).Quo(one, s2)
	invA2 := new(big.Float).SetPrec(prec).Quo(one, α2)
	diff := new(big.Float).SetPrec(prec).Sub(invA2, invS2)
	zBig := new(big.Float).SetPrec(prec).Quo(one, diff)
	z := NewBigComplexFromFloat(zBig, big.NewFloat(0).SetPrec(prec))

	fmt.Printf("SamplePz: q = %d, s = %f, α = %f, z = %s\n", ringQ.Modulus[0], s, alpha, zBig.Text('g', 10))

	// 2) Initialize accumulators a,b,d in COEFF domain
	aFld := NewFieldElemBig(n, prec) // starts zero
	bFld := NewFieldElemBig(n, prec)
	dFld := NewFieldElemBig(n, prec)
	// set a = s², d = s²
	s2c := NewBigComplexFromFloat(s2, big.NewFloat(0).SetPrec(prec))
	aFld.Coeffs[0] = s2c.Copy()
	dFld.Coeffs[0] = s2c.Copy()

	// 3) Accumulate Σ r̂ᵀr̂, Σ r̂ᵀê, Σ êᵀê, each multiplied by z
	aFld = FloatToEvalNegacyclic(aFld, prec)
	bFld = FloatToEvalNegacyclic(bFld, prec)
	dFld = FloatToEvalNegacyclic(dFld, prec)

	for j := 0; j < k; j++ {
		// back to coeff domain
		rPoly := ringQ.NewPoly()
		ringQ.InvNTT(Ttilde[0][j], rPoly)
		ePoly := ringQ.NewPoly()
		ringQ.InvNTT(Ttilde[1][j], ePoly)

		rF := NegacyclicEvaluatePoly(rPoly, ringQ, prec)
		// rF = ToEvalNegacyclic(rF, ringQ, prec) // rF is now in EVAL
		eF := NegacyclicEvaluatePoly(ePoly, ringQ, prec)
		// eF = ToEvalNegacyclic(eF, ringQ, prec) // eF is now in EVAL

		rT := HermitianTransposeFieldElem(rF)
		eT := HermitianTransposeFieldElem(eF)

		for i := 0; i < n; i++ {
			tmp := rT.Coeffs[i].Mul(rF.Coeffs[i]).Mul(z)
			aFld.Coeffs[i] = aFld.Coeffs[i].Sub(tmp)
			tmp = rT.Coeffs[i].Mul(eF.Coeffs[i]).Mul(z)
			bFld.Coeffs[i] = bFld.Coeffs[i].Sub(tmp)
			tmp = eT.Coeffs[i].Mul(eF.Coeffs[i]).Mul(z)
			dFld.Coeffs[i] = dFld.Coeffs[i].Sub(tmp)
		}
	}

	// ! —––– begin per‐slot sanity‐check de l’Alg 4 en big.Float —––––
	const bigPrec = 256 // précision pour tous les big.Float
	tolBig := new(big.Float).SetPrec(bigPrec).SetFloat64(1e-30)
	sBig := new(big.Float).SetPrec(bigPrec).SetFloat64(s)
	for i := 0; i < n; i++ {
		// 1) reconstruire sumRR, sumRE, sumEE en big.Float
		sumRR := new(big.Float).SetPrec(bigPrec).SetFloat64(0.0)
		sumRE := new(big.Float).SetPrec(bigPrec).SetFloat64(0.0)
		sumEE := new(big.Float).SetPrec(bigPrec).SetFloat64(0.0)

		for j := 0; j < k; j++ {
			// Tirer r_j et e_j en entier
			rpoly := ringQ.NewPoly()
			ringQ.InvNTT(Ttilde[0][j], rpoly)
			epoly := ringQ.NewPoly()
			ringQ.InvNTT(Ttilde[1][j], epoly)

			rawR := rpoly.Coeffs[0][i]
			rawE := epoly.Coeffs[0][i]
			rInt := float64(UnsignedToSigned(rawR, ringQ.Modulus[0])) // r en float64 temporaire
			eInt := float64(UnsignedToSigned(rawE, ringQ.Modulus[0])) // e en float64 temporaire

			// Convertir en big.Float, même précision
			rBig := new(big.Float).SetPrec(bigPrec).SetFloat64(rInt)
			eBig := new(big.Float).SetPrec(bigPrec).SetFloat64(eInt)

			// sumRR += r^2
			r2 := new(big.Float).SetPrec(bigPrec).Mul(rBig, rBig)
			sumRR.Add(sumRR, r2)

			// sumRE += r*e
			re := new(big.Float).SetPrec(bigPrec).Mul(rBig, eBig)
			sumRE.Add(sumRE, re)

			// sumEE += e^2
			e2 := new(big.Float).SetPrec(bigPrec).Mul(eBig, eBig)
			sumEE.Add(sumEE, e2)
		}

		// 2) calculer aiBig = s^2 - z*sumRR
		s2Big := new(big.Float).SetPrec(bigPrec).Mul(sBig, sBig)    // s^2
		zsumRR := new(big.Float).SetPrec(bigPrec).Mul(zBig, sumRR)  // z * sumRR
		aiBig := new(big.Float).SetPrec(bigPrec).Sub(s2Big, zsumRR) // s^2 - z*sumRR

		// 3) calculer biBig = - z*sumRE
		zsumRE := new(big.Float).SetPrec(bigPrec).Mul(zBig, sumRE) // z * sumRE
		biBig := new(big.Float).SetPrec(bigPrec).Neg(zsumRE)       // - (z * sumRE)

		// 4) calculer diBig = s^2 - z*sumEE
		zsumEE := new(big.Float).SetPrec(bigPrec).Mul(zBig, sumEE)  // z * sumEE
		diBig := new(big.Float).SetPrec(bigPrec).Sub(s2Big, zsumEE) // s^2 - z*sumEE

		// 5) vérifier les trois identités exactes (à tolérance tolBig)
		//    a) aiBig + z*sumRR == s^2  ==>  (aiBig + zsumRR) - s^2 == 0
		lhsA := new(big.Float).SetPrec(bigPrec).Add(aiBig, zsumRR) // aiBig + z*sumRR
		diffA := new(big.Float).SetPrec(bigPrec).Sub(lhsA, s2Big)  // (aiBig + zsumRR) - s^2

		//    b) biBig + z*sumRE == 0     ==>  (biBig + zsumRE) - 0 == 0
		lhsB := new(big.Float).SetPrec(bigPrec).Add(biBig, zsumRE) // biBig + z*sumRE
		diffB := new(big.Float).SetPrec(bigPrec).Set(lhsB)         // (biBig + zsumRE)

		//    c) diBig + z*sumEE == s^2  ==>  (diBig + zsumEE) - s^2 == 0
		lhsD := new(big.Float).SetPrec(bigPrec).Add(diBig, zsumEE) // diBig + z*sumEE
		diffD := new(big.Float).SetPrec(bigPrec).Sub(lhsD, s2Big)  // (diBig + zsumEE) - s^2

		// Absolu des écarts
		absA := new(big.Float).SetPrec(bigPrec).Abs(diffA)
		absB := new(big.Float).SetPrec(bigPrec).Abs(diffB)
		absD := new(big.Float).SetPrec(bigPrec).Abs(diffD)

		// Comparer à la tolérance tolBig
		if absA.Cmp(tolBig) > 0 || absB.Cmp(tolBig) > 0 || absD.Cmp(tolBig) > 0 {
			// Convertir en float64 pour le log (uniquement à des fins d’affichage)
			aF, _ := aiBig.Float64()
			bF, _ := biBig.Float64()
			dF, _ := diBig.Float64()
			log.Panicf(
				"SamplePz sanity‐check failed au slot %d : a_i=%.12g, b_i=%.12g, d_i=%.12g",
				i, aF, bF, dF,
			)
		}
	}

	// Si on arrive ici, tout est bon :
	fmt.Printf("SamplePz: sanity‐check passed\n")
	// ! —––– end sanity‐check —––––

	// 4) Sample Q ∈ ℤ^{k×n} from D_{√(s²−α²)} and CRT→NTT→qhat
	sigmaQ := math.Sqrt(s*s - alpha*alpha)
	fmt.Printf("SamplePz: σQ = %f\n", sigmaQ)
	qhat := make([]*ring.Poly, k)
	dgQ := NewDiscreteGaussian(sigmaQ)
	for j := 0; j < k; j++ {
		// draw a length‐n integer vector
		ints := make([]int64, n)
		for i := range ints {
			ints[i] = dgQ.Draw(0.0) // mean = 0
		}
		// encode & NTT
		P := ringQ.NewPoly()
		for lvl, qi := range ringQ.Modulus {
			mod := int64(qi)
			for i := 0; i < n; i++ {
				c := ints[i] % mod
				if c < 0 {
					c += mod
				}
				P.Coeffs[lvl][i] = uint64(c)
			}
		}
		ringQ.NTT(P, P)
		qhat[j] = P
	}
	//! --- SamplePz sanity‐check: empirical variance of the qhat samples ---
	var sum, sumSq float64
	count := 0
	for _, Phat := range qhat {
		// invert NTT → get level-0 coeffs
		coeff := ringQ.NewPoly()
		ringQ.InvNTT(Phat, coeff)
		for i := 0; i < n; i++ {
			v := int64(coeff.Coeffs[0][i])
			// center into [−q/2, q/2)
			if v > int64(ringQ.Modulus[0]/2) {
				v -= int64(ringQ.Modulus[0])
			}
			x := float64(v)
			sum += x
			sumSq += x * x
			count++
		}
	}
	mean := sum / float64(count)
	variance := sumSq/float64(count) - mean*mean
	fmt.Printf("SamplePz → q̂ empirical variance = %f  (target ≈ %f)\n",
		variance, sigmaQ*sigmaQ)
	//! --- end variance check ---

	// 5) Build centers c0, c1 ∈ K_{2n} (coeff‐domain) via CRT inner‐products
	scaleN := new(big.Float).SetPrec(prec).SetFloat64(-alpha * alpha)
	scaleD := new(big.Float).SetPrec(prec).Sub(s2, α2)
	scale := new(big.Float).SetPrec(prec).Quo(scaleN, scaleD)
	scaleC := NewBigComplexFromFloat(scale, big.NewFloat(0).SetPrec(prec))

	c0 := NewFieldElemBig(n, prec)
	c1 := NewFieldElemBig(n, prec)

	for j := 0; j < k; j++ {
		qPoly := ringQ.NewPoly()
		ringQ.InvNTT(qhat[j], qPoly)
		rPoly := ringQ.NewPoly()
		ringQ.InvNTT(Ttilde[0][j], rPoly)
		ePoly := ringQ.NewPoly()
		ringQ.InvNTT(Ttilde[1][j], ePoly)
		qF := NegacyclicEvaluatePoly(qPoly, ringQ, prec)
		rF := NegacyclicEvaluatePoly(rPoly, ringQ, prec)
		eF := NegacyclicEvaluatePoly(ePoly, ringQ, prec)
		rT := HermitianTransposeFieldElem(rF)
		eT := HermitianTransposeFieldElem(eF)
		for i := 0; i < n; i++ {
			c0.Coeffs[i] = c0.Coeffs[i].Add(rT.Coeffs[i].Mul(qF.Coeffs[i]))
			c1.Coeffs[i] = c1.Coeffs[i].Add(eT.Coeffs[i].Mul(qF.Coeffs[i]))
		}
	}
	c0 = FieldScalarMulBig(c0, scaleC)
	c1 = FieldScalarMulBig(c1, scaleC)
	c0.Domain = Eval
	c1.Domain = Eval
	// c0, c1 are now in EVAL domain, but we need them in COEFF
	c0 = FloatToCoeffNegacyclic(c0, prec)
	c1 = FloatToCoeffNegacyclic(c1, prec)
	fmt.Printf("SamplePz: c0 = %s, c1 = %s\n", c0.Coeffs[0].Real.Text('g', 10), c1.Coeffs[0].Real.Text('g', 10))

	// 6) Final 2×2 block‐recursive sample → p0,p1
	p0, p1 := Sample2zField(
		aFld, bFld, dFld,
		c0, c1,
		ringQ.N, ringQ.Modulus[0], prec,
	)
	fmt.Printf("SamplePz: p0 = %s, p1 = %s\n", p0.Coeffs[0].Real.Text('g', 10), p1.Coeffs[0].Real.Text('g', 10))
	// 7) Convert p0,p1 back to NTT‐domain polys
	out := make([]*ring.Poly, expectedLength)
	P0 := NegacyclicInterpolateElem(p0, ringQ)
	ringQ.NTT(P0, P0)
	P1 := NegacyclicInterpolateElem(p1, ringQ)
	ringQ.NTT(P1, P1)
	out[0], out[1] = P0, P1

	// 8) Append the k “q̂” polys
	for j := 0; j < k; j++ {
		out[j+2] = qhat[j]
	}
	return out
}

// makeSmallRing returns a *ring.Ring of degree N (a power of two) and modulus q.
func makeSmallRing(N int, q uint64) *ring.Ring {
	r, err := ring.NewRing(N, []uint64{q})
	if err != nil {
		panic(err)
	}
	return r
}
