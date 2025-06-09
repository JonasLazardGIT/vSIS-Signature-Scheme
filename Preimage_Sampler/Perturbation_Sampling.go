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
	for i := 0; i < n; i++ {
	}
	q1 = SampleFZBig(dCoeff, c1, modulus, prec)

	// 2) delta = q1 – c1  (both in COEFF), then lift to EVAL
	delta := FieldSubBig(q1, c1)

	delta = FloatToEvalNegacyclic(delta, prec) // delta is now in EVAL

	invD, norms := FieldInverseDiagWithNorm(d)
	invD = FieldScalarDiv(invD, norms)
	invD.Domain = Eval

	tmp := FieldMulBig(invD, delta)
	tmp.Domain = Eval
	scaledDelta := FieldMulBig(b, tmp)

	// scaledDelta.Coeffs is still Eval‐domain; fix the flag:
	scaledDelta.Domain = Eval
	scaledDelta = FloatToCoeffNegacyclic(scaledDelta, prec)

	// 5) c0′ = c0 + scaledDelta   (all COEFF)
	c0p := FieldAddBig(c0, scaledDelta)

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
	// fmt.Printf("Sample2zField: aPr = %s\n", aPr.Coeffs[0].Real.Text('g', 10))
	// fmt.Printf("Sample2zField: c0′ = %s\n", c0p.Coeffs[0].Real.Text('g', 10))

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

	fmt.Printf("SamplePz: q = %d, s = %f, α = %f, z = %s\n", ringQ.Modulus[0], s, alpha, zBig.Text('g', 10))

	// 2) Accumulate sums in evaluation domain
	va := NewFieldElemBig(n, prec)
	vb := NewFieldElemBig(n, prec)
	vd := NewFieldElemBig(n, prec)
	va.Domain = Eval
	vb.Domain = Eval
	vd.Domain = Eval

	for j := 0; j < k; j++ {
		// convert row polynomials back to coefficient form
		rPoly := ringQ.NewPoly()
		ringQ.InvNTT(Ttilde[0][j], rPoly)
		ePoly := ringQ.NewPoly()
		ringQ.InvNTT(Ttilde[1][j], ePoly)

		// evaluate as complex vectors (Eval domain)
		rF := NegacyclicEvaluatePoly(rPoly, ringQ, prec)
		eF := NegacyclicEvaluatePoly(ePoly, ringQ, prec)

		rT := HermitianTransposeFieldElem(rF)
		eT := HermitianTransposeFieldElem(eF)

		for i := 0; i < n; i++ {
			va.Coeffs[i] = va.Coeffs[i].Add(rT.Coeffs[i].Mul(rF.Coeffs[i]))
			vb.Coeffs[i] = vb.Coeffs[i].Add(eF.Coeffs[i].Mul(rT.Coeffs[i]))
			vd.Coeffs[i] = vd.Coeffs[i].Add(eT.Coeffs[i].Mul(eF.Coeffs[i]))
		}
	}
	// 3) Switch the sums to coefficient domain
	va = FloatToCoeffNegacyclic(va, prec)
	vb = FloatToCoeffNegacyclic(vb, prec)
	vd = FloatToCoeffNegacyclic(vd, prec)

	// scalarFactor = -z
	scalarFactor := new(big.Float).Copy(zBig)
	scalarC := NewBigComplexFromFloat(scalarFactor, big.NewFloat(0).SetPrec(prec))

	// build a,b,d in coefficient domain
	aFld := FieldScalarMulBig(va, scalarC)
	bFld := FieldScalarMulBig(vb, scalarC)
	dFld := FieldScalarMulBig(vd, scalarC)

	s2c := NewBigComplexFromFloat(s2, big.NewFloat(0).SetPrec(prec))
	aFld.Coeffs[0] = aFld.Coeffs[0].Add(s2c)
	dFld.Coeffs[0] = dFld.Coeffs[0].Add(s2c)

	// back to evaluation domain for further sampling
	aFld = FloatToEvalNegacyclic(aFld, prec)
	bFld = FloatToEvalNegacyclic(bFld, prec)
	dFld = FloatToEvalNegacyclic(dFld, prec)
	//!------------------------------------------------------------------
	//! 3-bis)  ‖α · [Tᵗ | I]‖  ≤  s   ︙  spectral-norm sanity check
	//!------------------------------------------------------------------
	{
		eps := 1e-12 // relative slack
		maxNorm := 0.0
		for i := 0; i < n; i++ {
			rr, _ := va.Coeffs[i].Real.Float64() // Σ r̂_i²  (real part)
			re, _ := vb.Coeffs[i].Real.Float64() // Σ r̂_i ê̂_i
			ee, _ := vd.Coeffs[i].Real.Float64() // Σ ê̂_i²
			rr += 1.0                            // + identity contribution
			ee += 1.0

			trace := rr + ee
			discr := trace*trace - 4*(rr*ee-re*re)
			if discr < 0 { // numerical guard
				discr = 0
			}
			lambdaMax := 0.5 * (trace + math.Sqrt(discr))
			spec := alpha * math.Sqrt(lambdaMax) // α‖·‖₂ for slot i
			if spec > maxNorm {
				maxNorm = spec
			}
		}
		if maxNorm > s*(1+eps) {
			log.Panicf("parameter s = %.6g is too small; need ≥ %.6g to satisfy ‖α[Tᵗ|I]‖≤s",
				s, maxNorm)
		} else {
			fmt.Printf("spectral-norm check: α‖[Tᵗ|I]‖ = %.6g  ≤  s = %.6g  ✔\n",
				maxNorm, s)
		}
	}
	//! ------------------------------------------------------------------
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
				c := SignedToUnsigned(ints[i], uint64(mod))
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
	standard_deviation := math.Sqrt(variance)
	fmt.Printf("SamplePz → q̂ empirical standard deviation = %f\n", standard_deviation)
	//! --- end variance check ---

	// ---------------------------
	//  5. Build centres in R_q
	//     c0 = –z · Σ_j r̂_j · Q_j
	//     c1 = –z · Σ_j ê̂_j · Q_j
	//
	// ---------------------------
	c0Poly := ringQ.NewPoly() // coefficient‐domain accumulator
	c1Poly := ringQ.NewPoly()
	rDotQEval := ringQ.NewPoly()
	eDotQEval := ringQ.NewPoly()
	tmp := ringQ.NewPoly()

	// (a) accumulate inner‐products in evaluation domain
	for j := 0; j < k; j++ {
		ringQ.MulCoeffsMontgomery(Ttilde[0][j], qhat[j], tmp)
		ringQ.Add(rDotQEval, tmp, rDotQEval)

		ringQ.MulCoeffsMontgomery(Ttilde[1][j], qhat[j], tmp)
		ringQ.Add(eDotQEval, tmp, eDotQEval)
	}

	// convert to coefficient domain
	ringQ.InvNTT(rDotQEval, c0Poly)
	ringQ.InvNTT(eDotQEval, c1Poly)

	// debug: show inner products before scaling
	centre := func(v uint64) int64 {
		q := int64(ringQ.Modulus[0])
		x := int64(v)
		if x > q/2 {
			x -= q
		}
		return x
	}
	r0 := centre(c0Poly.Coeffs[0][0])
	e0 := centre(c1Poly.Coeffs[0][0])
	fmt.Printf("SamplePz: inner r·Q(0)=%d, e·Q(0)=%d\n", r0, e0)

	// (b) scale each coefficient by –z and round to nearest integer mod q
	zF, _ := zBig.Float64()
	for i := 0; i < ringQ.N; i++ {
		// center into (−q/2,q/2]
		v0 := int64(c0Poly.Coeffs[0][i])
		if v0 > int64(ringQ.Modulus[0]/2) {
			v0 -= int64(ringQ.Modulus[0])
		}
		// multiply by –z and round
		s0 := int64(math.Round(-zF * float64(v0)))
		c0Poly.Coeffs[0][i] = SignedToUnsigned(s0, ringQ.Modulus[0])
		// verify modular reduction is exact
		expect := s0 % int64(ringQ.Modulus[0])
		if expect > int64(ringQ.Modulus[0])/2 {
			expect -= int64(ringQ.Modulus[0])
		}
		if expect < -int64(ringQ.Modulus[0])/2 {
			expect += int64(ringQ.Modulus[0])
		}
		if UnsignedToSigned(c0Poly.Coeffs[0][i], ringQ.Modulus[0]) != expect {
			log.Fatalf("rounding check failed for c0 coeff %d", i)
		}

		v1 := int64(c1Poly.Coeffs[0][i])
		if v1 > int64(ringQ.Modulus[0]/2) {
			v1 -= int64(ringQ.Modulus[0])
		}
		s1 := int64(math.Round(-zF * float64(v1)))
		c1Poly.Coeffs[0][i] = SignedToUnsigned(s1, ringQ.Modulus[0])
		expect1 := s1 % int64(ringQ.Modulus[0])
		if expect1 > int64(ringQ.Modulus[0])/2 {
			expect1 -= int64(ringQ.Modulus[0])
		}
		if expect1 < -int64(ringQ.Modulus[0])/2 {
			expect1 += int64(ringQ.Modulus[0])
		}
		if UnsignedToSigned(c1Poly.Coeffs[0][i], ringQ.Modulus[0]) != expect1 {
			log.Fatalf("rounding check failed for c1 coeff %d", i)
		}
	}

	c0Eval := ConvertFromPolyBig(ringQ, c0Poly, prec)
	c0Coeff := ToCoeffNegacyclic(c0Eval, ringQ, prec)
	c1Eval := ConvertFromPolyBig(ringQ, c1Poly, prec)
	c1Coeff := ToCoeffNegacyclic(c1Eval, ringQ, prec)

	fmt.Printf("SamplePz: c0 = %s, c1 = %s\n",
		c0Coeff.Coeffs[0].Real.Text('g', 10),
		c1Coeff.Coeffs[0].Real.Text('g', 10))

	// 6) Final 2×2 sampler
	p0, p1 := Sample2zField(
		aFld, bFld, dFld,
		c0Coeff, c1Coeff,
		ringQ.N, ringQ.Modulus[0], prec,
	)

	// fmt.Printf("SamplePz: p0 = %s, p1 = %s\n", p0.Coeffs[0].Real.Text('g', 10), p1.Coeffs[0].Real.Text('g', 10))
	// 7) Convert p0,p1 back to NTT‐domain polys
	out := make([]*ring.Poly, expectedLength)
	P0 := NegacyclicInterpolateElem(p0, ringQ)
	P1 := NegacyclicInterpolateElem(p1, ringQ)

	//! ------------------------------------------------------------------
	//! DEBUG -- sign diagnosis for p0 / p1  vs  (ê̂ᵀ·Q) / (r̂ᵀ·Q)
	//! ------------------------------------------------------------------
	{
		// helper to centre a coeff into (−q/2, q/2]
		centre := func(v uint64) int64 {
			q := int64(ringQ.Modulus[0])
			x := int64(v)
			if x > q/2 {
				x -= q
			}
			return x
		}

		// ❶  compute ê̂ᵀ·Q   and   r̂ᵀ·Q   in NTT then go back to COEFF
		eDotQEval := ringQ.NewPoly() // zero
		rDotQEval := ringQ.NewPoly() // zero
		tmp := ringQ.NewPoly()
		for j := 0; j < k; j++ {
			// eHat = Ttilde[1][j] ,  rHat = Ttilde[0][j]   (both NTT)
			ringQ.MulCoeffsMontgomery(Ttilde[1][j], qhat[j], tmp)
			ringQ.Add(eDotQEval, tmp, eDotQEval)

			ringQ.MulCoeffsMontgomery(Ttilde[0][j], qhat[j], tmp)
			ringQ.Add(rDotQEval, tmp, rDotQEval)
		}
		eDotQCoeff := ringQ.NewPoly()
		rDotQCoeff := ringQ.NewPoly()
		ringQ.InvNTT(eDotQEval, eDotQCoeff)
		ringQ.InvNTT(rDotQEval, rDotQCoeff)

		// ❷  pick the *constant* coefficient 0  for an easy sign test
		p0c0 := centre(P0.Coeffs[0][0])         // p₀(0)   (Coeff domain)
		p1c0 := centre(P1.Coeffs[0][0])         // p₁(0)
		eQc0 := centre(eDotQCoeff.Coeffs[0][0]) // (ê̂ᵀ·Q)(0)
		rQc0 := centre(rDotQCoeff.Coeffs[0][0]) // (r̂ᵀ·Q)(0)

		fmt.Printf("[SIGN-DIAG]  p0(0)=%d   r·Q(0)=%d   sum=%d\n", p0c0, rQc0, p0c0+rQc0)
		fmt.Printf("[SIGN-DIAG]  p1(0)=%d   e·Q(0)=%d   sum=%d\n", p1c0, eQc0, p1c0+eQc0)

	}

	//! ------------------------------------------------------------------

	ringQ.NTT(P0, P0)
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
