package Preimage_Sampler

import (
	"math"
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// transposeNegacyclic returns the negacyclic transpose of a polynomial.
func transposeNegacyclic(coeffs []int64) []int64 {
	n := len(coeffs)
	out := make([]int64, n)
	out[0] = coeffs[0]
	for i := 1; i < n; i++ {
		out[i] = -coeffs[n-i]
	}
	return out
}

// negacyclicConvolution multiplies a and b modulo x^n+1 over the integers.
func negacyclicConvolution(a, b []int64) []int64 {
	n := len(a)
	res := make([]int64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			k := i + j
			prod := a[i] * b[j]
			if k >= n {
				k -= n
				prod = -prod
			}
			res[k] += prod
		}
	}
	return res
}

// coeffsToIntSlice converts the real part of f to int64 coefficients.
func coeffsToIntSlice(f *CyclotomicFieldElem) []int64 {
	n := f.N
	out := make([]int64, n)
	for i := 0; i < n; i++ {
		val, _ := f.Coeffs[i].Real.Float64()
		out[i] = int64(math.Round(val))
	}
	return out
}

// computeABDUsingCode replicates the algebraic part of SamplePz.
func computeABDUsingCode(ringQ *ring.Ring, s, alpha float64, T [2][]*ring.Poly, prec uint) (a, b, d *CyclotomicFieldElem) {
	n := ringQ.N
	k := len(T[0])

	one := new(big.Float).SetPrec(prec).SetFloat64(1)
	s2 := new(big.Float).SetPrec(prec).SetFloat64(s * s)
	a2 := new(big.Float).SetPrec(prec).SetFloat64(alpha * alpha)
	invS2 := new(big.Float).SetPrec(prec).Quo(one, s2)
	invA2 := new(big.Float).SetPrec(prec).Quo(one, a2)
	diff := new(big.Float).SetPrec(prec).Sub(invA2, invS2)
	zBig := new(big.Float).SetPrec(prec).Quo(one, diff)

	va := NewFieldElemBig(n, prec)
	vb := NewFieldElemBig(n, prec)
	vd := NewFieldElemBig(n, prec)
	va.Domain, vb.Domain, vd.Domain = Eval, Eval, Eval

	for j := 0; j < k; j++ {
		rPoly := ringQ.NewPoly()
		ringQ.InvNTT(T[0][j], rPoly)
		ePoly := ringQ.NewPoly()
		ringQ.InvNTT(T[1][j], ePoly)
		rF := NegacyclicEvaluatePoly(rPoly, ringQ, prec)
		eF := NegacyclicEvaluatePoly(ePoly, ringQ, prec)
		rT := HermitianTransposeFieldElem(rF)
		eT := HermitianTransposeFieldElem(eF)
		va = FieldAddBig(va, FieldMulBig(rT, rF))
		vb = FieldAddBig(vb, FieldMulBig(eF, rT))
		vd = FieldAddBig(vd, FieldMulBig(eT, eF))
		va.Domain, vb.Domain, vd.Domain = Eval, Eval, Eval
	}

	va = FloatToCoeffNegacyclic(va, prec)
	vb = FloatToCoeffNegacyclic(vb, prec)
	vd = FloatToCoeffNegacyclic(vd, prec)

	scalar := new(big.Float).Neg(zBig)
	scalarC := NewBigComplexFromFloat(scalar, new(big.Float).SetPrec(prec).SetFloat64(0))

	a = FieldScalarMulBig(va, scalarC)
	b = FieldScalarMulBig(vb, scalarC)
	d = FieldScalarMulBig(vd, scalarC)

	s2c := NewBigComplexFromFloat(s2, new(big.Float).SetPrec(prec).SetFloat64(0))
	a.Coeffs[0] = a.Coeffs[0].Add(s2c)
	d.Coeffs[0] = d.Coeffs[0].Add(s2c)
	a.Domain, b.Domain, d.Domain = Coeff, Coeff, Coeff
	return
}

func TestPerturbationSamplingDeterministic(t *testing.T) {
	const (
		n    = 16
		q    = uint64(97)
		prec = uint(256)
	)
	ringQ := makeSmallRing(n, q)

	rCoeffs := []int64{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
	eCoeffs := []int64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

	rPoly := ringQ.NewPoly()
	ePoly := ringQ.NewPoly()
	for i := 0; i < n; i++ {
		rPoly.Coeffs[0][i] = SignedToUnsigned(rCoeffs[i], q)
		ePoly.Coeffs[0][i] = SignedToUnsigned(eCoeffs[i], q)
	}
	ringQ.NTT(rPoly, rPoly)
	ringQ.NTT(ePoly, ePoly)

	T := [2][]*ring.Poly{{rPoly}, {ePoly}}

	s := 12.0
	alpha := 2.0

	SetForceZero(true)
	defer SetForceZero(false)

	out := SamplePz(ringQ, s, alpha, T, 3, prec)
	if len(out) != 3 {
		t.Fatalf("expected 3 polys, got %d", len(out))
	}
	// All outputs should be zero polynomials
	for idx, p := range out {
		coeff := ringQ.NewPoly()
		ringQ.InvNTT(p, coeff)
		for i := 0; i < n; i++ {
			if coeff.Coeffs[0][i] != 0 {
				t.Fatalf("output %d coeff %d = %d, want 0", idx, i, coeff.Coeffs[0][i])
			}
		}
	}

	aCalc, bCalc, dCalc := computeABDUsingCode(ringQ, s, alpha, T, prec)

	aInts := coeffsToIntSlice(aCalc)
	bInts := coeffsToIntSlice(bCalc)
	dInts := coeffsToIntSlice(dCalc)

	rT := transposeNegacyclic(rCoeffs)
	eT := transposeNegacyclic(eCoeffs)
	rTr := negacyclicConvolution(rT, rCoeffs)
	rTe := negacyclicConvolution(rT, eCoeffs)
	eTe := negacyclicConvolution(eT, eCoeffs)

	s2 := int64(s * s)
	a2 := int64(alpha * alpha)

	aExp := make([]int64, n)
	bExp := make([]int64, n)
	dExp := make([]int64, n)
	for i := 0; i < n; i++ {
		aExp[i] = -a2 * rTr[i]
		dExp[i] = -a2 * eTe[i]
		bExp[i] = -a2 * rTe[i]
	}
	aExp[0] += s2
	dExp[0] += s2

	for i := 0; i < n; i++ {
		if aInts[i] != aExp[i] {
			t.Fatalf("a[%d]=%d want %d", i, aInts[i], aExp[i])
		}
		if bInts[i] != bExp[i] {
			t.Fatalf("b[%d]=%d want %d", i, bInts[i], bExp[i])
		}
		if dInts[i] != dExp[i] {
			t.Fatalf("d[%d]=%d want %d", i, dInts[i], dExp[i])
		}
	}
}
