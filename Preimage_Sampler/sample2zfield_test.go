package Preimage_Sampler

import "testing"

func TestSample2zFieldDebug(t *testing.T) {
	const (
		n    = 1
		q    = uint64(97)
		prec = uint(64)
	)

	aCoeff := NewFieldElemBig(n, prec)
	bCoeff := NewFieldElemBig(n, prec)
	dCoeff := NewFieldElemBig(n, prec)
	c0 := NewFieldElemBig(n, prec)
	c1 := NewFieldElemBig(n, prec)

	for i := 0; i < n; i++ {
		aCoeff.Coeffs[i] = NewBigComplex(2, 0, prec)
		bCoeff.Coeffs[i] = NewBigComplex(1, 0, prec)
		dCoeff.Coeffs[i] = NewBigComplex(1, 0, prec)
		c0.Coeffs[i] = NewBigComplex(0, 0, prec)
		c1.Coeffs[i] = NewBigComplex(0, 0, prec)
	}
	a := FloatToEvalNegacyclic(aCoeff, prec)
	b := FloatToEvalNegacyclic(bCoeff, prec)
	d := FloatToEvalNegacyclic(dCoeff, prec)
	c0.Domain, c1.Domain = Coeff, Coeff

	SetForceZero(true)
	defer SetForceZero(false)

	for i := 0; i < 100; i++ {
		q0, q1 := Sample2zField(a, b, d, c0, c1, n, q, prec)
		if q0 == nil || q1 == nil {
			t.Fatalf("nil output on trial %d", i)
		}
	}
}
