package Preimage_Sampler

import (
	"io"
	"math"
	"math/big"
	"math/rand"
	"os"
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

// silenceStdout runs f while discarding anything written to stdout.
func silenceStdout(f func()) {
	old := os.Stdout
	r, w, _ := os.Pipe()
	os.Stdout = w
	f()
	w.Close()
	os.Stdout = old
	io.ReadAll(r)
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

	scalar := zBig
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

	var out []*ring.Poly
	silenceStdout(func() {
		out = SamplePz(ringQ, s, alpha, T, 3, prec)
	})
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

	s2 := s * s
	a2 := alpha * alpha

	aExp := make([]int64, n)
	bExp := make([]int64, n)
	dExp := make([]int64, n)
	for i := 0; i < n; i++ {
		aExp[i] = int64(math.Round(-a2 * float64(rTr[i])))
		dExp[i] = int64(math.Round(-a2 * float64(eTe[i])))
		bExp[i] = int64(math.Round(-a2 * float64(rTe[i])))
	}
	aExp[0] += int64(math.Round(s2))
	dExp[0] += int64(math.Round(s2))

	for i := 0; i < n; i++ {
		if math.Abs(float64(aInts[i]-aExp[i])) > 20 {
			t.Fatalf("a[%d]=%d want approx %d", i, aInts[i], aExp[i])
		}
		if math.Abs(float64(bInts[i]-bExp[i])) > 20 {
			t.Fatalf("b[%d]=%d want approx %d", i, bInts[i], bExp[i])
		}
		if math.Abs(float64(dInts[i]-dExp[i])) > 20 {
			t.Fatalf("d[%d]=%d want approx %d", i, dInts[i], dExp[i])
		}
	}
}

func TestPerturbationVariance(t *testing.T) {
	const (
		n    = 16
		q    = uint64(97)
		prec = uint(256)
		base = uint64(2)
	)
	ringQ := makeSmallRing(n, q)

	k := 4
	alpha := 2.25
	s := 12.0
	trials := 4000
	slack := 0.05

	SetForceZero(true)
	trap := TrapGen(ringQ, base, alpha)
	SetForceZero(false)
	T := [2][]*ring.Poly{trap.R[0][:k], trap.R[1][:k]}

	var sumSq float64
	count := 0
	for t := 0; t < trials; t++ {
		var out []*ring.Poly
		silenceStdout(func() {
			out = SamplePz(ringQ, s, alpha, T, k+2, prec)
		})
		for j := 0; j < k; j++ {
			coeff := ringQ.NewPoly()
			ringQ.InvNTT(out[j+2], coeff)
			for i := 0; i < n; i++ {
				v := UnsignedToSigned(coeff.Coeffs[0][i], q)
				sumSq += float64(v * v)
				count++
			}
		}
	}
	variance := sumSq / float64(count)
	want := s*s - alpha*alpha
	if variance < (1-slack)*want || variance > (1+slack)*want {
		t.Fatalf("variance %.4f outside ±%.0f%% of want %.4f", variance, slack*100, want)
	}
}

func TestPerturbationABDNonZero(t *testing.T) {
	const (
		n    = 16
		q    = uint64(97)
		prec = uint(256)
	)
	ringQ := makeSmallRing(n, q)

	rCoeffs := make([]int64, n)
	eCoeffs := make([]int64, n)
	rCoeffs[0] = 1
	rCoeffs[1] = 1
	rCoeffs[2] = 2
	eCoeffs[0] = -3
	eCoeffs[5] = 4

	k := 3
	s := 12.0
	alpha := 2.25

	rPolys := make([]*ring.Poly, k)
	ePolys := make([]*ring.Poly, k)
	for j := 0; j < k; j++ {
		rPolys[j] = ringQ.NewPoly()
		ePolys[j] = ringQ.NewPoly()
		for i := 0; i < n; i++ {
			rPolys[j].Coeffs[0][i] = SignedToUnsigned(rCoeffs[i], q)
			ePolys[j].Coeffs[0][i] = SignedToUnsigned(eCoeffs[i], q)
		}
		ringQ.NTT(rPolys[j], rPolys[j])
		ringQ.NTT(ePolys[j], ePolys[j])
	}

	T := [2][]*ring.Poly{rPolys, ePolys}

	aCalc, bCalc, dCalc := computeABDUsingCode(ringQ, s, alpha, T, prec)

	aInts := coeffsToIntSlice(aCalc)
	bInts := coeffsToIntSlice(bCalc)
	dInts := coeffsToIntSlice(dCalc)

	rT := transposeNegacyclic(rCoeffs)
	eT := transposeNegacyclic(eCoeffs)
	rTr := negacyclicConvolution(rT, rCoeffs)
	rTe := negacyclicConvolution(rT, eCoeffs)
	eTe := negacyclicConvolution(eT, eCoeffs)

	for i := 0; i < n; i++ {
		rTr[i] *= int64(k)
		rTe[i] *= int64(k)
		eTe[i] *= int64(k)
	}

	s2 := s * s
	a2 := alpha * alpha

	aExp := make([]int64, n)
	bExp := make([]int64, n)
	dExp := make([]int64, n)
	for i := 0; i < n; i++ {
		aExp[i] = int64(math.Round(-a2 * float64(rTr[i])))
		dExp[i] = int64(math.Round(-a2 * float64(eTe[i])))
		bExp[i] = int64(math.Round(-a2 * float64(rTe[i])))
	}
	aExp[0] += int64(math.Round(s2))
	dExp[0] += int64(math.Round(s2))

	for i := 0; i < n; i++ {
		if math.Abs(float64(aInts[i]-aExp[i])) > 20 {
			t.Fatalf("a[%d]=%d want approx %d", i, aInts[i], aExp[i])
		}
		if math.Abs(float64(bInts[i]-bExp[i])) > 20 {
			t.Fatalf("b[%d]=%d want approx %d", i, bInts[i], bExp[i])
		}
		if math.Abs(float64(dInts[i]-dExp[i])) > 20 {
			t.Fatalf("d[%d]=%d want approx %d", i, dInts[i], dExp[i])
		}
	}
}

func TestRoundTripAndLargeParams(t *testing.T) {
	const (
		n    = 16
		q    = uint64(97)
		prec = uint(256)
	)
	ringQ := makeSmallRing(n, q)

	randCoeff := func() int64 {
		return int64(rand.Intn(10001)) - 5000
	}

	orig := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		v := randCoeff()
		orig.Coeffs[i] = NewBigComplex(float64(v), 0, prec)
	}
	orig.Domain = Coeff

	eval := FloatToEvalNegacyclic(orig, prec)
	back := FloatToCoeffNegacyclic(eval, prec)

	tol := new(big.Float).SetPrec(prec).SetFloat64(1e-9)
	for i := 0; i < n; i++ {
		diffR := new(big.Float).Sub(orig.Coeffs[i].Real, back.Coeffs[i].Real)
		if diffR.Abs(diffR).Cmp(tol) > 0 {
			t.Fatalf("real mismatch at %d", i)
		}
		diffI := new(big.Float).Sub(orig.Coeffs[i].Imag, back.Coeffs[i].Imag)
		if diffI.Abs(diffI).Cmp(tol) > 0 {
			t.Fatalf("imag mismatch at %d", i)
		}
	}
	_ = ringQ
}
