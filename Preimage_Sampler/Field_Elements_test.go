// bigcomplex_test.go
package Preimage_Sampler

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// ---------------------------------------------------------------------------------------------------------------------
// 1) Helpers
// ---------------------------------------------------------------------------------------------------------------------

// approxEqual reports whether |a-b| ≤ tol.
func approxEqual(a, b, tol float64) bool {
	return math.Abs(a-b) <= tol
}

// randomUintPoly fills a *ring.Poly with random coefficients in [0, q-1].
func randomUintPoly(ringQ *ring.Ring) *ring.Poly {
	P := ringQ.NewPoly()
	mod := ringQ.Modulus[0]
	for i := 0; i < ringQ.N; i++ {
		P.Coeffs[0][i] = uint64(rand.Intn(int(mod)))
	}
	return P
}

// ---------------------------------------------------------------------------------------------------------------------
// 2) BigComplex arithmetic tests
// ---------------------------------------------------------------------------------------------------------------------

func TestBigComplexBasicOperations(t *testing.T) {
	const prec = uint(64)
	// z = (3 + 4i), w = (-1 + 2i)
	z := NewBigComplex(3.0, 4.0, prec)
	w := NewBigComplex(-1.0, 2.0, prec)

	// Add: z + w = (2 + 6i)
	sum := z.Add(w)
	r1, _ := sum.Real.Float64()
	i1, _ := sum.Imag.Float64()
	if !approxEqual(r1, 2.0, 1e-15) || !approxEqual(i1, 6.0, 1e-15) {
		t.Fatalf("Add: got %v + %vi, want 2 + 6i", r1, i1)
	}

	// Sub: z - w = (4 + 2i)
	diff := z.Sub(w)
	r2, _ := diff.Real.Float64()
	i2, _ := diff.Imag.Float64()
	if !approxEqual(r2, 4.0, 1e-15) || !approxEqual(i2, 2.0, 1e-15) {
		t.Fatalf("Sub: got %v + %vi, want 4 + 2i", r2, i2)
	}

	// Mul: z * w = (3·-1 - 4·2) + (3·2 + 4·-1)i = (-11 + 2i)
	prod := z.Mul(w)
	r3, _ := prod.Real.Float64()
	i3, _ := prod.Imag.Float64()
	if !approxEqual(r3, -11.0, 1e-15) || !approxEqual(i3, 2.0, 1e-15) {
		t.Fatalf("Mul: got %v + %vi, want -11 + 2i", r3, i3)
	}

	// Conj: conj(z) = (3 - 4i)
	conj := z.Conj()
	r4, _ := conj.Real.Float64()
	i4, _ := conj.Imag.Float64()
	if !approxEqual(r4, 3.0, 1e-15) || !approxEqual(i4, -4.0, 1e-15) {
		t.Fatalf("Conj: got %v + %vi, want 3 - 4i", r4, i4)
	}

	// AbsSquared: |z|^2 = 3^2 + 4^2 = 25
	mag2 := z.AbsSquared()
	m2, _ := mag2.Float64()
	if !approxEqual(m2, 25.0, 1e-15) {
		t.Fatalf("AbsSquared: got %v, want 25", m2)
	}

	// Inv: 1/(3+4i) = (3 - 4i)/(3^2+4^2) = (3/25) - (4/25)i
	inv := z.Inv()
	r5, _ := inv.Real.Float64()
	i5, _ := inv.Imag.Float64()
	if !approxEqual(r5, 3.0/25.0, 1e-15) || !approxEqual(i5, -4.0/25.0, 1e-15) {
		t.Fatalf("Inv: got %v + %vi, want 0.12 - 0.16i", r5, i5)
	}

	// Copy: modifying copy does not affect original
	cp := z.Copy()
	// change cp’s real part to zero:
	cp.Real.SetFloat64(0.0)
	r6, _ := z.Real.Float64()
	if !approxEqual(r6, 3.0, 1e-15) {
		t.Fatalf("Copy: modifying copy changed original; original Real = %v", r6)
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 3) CyclotomicFieldElem ↔ ring.Poly (FFTBig/IFFTBig) round-trip
// ---------------------------------------------------------------------------------------------------------------------

func TestConvertFromPolyBigAndConvertToPolyBig(t *testing.T) {
	rand.Seed(42)
	const prec = uint(64)
	const q = uint64(97) // small prime ≡ 1 mod 2n for n=16
	const n = 16

	ringQ := makeSmallRing(n, q)

	for trial := 0; trial < 5; trial++ {
		// pick a small random polynomial
		P := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}

		// Convert to high-precision field element in Eval domain
		fEval := ConvertFromPolyBig(ringQ, P, prec)
		if fEval.Domain != Eval {
			t.Fatalf("ConvertFromPolyBig: Domain = %v, want Eval", fEval.Domain)
		}

		// Convert back to *ring.Poly
		P2 := ConvertToPolyBig(fEval, ringQ)

		// Compare P and P2 coefficient-by-coefficient
		for i := 0; i < n; i++ {
			if P.Coeffs[0][i] != P2.Coeffs[0][i] {
				t.Fatalf("Round-trip: P[%d]=%d, P2[%d]=%d", i, P.Coeffs[0][i], i, P2.Coeffs[0][i])
			}
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 4) NegacyclicEvaluatePoly / NegacyclicInterpolateElem round-trip
// ---------------------------------------------------------------------------------------------------------------------

func TestNegacyclicEvalInterpolate(t *testing.T) {
	rand.Seed(135)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	for trial := 0; trial < 5; trial++ {
		// pick a random polynomial in R = Z_q[x]/(x^n + 1)
		P := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}

		// Evaluate via negacyclic FFT
		fEval := NegacyclicEvaluatePoly(P, ringQ, prec)
		if fEval.Domain != Eval {
			t.Fatalf("NegacyclicEvaluatePoly: expected Domain=Eval, got %v", fEval.Domain)
		}

		// Interpolate back
		P2 := NegacyclicInterpolateElem(fEval, ringQ)

		// Compare P and P2
		for i := 0; i < n; i++ {
			if P.Coeffs[0][i] != P2.Coeffs[0][i] {
				t.Fatalf("Negacyclic round-trip mismatch at index %d: got %d, want %d", i, P2.Coeffs[0][i], P.Coeffs[0][i])
			}
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 5) ToEvalNegacyclic / ToCoeffNegacyclic round-trip
// ---------------------------------------------------------------------------------------------------------------------

func TestToEvalNegacyclicAndBack(t *testing.T) {
	rand.Seed(246)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	for trial := 0; trial < 5; trial++ {
		// build a random CyclotomicFieldElem in Coeff domain: lift a random ring.Poly
		P := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}
		// create a field-element whose Coeffs[i] = BigComplex(P[i], 0)
		e := NewFieldElemBig(n, prec)
		for i := 0; i < n; i++ {
			reBF := new(big.Float).SetPrec(prec).SetFloat64(float64(P.Coeffs[0][i]))
			e.Coeffs[i] = NewBigComplexFromFloat(reBF, new(big.Float).SetPrec(prec).SetFloat64(0))
		}
		e.Domain = Coeff

		// Convert to Eval
		fEval := ToEvalNegacyclic(e, ringQ, prec)
		if fEval.Domain != Eval {
			t.Fatalf("ToEvalNegacyclic: Domain = %v, want Eval", fEval.Domain)
		}

		// Convert back to Coeff
		e2 := ToCoeffNegacyclic(fEval, ringQ, prec)
		if e2.Domain != Coeff {
			t.Fatalf("ToCoeffNegacyclic: Domain = %v, want Coeff", e2.Domain)
		}

		// Compare e.Coeffs[i].Real to e2.Coeffs[i].Real
		for i := 0; i < n; i++ {
			orig, _ := e.Coeffs[i].Real.Float64()
			back, _ := e2.Coeffs[i].Real.Float64()
			if !approxEqual(orig, back, 1e-9) {
				t.Fatalf("ToEval/ToCoeff round-trip: index %d, got %v, want %v", i, back, orig)
			}
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 6) FieldAddBig, FieldSubBig, FieldMulBig, FieldScalarMulBig, FieldScalarDiv round-trip
// ---------------------------------------------------------------------------------------------------------------------

func TestFieldArithmeticSmall(t *testing.T) {
	const prec = uint(64)
	const n = 4

	// build two simple field elements in Coeff domain:
	e1 := NewFieldElemBig(n, prec)
	e2 := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		e1.Coeffs[i] = NewBigComplex(float64(i+1), 0, prec)      // [1,2,3,4]
		e2.Coeffs[i] = NewBigComplex(float64((i+1)*10), 0, prec) // [10,20,30,40]
	}
	e1.Domain = Coeff
	e2.Domain = Coeff

	// Add: e1 + e2 = [11,22,33,44]
	sum := FieldAddBig(e1, e2)
	for i := 0; i < n; i++ {
		v, _ := sum.Coeffs[i].Real.Float64()
		if !approxEqual(v, float64((i+1)+(i+1)*10), 1e-15) {
			t.Fatalf("FieldAddBig: index %d: got %v, want %v", i, v, float64((i+1)+(i+1)*10))
		}
	}

	// Sub: e2 - e1 = [9,18,27,36]
	diff := FieldSubBig(e2, e1)
	for i := 0; i < n; i++ {
		v, _ := diff.Coeffs[i].Real.Float64()
		want := float64((i+1)*10 - (i + 1))
		if !approxEqual(v, want, 1e-15) {
			t.Fatalf("FieldSubBig: index %d: got %v, want %v", i, v, want)
		}
	}

	// Mul: point-wise: [1*10, 2*20, 3*30, 4*40] = [10,40,90,160]
	mul := FieldMulBig(e1, e2)
	for i := 0; i < n; i++ {
		v, _ := mul.Coeffs[i].Real.Float64()
		want := float64(i+1) * float64((i+1)*10)
		if !approxEqual(v, want, 1e-15) {
			t.Fatalf("FieldMulBig: index %d: got %v, want %v", i, v, want)
		}
	}

	// Scalar multiply: multiply e1 by scalar s = (2+3i)
	s := NewBigComplex(2, 3, prec)
	scmul := FieldScalarMulBig(e1, s)
	for i := 0; i < n; i++ {
		// (i+1)*(2+3i) = ((i+1)*2) + ((i+1)*3)i:
		r, _ := scmul.Coeffs[i].Real.Float64()
		im, _ := scmul.Coeffs[i].Imag.Float64()
		wantR := float64(i+1) * 2.0
		wantI := float64(i+1) * 3.0
		if !approxEqual(r, wantR, 1e-15) || !approxEqual(im, wantI, 1e-15) {
			t.Fatalf("FieldScalarMulBig: index %d: got %v + %vi, want %v + %vi",
				i, r, im, wantR, wantI)
		}
	}

	// For FieldScalarDiv we need a list of norms.  Let’s pick a simple element d = [(2+0i),(3+0i),(4+0i),(5+0i)] in Eval domain:
	d := NewFieldElemBig(n, prec)
	d.Domain = Eval
	for i := 0; i < n; i++ {
		d.Coeffs[i] = NewBigComplex(float64(i+2), 0, prec) // [2,3,4,5]
	}

	// If your goal is simply to get 1 + 0i in every slot, use:
	inv := NewFieldElemBig(n, prec)
	norms := make([]*big.Float, n)
	for i := 0; i < n; i++ {
		// numerator = d[i], denominator = d[i]; d[i]/d[i] = 1
		inv.Coeffs[i] = d.Coeffs[i].Copy()
		norms[i] = new(big.Float).SetPrec(prec).SetFloat64(float64(i + 2)) // same as d[i].Real
	}
	dived := FieldScalarDiv(inv, norms)
	// Now dived.Coeffs[i] = d[i] / d[i] = 1 + 0i
	for i := 0; i < n; i++ {
		r, _ := dived.Coeffs[i].Real.Float64()
		im, _ := dived.Coeffs[i].Imag.Float64()
		if !approxEqual(r, 1.0, 1e-15) || !approxEqual(im, 0.0, 1e-15) {
			t.Fatalf("FieldScalarDiv: index %d: got %v + %vi, want 1 + 0i", i, r, im)
		}
	}

	// HermitianTranspose (Coefficient domain):
	e := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		e.Coeffs[i] = NewBigComplex(float64(i+1), float64(i+2), prec)
	}
	e.Domain = Coeff
	ht := HermitianTransposeFieldElem(e)

	for i := 0; i < n; i++ {
		var wantR, wantI float64
		if i == 0 {
			// out[0] = e.Coeffs[0] unchanged
			wantR, _ = e.Coeffs[0].Real.Float64()
			wantI, _ = e.Coeffs[0].Imag.Float64()
		} else {
			// out[i] = − e.Coeffs[n−i]
			orig := e.Coeffs[n-i]
			or, _ := orig.Real.Float64()
			oi, _ := orig.Imag.Float64()
			wantR = -or
			wantI = -oi
		}

		rHT, _ := ht.Coeffs[i].Real.Float64()
		iHT, _ := ht.Coeffs[i].Imag.Float64()
		if !approxEqual(rHT, wantR, 1e-15) || !approxEqual(iHT, wantI, 1e-15) {
			t.Fatalf("HermitianTranspose (Coeff): index %d: got %v + %vi, want %v + %vi",
				i, rHT, iHT, wantR, wantI)
		}
	}

	// Test PstrideBig: split even/odd.  For e = [1,2,3,4], f0=[1,3], f1=[2,4]
	e2 = NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		e2.Coeffs[i] = NewBigComplex(float64(i+1), 0, prec)
	}
	fe, fo := PstrideBig(e2)
	if fe.N != n/2 || fo.N != n/2 {
		t.Fatalf("PstrideBig: wrong lengths: got %d,%d want %d,%d", fe.N, fo.N, n/2, n/2)
	}
	for i := 0; i < n/2; i++ {
		r0, _ := fe.Coeffs[i].Real.Float64()
		r1, _ := fo.Coeffs[i].Real.Float64()
		want0 := float64(2*i + 1) // indices 0,2,4,6...
		want1 := float64(2*i + 2) // indices 1,3,5,7...
		if !approxEqual(r0, want0, 1e-15) || !approxEqual(r1, want1, 1e-15) {
			t.Fatalf("PstrideBig: slot %d: got f0=%v, f1=%v; want %v, %v",
				i, r0, r1, want0, want1)
		}
	}

	// SubScalar / AddScalar: subtract/add a constant
	g := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		g.Coeffs[i] = NewBigComplex(float64(5), float64(1), prec)
	}
	g2 := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		g2.Coeffs[i] = NewBigComplex(float64(2), float64(3), prec)
	}
	g.SubScalar(g2.Coeffs[0]) // now each slot = (5-2) + (1-3)i = 3 -2i
	for i := 0; i < n; i++ {
		rv, _ := g.Coeffs[i].Real.Float64()
		iv, _ := g.Coeffs[i].Imag.Float64()
		if !approxEqual(rv, 3.0, 1e-15) || !approxEqual(iv, -2.0, 1e-15) {
			t.Fatalf("SubScalar: slot %d: got %v+%vi, want 3 - 2i", i, rv, iv)
		}
	}
	g.AddScalar(g2.Coeffs[0]) // back to (5,1)
	for i := 0; i < n; i++ {
		rv, _ := g.Coeffs[i].Real.Float64()
		iv, _ := g.Coeffs[i].Imag.Float64()
		if !approxEqual(rv, 5.0, 1e-15) || !approxEqual(iv, 1.0, 1e-15) {
			t.Fatalf("AddScalar: slot %d: got %v+%vi, want 5 + 1i", i, rv, iv)
		}
	}

	// Copy and Conj methods
	h := e2.Copy()
	for i := 0; i < n; i++ {
		orig, _ := e2.Coeffs[i].Real.Float64()
		hv, _ := h.Coeffs[i].Real.Float64()
		if !approxEqual(orig, hv, 1e-15) {
			t.Fatalf("Copy: element %d changed: got %v, want %v", i, hv, orig)
		}
	}
	h2 := e2.Conj()
	for i := 0; i < n; i++ {
		rv, _ := h2.Coeffs[i].Real.Float64()
		iv, _ := h2.Coeffs[i].Imag.Float64()
		wantR := float64(i + 1)
		wantI := 0.0 // original imag was zero
		if !approxEqual(rv, wantR, 1e-15) || !approxEqual(iv, wantI, 1e-15) {
			t.Fatalf("Conj (elemwise) at %d: got %v+%vi, want %v+%vi", i, rv, iv, wantR, wantI)
		}
	}

	// SetCoeffs: copy coefficients from e1 into a fresh element
	e3 := NewFieldElemBig(n, prec)
	e3.SetCoeffs(e1)
	for i := 0; i < n; i++ {
		rv, _ := e3.Coeffs[i].Real.Float64()
		want := float64(i + 1)
		if i < len(e1.Coeffs) {
			want = float64(i + 1)
		}
		if !approxEqual(rv, want, 1e-15) {
			t.Fatalf("SetCoeffs: index %d: got %v, want %v", i, rv, want)
		}
	}

	// ExtractEven / ExtractOdd / InversePermute round-trip
	f := NewFieldElemBig(4, prec)
	// f = [10, 20, 30, 40]
	for i := 0; i < 4; i++ {
		f.Coeffs[i] = NewBigComplex(float64((i+1)*10), 0, prec)
	}
	fe2, fo2 := PstrideBig(f)
	// fe2 = [10,30], fo2 = [20,40]
	interleaved := NewFieldElemBig(4, prec)
	for i := 0; i < 2; i++ {
		interleaved.Coeffs[i] = fe2.Coeffs[i]
		interleaved.Coeffs[i+2] = fo2.Coeffs[i]
	}
	InversePermuteFieldElem(interleaved)
	// Now interleaved should be [10,20,30,40]
	for i := 0; i < 4; i++ {
		rv, _ := interleaved.Coeffs[i].Real.Float64()
		want := float64((i + 1) * 10)
		if !approxEqual(rv, want, 1e-15) {
			t.Fatalf("InversePermuteFieldElem: index %d: got %v, want %v", i, rv, want)
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 7) ConvertFrom/ConvertTo using FFT vs. NegacyclicEvaluate/NegacyclicInterpolate consistency
// ---------------------------------------------------------------------------------------------------------------------

func TestConsistency_FFT_vs_Negacyclic(t *testing.T) {
	rand.Seed(999)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	// Pick a random polynomial P
	P := ringQ.NewPoly()
	for i := 0; i < n; i++ {
		P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
	}

	// 1) Evaluate P via ConvertFromPolyBig / ConvertToPolyBig (regular FFT)
	fEvalFFT := ConvertFromPolyBig(ringQ, P, prec)
	P2 := ConvertToPolyBig(fEvalFFT, ringQ)

	// 2) Evaluate P via NegacyclicEvaluatePoly / NegacyclicInterpolateElem
	fEvalN := NegacyclicEvaluatePoly(P, ringQ, prec)
	P3 := NegacyclicInterpolateElem(fEvalN, ringQ)

	// P2 and P3 both represent “P mod (x^n+1)” under two different transforms.
	// Since P is arbitrary, and our ring is truly (x^n+1), we expect P2 ≡ P3 ≡ P as coefficient vectors.
	for i := 0; i < n; i++ {
		if P.Coeffs[0][i] != P2.Coeffs[0][i] {
			t.Fatalf("FFT round-trip mismatch at %d: P2[%d]=%d, P[%d]=%d", i, i, P2.Coeffs[0][i], i, P.Coeffs[0][i])
		}
		if P.Coeffs[0][i] != P3.Coeffs[0][i] {
			t.Fatalf("Negacyclic round-trip mismatch at %d: P3[%d]=%d, P[%d]=%d", i, i, P3.Coeffs[0][i], i, P.Coeffs[0][i])
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 8) SwitchToEval / SwitchToCoeff consistency (using plain FFTBig pipeline)
// ---------------------------------------------------------------------------------------------------------------------

func TestSwitchToEvalAndBack(t *testing.T) {
	rand.Seed(1234)
	const prec = uint(256)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	// Build a random CyclotomicFieldElem in Coeff domain:
	coeffElem := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		// set e.Coeffs[i] = BigComplex(random integer mod q, 0)
		x := float64(rand.Intn(int(q)))
		coeffElem.Coeffs[i] = NewBigComplex(x, 0, prec)
	}
	coeffElem.Domain = Coeff

	// Switch to Eval (using FFTBig pipeline, via ConvertFromPolyBig)
	// First convert coeffElem→ring.Poly
	P := ringQ.NewPoly()
	for i := 0; i < n; i++ {
		val, _ := coeffElem.Coeffs[i].Real.Float64()
		ri := int64(math.Round(val)) % int64(q)
		if ri < 0 {
			ri += int64(q)
		}
		P.Coeffs[0][i] = uint64(ri)
	}
	// now fEval := ConvertFromPolyBig
	fEval := ConvertFromPolyBig(ringQ, P, prec)
	if fEval.Domain != Eval {
		t.Fatalf("ConvertFromPolyBig: expected Domain=Eval, got %v", fEval.Domain)
	}

	// Overwrite coeffElem with fEval.Coeffs and set its Domain=Eval
	coeffElem.Coeffs = fEval.Coeffs
	coeffElem.Domain = Eval

	// Now switch back to Coeff using the regular (non-negacyclic) inverse FFT:
	coeffPoly := ConvertToPolyBig(coeffElem, ringQ)

	// Check that coeffPoly.Coeffs[i].Real ≈ original P[i]
	for i := 0; i < n; i++ {
		val := coeffPoly.Coeffs[0][i]
		want := float64(P.Coeffs[0][i])
		if !approxEqual(float64(val), want, 1e-9) {
			t.Fatalf("SwitchToCoeff: index %d: got %v, want %v", i, val, want)
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 9) SwitchToEval (regular FFT) vs. ToEvalNegacyclic: cross-comparison
// ---------------------------------------------------------------------------------------------------------------------

func TestSwitchToEvalVsToEvalNegacyclic(t *testing.T) {
	rand.Seed(777)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	// pick a random CyclotomicFieldElem in Coeff domain
	elemC := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		x := float64(rand.Intn(int(q)))
		elemC.Coeffs[i] = NewBigComplex(x, 0, prec)
	}
	elemC.Domain = Coeff

	// Convert to Eval via full-length FFT (ConvertFromPolyBig)
	// first build a *ring.Poly
	P := ringQ.NewPoly()
	for i := 0; i < n; i++ {
		v, _ := elemC.Coeffs[i].Real.Float64()
		ri := int64(math.Round(v)) % int64(q)
		if ri < 0 {
			ri += int64(q)
		}
		P.Coeffs[0][i] = uint64(ri)
	}
	fEvalFFT := ConvertFromPolyBig(ringQ, P, prec)

	// Convert to Eval via Negacyclic (ToEvalNegacyclic)
	fEvalNeg := ToEvalNegacyclic(elemC.Copy(), ringQ, prec)

	// Now fEvalFFT[k] corresponds to p(ω^k), where ω = e^{-2πi/(2n)}.
	// fEvalNeg[k] corresponds to p(ω^{2k+1}).  They are different transforms,
	// so we only check that *both* are in Domain=Eval and have length n.
	if fEvalFFT.Domain != Eval {
		t.Fatalf("ConvertFromPolyBig: Domain=%v, want Eval", fEvalFFT.Domain)
	}
	if fEvalNeg.Domain != Eval {
		t.Fatalf("ToEvalNegacyclic: Domain=%v, want Eval", fEvalNeg.Domain)
	}
	if fEvalFFT.N != n || fEvalNeg.N != n {
		t.Fatalf("Eval lengths: got %d and %d, want %d", fEvalFFT.N, fEvalNeg.N, n)
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 10) InversePermuteFieldElem alone
// ---------------------------------------------------------------------------------------------------------------------

func TestInversePermuteFieldElemAlone(t *testing.T) {
	const prec = uint(64)
	const n = 6 // must be even
	e := NewFieldElemBig(n, prec)
	// put [0,1,2,3,4,5], then after InversePermute → [0,3,1,4,2,5]
	for i := 0; i < n; i++ {
		e.Coeffs[i] = NewBigComplex(float64(i), 0, prec)
	}
	InversePermuteFieldElem(e)
	expected := []float64{0, 3, 1, 4, 2, 5}
	for i := 0; i < n; i++ {
		v, _ := e.Coeffs[i].Real.Float64()
		if !approxEqual(v, expected[i], 1e-15) {
			t.Fatalf("InversePermuteFieldElem: index %d: got %v, want %v", i, v, expected[i])
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 11) HermitianTranspose in Eval domain vs. direct Coeff transpose
// ---------------------------------------------------------------------------------------------------------------------

func TestHermitianTransposeEvalConsistency(t *testing.T) {
	rand.Seed(2024)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	for trial := 0; trial < 5; trial++ {
		// 1) Pick a random real-coeff polynomial P in R = Z_q[x]/(x^n+1)
		P := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}

		// 2) Compute its direct Coeff-domain transpose Q_coeff:
		//    Q_coeff[0] = P[0], Q_coeff[i] = q - P[n-i]  (mod q), for i=1..n-1
		Q_coeff := ringQ.NewPoly()
		Q_coeff.Coeffs[0][0] = P.Coeffs[0][0] % q
		for i := 1; i < n; i++ {
			qi := P.Coeffs[0][n-i] % q
			if qi != 0 {
				Q_coeff.Coeffs[0][i] = (q - qi) % q
			} else {
				Q_coeff.Coeffs[0][i] = 0
			}
		}

		// 3) Convert P to a FieldElem in Eval:
		coeffElem := NewFieldElemBig(n, prec)
		for i := 0; i < n; i++ {
			val := float64(P.Coeffs[0][i])
			coeffElem.Coeffs[i] = NewBigComplex(val, 0, prec)
		}
		coeffElem.Domain = Coeff
		fEval := ToEvalNegacyclic(coeffElem, ringQ, prec) // Domain=Eval

		// 4) Apply HermitianTransposeFieldElem in Eval:
		htEval := HermitianTransposeFieldElem(fEval)
		if htEval.Domain != Eval {
			t.Fatalf("HermitianTransposeFieldElem did not return Eval-domain; got %v", htEval.Domain)
		}

		// 5) Convert htEval back to Coeff:
		htCoeff := ToCoeffNegacyclic(htEval, ringQ, prec)

		// 6) Compare htCoeff to Q_coeff
		for i := 0; i < n; i++ {
			got, _ := htCoeff.Coeffs[i].Real.Float64()
			want := float64(Q_coeff.Coeffs[0][i])
			if got != want {
				t.Fatalf("Trial %d: HermitianTranspose Eval mismatch at %d: got %f, want %f", trial, i, got, want)
			}
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 12) Test FieldInverseDiagWithNorm & FieldScalarDiv (invD check) in isolation
// ---------------------------------------------------------------------------------------------------------------------

func TestInverseDiagAndScalarDiv(t *testing.T) {
	rand.Seed(12345)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	for trial := 0; trial < 5; trial++ {
		// 1) Pick a random polynomial P in R = Z_q[x]/(x^n+1)
		P := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}

		// 2) Convert P to a FieldElem in Eval:
		coeffElem := NewFieldElemBig(n, prec)
		for i := 0; i < n; i++ {
			val := float64(P.Coeffs[0][i])
			coeffElem.Coeffs[i] = NewBigComplex(val, 0, prec)
		}
		coeffElem.Domain = Coeff
		fEval := ToEvalNegacyclic(coeffElem, ringQ, prec) // Domain=Eval

		// 3) Compute invD, norms := FieldInverseDiagWithNorm(fEval),
		//    then invEval := FieldScalarDiv(invD, norms)
		invD, norms := FieldInverseDiagWithNorm(fEval)
		invEval := FieldScalarDiv(invD, norms)
		invEval.Domain = Eval

		// 4) Check inv: fEval * invEval = 1 (pointwise)
		prod := FieldMulBig(fEval, invEval)
		for i := 0; i < n; i++ {
			r, _ := prod.Coeffs[i].Real.Float64()
			im, _ := prod.Coeffs[i].Imag.Float64()
			if !approxEqual(r, 1.0, 1e-9) || !approxEqual(im, 0.0, 1e-9) {
				t.Fatalf("Trial %d: inversion failed at slot %d: got %v+%vi, want 1+0i", trial, i, r, im)
			}
		}

		// 5) Check norms: norms[i] = fEval.Coeffs[i].Real^2 + fEval.Coeffs[i].Imag^2
		for i := 0; i < n; i++ {
			// recompute abs^2 of fEval
			rf, _ := fEval.Coeffs[i].Real.Float64()
			ifm, _ := fEval.Coeffs[i].Imag.Float64()
			wantNorm := rf*rf + ifm*ifm

			gotNorm, _ := norms[i].Float64()
			if !approxEqual(gotNorm, wantNorm, 1e-9) {
				t.Fatalf("Trial %d: norms mismatch at slot %d: got %v, want %v", trial, i, gotNorm, wantNorm)
			}
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 13) Test RingFreeToCoeffNegacyclic / ToCoeffNegacyclic inversion correctness
// ---------------------------------------------------------------------------------------------------------------------

func TestRingFreeToCoeffNegacyclicInversion(t *testing.T) {
	rand.Seed(555)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16

	ringQ := makeSmallRing(n, q)

	for trial := 0; trial < 5; trial++ {
		// 1) Pick a random polynomial P in R = Z_q[x]/(x^n+1)
		P := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}

		// 2) Convert P to a FieldElem in Eval via ToEvalNegacyclic
		coeffElem := NewFieldElemBig(n, prec)
		for i := 0; i < n; i++ {
			val := float64(P.Coeffs[0][i])
			coeffElem.Coeffs[i] = NewBigComplex(val, 0, prec)
		}
		coeffElem.Domain = Coeff
		fEval := ToEvalNegacyclic(coeffElem, ringQ, prec)
		// for i := 0; i < n; i++ {
		// 	fmt.Printf("fEval[%d] = %v\n", i, fEval.Coeffs[i])
		// }

		// 3) Now convert fEval back via ToCoeffNegacyclic
		back := ToCoeffNegacyclic(fEval, ringQ, prec)

		// for i := 0; i < n; i++ {
		// 	fmt.Printf("back.Coeffs[%d] = %v\n", i, back.Coeffs[i])
		// }
		// Compare back.Coeffs[i].Real to original P.Coeffs[0][i]
		for i := 0; i < n; i++ {
			got, _ := back.Coeffs[i].Real.Float64()
			want := float64(P.Coeffs[0][i])
			if !approxEqual(got, want, 1e-9) {
				t.Fatalf("Trial %d: ToCoeffNegacyclic mismatch at %d: got %v, want %v", trial, i, got, want)
			}
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 14) ComplexFieldLoop: full Eval arithmetic vs. naive coefficient arithmetic, extended
// ---------------------------------------------------------------------------------------------------------------------

func TestComplexFieldLoop(t *testing.T) {
	rand.Seed(2025)
	const prec = uint(64)
	const q = uint64(97)
	const n = 16
	ringQ := makeSmallRing(n, q)

	trials := 10
	for trial := 0; trial < trials; trial++ {
		// 1) Pick a random polynomial P in R = ℤ₉₇[x]/(xⁿ+1)
		P := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			P.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}

		// 2) Build a CyclotomicFieldElem “coeffElem” in Coeff domain from P
		coeffElem := NewFieldElemBig(n, prec)
		coeffElem.Domain = Coeff
		for i := 0; i < n; i++ {
			coeffElem.Coeffs[i] = NewBigComplex(float64(P.Coeffs[0][i]), 0, prec)
		}

		// 3) Convert coeffElem → Eval domain using negacyclic transform
		fEval := ToEvalNegacyclic(coeffElem, ringQ, prec)
		if fEval.Domain != Eval {
			t.Fatalf("trial %d: ToEvalNegacyclic: expected Domain=Eval, got %v", trial, fEval.Domain)
		}

		// 4) Inversion in Eval domain: compute invNumerator = conj(fEval[i]), norms = |fEval[i]|²
		invNumerator, norms := FieldInverseDiagWithNorm(fEval)
		fInv := FieldScalarDiv(invNumerator, norms)

		// 5) Check fEval * fInv = 1 in every slot
		prod := FieldMulBig(fEval, fInv)
		for i := 0; i < n; i++ {
			r, _ := prod.Coeffs[i].Real.Float64()
			im, _ := prod.Coeffs[i].Imag.Float64()
			if !approxEqual(r, 1.0, 1e-9) || !approxEqual(im, 0.0, 1e-9) {
				t.Fatalf("trial %d, inversion: slot %d: got %v+%vi, want 1+0i", trial, i, r, im)
			}
		}

		// 6) Check addition: 2·P mod (xⁿ+1) vs (fEval + fEval) → ToCoeffNegacyclic
		//    Naive coefficient‐domain doubling:
		Q := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			Q.Coeffs[0][i] = (P.Coeffs[0][i] * 2) % q
		}
		//    Eval‐domain doubling:
		f2Eval := FieldAddBig(fEval, fEval)
		f2Eval.Domain = Eval
		polyF2 := ToCoeffNegacyclic(f2Eval, ringQ, prec)
		for i := 0; i < n; i++ {
			got := polyF2.Coeffs[i].Real
			gotFloat, _ := got.Float64()
			want := Q.Coeffs[0][i]
			if gotFloat != float64(want) {
				t.Fatalf("trial %d, add round‐trip: index %d: got %v, want %d", trial, i, gotFloat, want)
			}
		}

		// 7) Check multiplication: P·P mod (xⁿ+1) vs (fEval * fEval) → ToCoeffNegacyclic
		//    Naive coefficient‐domain convolution modulo xⁿ+1:
		R := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				prodCoeff := (P.Coeffs[0][i] * P.Coeffs[0][j]) % q
				idx := (i + j) % n
				if i+j >= n {
					// xⁿ ≡ −1
					prodCoeff = (q - prodCoeff) % q
				}
				R.Coeffs[0][idx] = (R.Coeffs[0][idx] + prodCoeff) % q
			}
		}
		//    Eval‐domain squaring:
		fSquaredEval := FieldMulBig(fEval, fEval)
		fSquaredEval.Domain = Eval
		polySquared := ToCoeffNegacyclic(fSquaredEval, ringQ, prec)
		for i := 0; i < n; i++ {
			got, _ := polySquared.Coeffs[i].Real.Float64()
			want := float64(R.Coeffs[0][i])
			if got != want {
				t.Fatalf("trial %d, mul round‐trip: index %d: got %v, want %v", trial, i, got, want)
			}
		}
	}

}

func TestNTTInvNTTNoScale(t *testing.T) {
	const n = 16
	const q = uint64(97)

	// Build a ring of degree 16, modulus 97 (97 ≡ 1 mod 2·16)
	ringQ := makeSmallRing(n, q)

	// Create a polynomial P with deterministic coeffs [1,2,3,…,16]
	P := ringQ.NewPoly()
	for i := 0; i < n; i++ {
		P.Coeffs[0][i] = uint64(i + 1)
	}

	// Copy P to Pcopy for later comparison
	Pcopy := ringQ.NewPoly()
	for i := 0; i < n; i++ {
		Pcopy.Coeffs[0][i] = P.Coeffs[0][i]
	}

	// Apply NTT then InvNTT in place
	ringQ.NTT(P, P)
	ringQ.InvNTT(P, P)

	// After NTT→InvNTT, P should match Pcopy exactly (no scaling factor)
	for i := 0; i < n; i++ {
		if P.Coeffs[0][i] != Pcopy.Coeffs[0][i] {
			t.Fatalf("NTT/InvNTT mismatch at index %d: got %d, want %d",
				i, P.Coeffs[0][i], Pcopy.Coeffs[0][i])
		}
	}
}

// helper: generate a random coefficient vector with prec bits
func randomCoeffElem(m int, prec uint) *CyclotomicFieldElem {
	out := NewFieldElemBig(m, prec)
	for i := 0; i < m; i++ {
		val := (rand.Float64()*2 - 1) * 1000 // ±1000 range
		out.Coeffs[i] = NewBigComplex(val, 0, prec)
	}
	out.Domain = Coeff
	return out
}

func TestFloatFFTRoundTrip(t *testing.T) {
	rand.Seed(time.Now().UnixNano())

	const (
		n    = 8   // ring dimension
		prec = 128 // test precision in bits
	)

	orig := randomCoeffElem(n, prec)

	// coeff → eval → coeff
	eval := FloatToEvalNegacyclic(orig, prec)
	back := FloatToCoeffNegacyclic(eval, prec)
	for i := 0; i < n; i++ {
		fmt.Printf("orig[%d] = %s, back[%d] = %s\n", i, back.Coeffs[i].Real.Text('g', 10), i, orig.Coeffs[i].Real.Text('g', 10))
	}
	// tolerance = 2^(−prec+10)
	eps := new(big.Float).SetPrec(prec).SetFloat64(1)
	eps.SetMantExp(eps, int(-prec+10))

	for i := 0; i < n; i++ {
		diff := new(big.Float).Sub(orig.Coeffs[i].Real, back.Coeffs[i].Real)
		if diff.Abs(diff).Cmp(eps) > 0 {
			t.Fatalf("mismatch at coeff %d: want %s, got %s",
				i,
				orig.Coeffs[i].Real.Text('g', 10),
				back.Coeffs[i].Real.Text('g', 10))
		}
		if back.Coeffs[i].Imag.Sign() != 0 {
			t.Fatalf("imaginary residue at coeff %d: %s",
				i, back.Coeffs[i].Imag.Text('g', 10))
		}
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 15) Eval-domain transpose should reverse the evaluation slots
func TestTransposeEvalDomain(t *testing.T) {
	const (
		n    = 16
		prec = 128
	)
	rand.Seed(4242)

	for trial := 0; trial < 5; trial++ {
		r := randomCoeffElem(n, prec)           // Coeff
		rEval := FloatToEvalNegacyclic(r, prec) // Eval
		rT := HermitianTransposeFieldElem(rEval)

		if rT.Domain != Eval {
			t.Fatalf("trial %d: transpose output not in Eval domain", trial)
		}
		for i := 0; i < n; i++ {
			exp := rEval.Coeffs[n-i-1]
			got := rT.Coeffs[i]
			rExp, _ := exp.Real.Float64()
			iExp, _ := exp.Imag.Float64()
			rGot, _ := got.Real.Float64()
			iGot, _ := got.Imag.Float64()
			if !approxEqual(rGot, rExp, 1e-12) || !approxEqual(iGot, iExp, 1e-12) {
				t.Fatalf("trial %d: slot %d: got %v+%vi, want %v+%vi", trial, i, rGot, iGot, rExp, iExp)
			}
		}
	}

}

// ---------------------------------------------------------------------------------------------------------------------
// 16) Scalar (n = 1) sanity-check for Sample2zField
// ---------------------------------------------------------------------------------------------------------------------

func TestSample2zFieldScalarCase(t *testing.T) {
	const (
		prec = 128
		q    = uint64(65537)
	)
	n := 1

	// a = d = s² = 25, b = 0  (all Eval)
	a := NewFieldElemBig(n, prec)
	d := NewFieldElemBig(n, prec)
	b := NewFieldElemBig(n, prec)
	s2 := NewBigComplex(25, 0, prec)
	a.Coeffs[0], d.Coeffs[0] = s2.Copy(), s2.Copy()
	a.Domain, b.Domain, d.Domain = Eval, Eval, Eval

	// centres c0 = c1 = 0 (Coeff)
	c0 := NewFieldElemBig(n, prec)
	c1 := NewFieldElemBig(n, prec)
	c0.Domain, c1.Domain = Coeff, Coeff

	q0, q1 := Sample2zField(a, b, d, c0, c1, n, q, prec)
	if q0 == nil || q1 == nil {
		t.Fatalf("Sample2zField returned nil outputs")
	}
	// variance used to draw q0 is aPr = 25 > 0; check printed in function already
	val, _ := q1.Coeffs[0].Real.Float64()
	if math.Abs(val) > 1000 { // 6σ guardrail for std≈5
		t.Fatalf("Sample2zField scalar: q1 sample too large: %v", val)
	}
}

// ---------------------------------------------------------------------------------------------------------------------
// 17) Even/Odd split and recomposition in coefficient domain
// ---------------------------------------------------------------------------------------------------------------------

func TestEvenOddSplitRoundTrip(t *testing.T) {
	const (
		n    = 16
		prec = 64
	)
	rand.Seed(606060)

	for trial := 0; trial < 5; trial++ {
		p := randomCoeffElem(n, prec) // Coeff
		fe := p.Copy().ExtractEven()  // still Coeff
		fo := p.Copy().ExtractOdd()

		// Recompose:  g(x) = fe(x²) + x·fo(x²)
		recomp := NewFieldElemBig(n, prec)
		for i := 0; i < n/2; i++ {
			recomp.Coeffs[2*i] = fe.Coeffs[i].Copy()
			recomp.Coeffs[2*i+1] = fo.Coeffs[i].Copy()
		}
		recomp.Domain = Coeff

		for i := 0; i < n; i++ {
			orig, _ := p.Coeffs[i].Real.Float64()
			back, _ := recomp.Coeffs[i].Real.Float64()
			if !approxEqual(orig, back, 1e-9) {
				t.Fatalf("trial %d: coeff mismatch at %d: got %v, want %v", trial, i, back, orig)
			}
		}
	}
}

// InverseFieldElem computes the coefficient-wise inverse of f.
func InverseFieldElem(f *CyclotomicFieldElem) *CyclotomicFieldElem {
	n := f.N
	prec := f.Coeffs[0].Real.Prec()
	out := NewFieldElemBig(n, prec)
	out.Domain = f.Domain
	for i := 0; i < n; i++ {
		out.Coeffs[i] = f.Coeffs[i].Inv()
	}
	return out
}

// shiftRightCoeff performs the coefficient-domain right shift used in PALISADE.
func shiftRightCoeff(e *CyclotomicFieldElem) *CyclotomicFieldElem {
	if e.Domain != Coeff {
		panic("shiftRightCoeff: input must be in Coeff domain")
	}
	n := e.N
	prec := e.Coeffs[0].Real.Prec()
	out := NewFieldElemBig(n, prec)
	for i := 0; i < n-1; i++ {
		out.Coeffs[i+1] = e.Coeffs[i].Copy()
	}
	// out[0] = -e[n-1]
	negR := new(big.Float).SetPrec(prec).Neg(e.Coeffs[n-1].Real)
	negI := new(big.Float).SetPrec(prec).Neg(e.Coeffs[n-1].Imag)
	out.Coeffs[0] = &BigComplex{Real: negR, Imag: negI}
	out.Domain = Coeff
	return out
}

// -----------------------------------------------------------------------------
// Basic Field2n-like operations
// -----------------------------------------------------------------------------

func TestPALGetFormat(t *testing.T) {
	const prec = uint(64)
	f := NewFieldElemBig(2, prec)
	if f.Domain != Coeff {
		t.Fatalf("GetFormat mismatch: got %v", f.Domain)
	}
}

func TestPALInverse(t *testing.T) {
	const prec = uint(64)
	f := NewFieldElemBig(2, prec)
	f.Domain = Eval
	f.Coeffs[0] = NewBigComplex(2, 1, prec)
	f.Coeffs[1] = NewBigComplex(-4, -2, prec)

	inv := InverseFieldElem(f)
	want0 := NewBigComplex(0.4, -0.2, prec)
	want1 := NewBigComplex(-0.2, 0.1, prec)

	r0, _ := inv.Coeffs[0].Real.Float64()
	i0, _ := inv.Coeffs[0].Imag.Float64()
	wr0, _ := want0.Real.Float64()
	wi0, _ := want0.Imag.Float64()
	if !approxEqual(r0, wr0, 1e-12) || !approxEqual(i0, wi0, 1e-12) {
		t.Fatalf("inverse slot0: got %v+%vi, want %v+%vi", r0, i0, wr0, wi0)
	}
	r1, _ := inv.Coeffs[1].Real.Float64()
	i1, _ := inv.Coeffs[1].Imag.Float64()
	wr1, _ := want1.Real.Float64()
	wi1, _ := want1.Imag.Float64()
	if !approxEqual(r1, wr1, 1e-12) || !approxEqual(i1, wi1, 1e-12) {
		t.Fatalf("inverse slot1: got %v+%vi, want %v+%vi", r1, i1, wr1, wi1)
	}
}

func TestPALPlus(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(2, prec)
	b := NewFieldElemBig(2, prec)
	a.Domain, b.Domain = Eval, Eval
	a.Coeffs[0] = NewBigComplex(2, 1, prec)
	a.Coeffs[1] = NewBigComplex(-4, 2, prec)
	b.Coeffs[0] = NewBigComplex(3, -0.1, prec)
	b.Coeffs[1] = NewBigComplex(-4, 3.2, prec)

	sum := FieldAddBig(a, b)
	want := []*BigComplex{
		NewBigComplex(5, 0.9, prec),
		NewBigComplex(-8, 5.2, prec),
	}
	for i := 0; i < 2; i++ {
		r, _ := sum.Coeffs[i].Real.Float64()
		im, _ := sum.Coeffs[i].Imag.Float64()
		wr, _ := want[i].Real.Float64()
		wi, _ := want[i].Imag.Float64()
		if !approxEqual(r, wr, 1e-12) || !approxEqual(im, wi, 1e-12) {
			t.Fatalf("plus[%d]: got %v+%vi, want %v+%vi", i, r, im, wr, wi)
		}
	}
}

func TestPALScalarPlus(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(2, prec)
	a.Domain = Coeff
	a.Coeffs[0] = NewBigComplex(2, 0, prec)
	a.Coeffs[1] = NewBigComplex(-4, 0, prec)

	scalar := NewBigComplex(3.2, 0, prec)
	scalarSum := a.Copy()
	scalarSum.Coeffs[0] = scalarSum.Coeffs[0].Add(scalar)

	want := []*BigComplex{NewBigComplex(5.2, 0, prec), NewBigComplex(-4, 0, prec)}
	for i := 0; i < 2; i++ {
		r, _ := scalarSum.Coeffs[i].Real.Float64()
		im, _ := scalarSum.Coeffs[i].Imag.Float64()
		wr, _ := want[i].Real.Float64()
		wi, _ := want[i].Imag.Float64()
		if !approxEqual(r, wr, 1e-12) || !approxEqual(im, wi, 1e-12) {
			t.Fatalf("scalar_plus[%d]: got %v+%vi, want %v+%vi", i, r, im, wr, wi)
		}
	}
}

func TestPALMinus(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(2, prec)
	b := NewFieldElemBig(2, prec)
	a.Domain, b.Domain = Eval, Eval
	a.Coeffs[0] = NewBigComplex(2, 1, prec)
	a.Coeffs[1] = NewBigComplex(-4, 2, prec)
	b.Coeffs[0] = NewBigComplex(3, -0.1, prec)
	b.Coeffs[1] = NewBigComplex(-4, 3.2, prec)

	diff := FieldSubBig(a, b)
	want := []*BigComplex{
		NewBigComplex(-1, 1.1, prec),
		NewBigComplex(0, -1.2, prec),
	}
	for i := 0; i < 2; i++ {
		r, _ := diff.Coeffs[i].Real.Float64()
		im, _ := diff.Coeffs[i].Imag.Float64()
		wr, _ := want[i].Real.Float64()
		wi, _ := want[i].Imag.Float64()
		if !approxEqual(r, wr, 1e-9) || !approxEqual(im, wi, 1e-9) {
			t.Fatalf("minus[%d]: got %v+%vi, want %v+%vi", i, r, im, wr, wi)
		}
	}
}

func TestPALTimes(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(2, prec)
	b := NewFieldElemBig(2, prec)
	a.Domain, b.Domain = Eval, Eval
	a.Coeffs[0] = NewBigComplex(4, 3, prec)
	a.Coeffs[1] = NewBigComplex(6, -3, prec)
	b.Coeffs[0] = NewBigComplex(4, -3, prec)
	b.Coeffs[1] = NewBigComplex(4, -2.8, prec)

	prod := FieldMulBig(a, b)
	want := []*BigComplex{
		NewBigComplex(25, 0, prec),
		NewBigComplex(15.6, -28.8, prec),
	}
	for i := 0; i < 2; i++ {
		r, _ := prod.Coeffs[i].Real.Float64()
		im, _ := prod.Coeffs[i].Imag.Float64()
		wr, _ := want[i].Real.Float64()
		wi, _ := want[i].Imag.Float64()
		if !approxEqual(r, wr, 1e-9) || !approxEqual(im, wi, 1e-9) {
			t.Fatalf("times[%d]: got %v+%vi, want %v+%vi", i, r, im, wr, wi)
		}
	}
}

func TestPALTimesWithSwitch(t *testing.T) {
	const prec = uint(64)
	const n = 4
	rand.Seed(42)
	a := NewFieldElemBig(n, prec)
	b := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		a.Coeffs[i] = NewBigComplex(1, 0, prec)
	}
	b.Coeffs[0] = NewBigComplex(1, 0, prec)
	b.Coeffs[1] = NewBigComplex(0, 0, prec)
	b.Coeffs[2] = NewBigComplex(1, 0, prec)
	b.Coeffs[3] = NewBigComplex(0, 0, prec)
	a.Domain = Coeff
	b.Domain = Coeff

	aEval := FloatToEvalNegacyclic(a, prec)
	bEval := FloatToEvalNegacyclic(b, prec)
	prod := FieldMulBig(aEval, bEval)
	prod.Domain = Eval
	prodCoeff := FloatToCoeffNegacyclic(prod, prec)

	want := []*BigComplex{
		NewBigComplex(0, 0, prec),
		NewBigComplex(0, 0, prec),
		NewBigComplex(2, 0, prec),
		NewBigComplex(2, 0, prec),
	}
	for i := 0; i < n; i++ {
		r, _ := prodCoeff.Coeffs[i].Real.Float64()
		if !approxEqual(r, mustFloat64(want[i].Real), 1e-9) {
			t.Fatalf("times_with_switch[%d]: got %v, want %v", i, r, mustFloat64(want[i].Real))
		}
	}
}

func mustFloat64(x *big.Float) float64 {
	v, _ := x.Float64()
	return v
}

func TestPALShiftRight(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(4, prec)
	a.Domain = Coeff
	a.Coeffs[0] = NewBigComplex(4, 0, prec)
	a.Coeffs[1] = NewBigComplex(3, 0, prec)
	a.Coeffs[2] = NewBigComplex(2, 0, prec)
	a.Coeffs[3] = NewBigComplex(1, 0, prec)

	out := shiftRightCoeff(a)
	want := []*BigComplex{
		NewBigComplex(-1, 0, prec),
		NewBigComplex(4, 0, prec),
		NewBigComplex(3, 0, prec),
		NewBigComplex(2, 0, prec),
	}
	for i := 0; i < 4; i++ {
		r, _ := out.Coeffs[i].Real.Float64()
		if !approxEqual(r, mustFloat64(want[i].Real), 1e-12) {
			t.Fatalf("shift_right[%d]: got %v, want %v", i, r, mustFloat64(want[i].Real))
		}
	}
}

func TestPALTransposeCoeff(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(4, prec)
	a.Domain = Coeff
	a.Coeffs[0] = NewBigComplex(4, 0, prec)
	a.Coeffs[1] = NewBigComplex(3, 0, prec)
	a.Coeffs[2] = NewBigComplex(2, 0, prec)
	a.Coeffs[3] = NewBigComplex(1, 0, prec)

	t2 := HermitianTransposeFieldElem(a)
	want := []*BigComplex{
		NewBigComplex(4, 0, prec),
		NewBigComplex(-1, 0, prec),
		NewBigComplex(-2, 0, prec),
		NewBigComplex(-3, 0, prec),
	}
	for i := 0; i < 4; i++ {
		r, _ := t2.Coeffs[i].Real.Float64()
		if !approxEqual(r, mustFloat64(want[i].Real), 1e-12) {
			t.Fatalf("transpose_coeff[%d]: got %v, want %v", i, r, mustFloat64(want[i].Real))
		}
	}
}

func TestPALTransposeEval(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(4, prec)
	a.Domain = Coeff
	a.Coeffs[0] = NewBigComplex(4, 0, prec)
	a.Coeffs[1] = NewBigComplex(3, 0, prec)
	a.Coeffs[2] = NewBigComplex(2, 0, prec)
	a.Coeffs[3] = NewBigComplex(1, 0, prec)
	aEval := FloatToEvalNegacyclic(a, prec)
	tEval := HermitianTransposeFieldElem(aEval)
	back := FloatToCoeffNegacyclic(tEval, prec)

	want := []*BigComplex{
		NewBigComplex(4, 0, prec),
		NewBigComplex(-1, 0, prec),
		NewBigComplex(-2, 0, prec),
		NewBigComplex(-3, 0, prec),
	}
	for i := 0; i < 4; i++ {
		r, _ := back.Coeffs[i].Real.Float64()
		if !approxEqual(r, mustFloat64(want[i].Real), 1e-9) {
			t.Fatalf("transpose_eval[%d]: got %v, want %v", i, r, mustFloat64(want[i].Real))
		}
	}
}

func TestPALExtractOddEven(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(4, prec)
	a.Domain = Coeff
	a.Coeffs[0] = NewBigComplex(4, 0, prec)
	a.Coeffs[1] = NewBigComplex(3, 0, prec)
	a.Coeffs[2] = NewBigComplex(2, 0, prec)
	a.Coeffs[3] = NewBigComplex(1, 0, prec)
	odd := a.ExtractOdd()
	even := a.ExtractEven()

	wantOdd := []*BigComplex{NewBigComplex(3, 0, prec), NewBigComplex(1, 0, prec)}
	wantEven := []*BigComplex{NewBigComplex(4, 0, prec), NewBigComplex(2, 0, prec)}
	for i := 0; i < 2; i++ {
		r, _ := odd.Coeffs[i].Real.Float64()
		if !approxEqual(r, mustFloat64(wantOdd[i].Real), 1e-12) {
			t.Fatalf("extract_odd[%d]: got %v, want %v", i, r, mustFloat64(wantOdd[i].Real))
		}
		r, _ = even.Coeffs[i].Real.Float64()
		if !approxEqual(r, mustFloat64(wantEven[i].Real), 1e-12) {
			t.Fatalf("extract_even[%d]: got %v, want %v", i, r, mustFloat64(wantEven[i].Real))
		}
	}
}

func TestPALInversePermute(t *testing.T) {
	const prec = uint(64)
	e := NewFieldElemBig(6, prec)
	e.Domain = Coeff
	for i := 0; i < 6; i++ {
		e.Coeffs[i] = NewBigComplex(float64(i), 0, prec)
	}
	InversePermuteFieldElem(e)
	expect := []float64{0, 3, 1, 4, 2, 5}
	for i := 0; i < 6; i++ {
		r, _ := e.Coeffs[i].Real.Float64()
		if !approxEqual(r, expect[i], 1e-12) {
			t.Fatalf("inverse_permute[%d]: got %v, want %v", i, r, expect[i])
		}
	}
}

func TestPALScalarMult(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(4, prec)
	a.Domain = Eval
	a.Coeffs[0] = NewBigComplex(1, -1, prec)
	a.Coeffs[1] = NewBigComplex(3, -2, prec)
	a.Coeffs[2] = NewBigComplex(2, -3, prec)
	a.Coeffs[3] = NewBigComplex(4, -4, prec)

	scalar := NewBigComplex(3, 0, prec)
	out := FieldScalarMulBig(a, scalar)
	want := []*BigComplex{
		NewBigComplex(3, -3, prec),
		NewBigComplex(9, -6, prec),
		NewBigComplex(6, -9, prec),
		NewBigComplex(12, -12, prec),
	}
	for i := 0; i < 4; i++ {
		r, _ := out.Coeffs[i].Real.Float64()
		im, _ := out.Coeffs[i].Imag.Float64()
		wr, _ := want[i].Real.Float64()
		wi, _ := want[i].Imag.Float64()
		if !approxEqual(r, wr, 1e-12) || !approxEqual(im, wi, 1e-12) {
			t.Fatalf("scalar_mult[%d]: got %v+%vi, want %v+%vi", i, r, im, wr, wi)
		}
	}
}

func TestPALCoeffToEval(t *testing.T) {
	const prec = uint(64)
	a := NewFieldElemBig(8, prec)
	a.Domain = Coeff
	vals := []float64{4, 5, 5, 4.2, 5, 7.1, 6, 3}
	for i, v := range vals {
		a.Coeffs[i] = NewBigComplex(v, 0, prec)
	}
	aEval := FloatToEvalNegacyclic(a, prec)
	if aEval.Domain != Eval {
		t.Fatalf("expected Eval domain")
	}
	back := FloatToCoeffNegacyclic(aEval, prec)
	for i := 0; i < 8; i++ {
		orig, _ := a.Coeffs[i].Real.Float64()
		got, _ := back.Coeffs[i].Real.Float64()
		if !approxEqual(orig, got, 1e-9) {
			t.Fatalf("coeff_eval[%d]: got %v, want %v", i, got, orig)
		}
	}
}

func TestPALEvalToCoeff(t *testing.T) {
	const prec = uint(64)
	b := NewFieldElemBig(8, prec)
	b.Domain = Eval
	valsR := []float64{4.03087, 8.15172, 1.26249, 2.55492, 2.55492, 1.26249, 8.15172, 4.03087}
	valsI := []float64{26.2795, 5.84489, 0.288539, 0.723132, -0.723132, -0.288539, -5.84489, -26.2795}
	for i := 0; i < 8; i++ {
		b.Coeffs[i] = NewBigComplex(valsR[i], valsI[i], prec)
	}
	bCoeff := FloatToCoeffNegacyclic(b, prec)
	if bCoeff.Domain != Coeff {
		t.Fatalf("expected Coeff domain")
	}
	back := FloatToEvalNegacyclic(bCoeff, prec)
	for i := 0; i < 8; i++ {
		r, _ := back.Coeffs[i].Real.Float64()
		im, _ := back.Coeffs[i].Imag.Float64()
		wr := valsR[i]
		wi := valsI[i]
		if !approxEqual(r, wr, 1e-6) || !approxEqual(im, wi, 1e-6) {
			t.Fatalf("eval_coeff[%d]: got %v+%vi, want %v+%vi", i, r, im, wr, wi)
		}
	}
}
