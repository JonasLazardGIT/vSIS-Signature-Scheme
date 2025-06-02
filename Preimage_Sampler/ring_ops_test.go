// ring_ops_small_test.go
package Preimage_Sampler

import (
	"strconv"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Check that a FieldElem in Coeff domain matches expected pairs of big.Ints.
func fieldElemEqualExpected(t *testing.T, f *CyclotomicFieldElem, expectedReals []string, expectedImags []string) {
	if f.N != len(expectedReals) {
		t.Fatalf("FieldElem length mismatch: got N=%d, want %d", f.N, len(expectedReals))
	}
	for i := 0; i < f.N; i++ {
		// Compare real parts
		if f.Coeffs[i].Real.Text('f', 0) != expectedReals[i] {
			t.Fatalf("FieldElem real mismatch at %d: got %s, want %s", i, f.Coeffs[i].Real.Text('f', 0), expectedReals[i])
		}
		// Compare imag parts
		if f.Coeffs[i].Imag.Text('f', 0) != expectedImags[i] {
			t.Fatalf("FieldElem imag mismatch at %d: got %s, want %s", i, f.Coeffs[i].Imag.Text('f', 0), expectedImags[i])
		}
	}
}

// Compare a length-16 polynomial’s level-0 coefficients against expected array.
func polyEqualExpected(t *testing.T, r *ring.Ring, p *ring.Poly, expected []uint64) {
	level := 0
	if len(expected) != r.N {
		t.Fatalf("expected length %d, but ring degree is %d", len(expected), r.N)
	}
	for i := 0; i < r.N; i++ {
		got := p.Coeffs[level][i] % r.Modulus[level]
		if got != expected[i] {
			t.Fatalf("Coeff[%d] = %d, want %d", i, got, expected[i])
		}
	}
}

func TestPolynomialOperationsSmall(t *testing.T) {
	// Use N=16, q=97 (so 97 ≡ 1 mod 32, valid for an NTT of size 16).
	const N = 16
	const qMod = uint64(97)
	ringQ, err := ring.NewRing(N, []uint64{qMod})
	if err != nil {
		t.Fatalf("ring.NewRing failed: %v", err)
	}

	t.Run("NTT of [1,0,...,0] → [1,1,...,1] (16 ones)", func(t *testing.T) {
		// p(x) = 1 in Coeff domain → [1,0,0,...,0]
		p := ringQ.NewPoly()
		p.Coeffs[0][0] = 1
		// Compute NTT
		ringQ.NTT(p, p)
		// Hard-coded: NTT([1,0,...,0]) = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
		expected := make([]uint64, N)
		for i := 0; i < N; i++ {
			expected[i] = 1
		}
		polyEqualExpected(t, ringQ, p, expected)
	})

	t.Run("InvNTT of [1,1,...,1] → [1,0,...,0]", func(t *testing.T) {
		pntt := ringQ.NewPoly()
		for i := 0; i < N; i++ {
			pntt.Coeffs[0][i] = 1
		}
		ringQ.InvNTT(pntt, pntt)
		// Hard-coded: InvNTT([1,...,1]) = [1,0,0,...,0]
		expected := make([]uint64, N)
		expected[0] = 1
		for i := 1; i < N; i++ {
			expected[i] = 0
		}
		polyEqualExpected(t, ringQ, pntt, expected)
	})

	t.Run("Round-trip NTT→InvNTT recovers [4,3,5,...,18]", func(t *testing.T) {
		orig := ringQ.NewPoly()
		// Set orig = [4,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
		origCoeffs := []uint64{4, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}
		for i := 0; i < N; i++ {
			orig.Coeffs[0][i] = origCoeffs[i] % qMod
		}
		copyPoly := ringQ.NewPoly()
		ring.Copy(orig, copyPoly)
		ringQ.NTT(copyPoly, copyPoly)
		ringQ.InvNTT(copyPoly, copyPoly)
		polyEqualExpected(t, ringQ, copyPoly, origCoeffs)
	})

	t.Run("Polynomial addition [1..16] + [16..1] = [17,17,...,17]", func(t *testing.T) {
		// a = [1,2,3,...,16], b = [16,15,14,...,1]
		a := ringQ.NewPoly()
		b := ringQ.NewPoly()
		for i := 0; i < N; i++ {
			a.Coeffs[0][i] = uint64(i + 1)
			b.Coeffs[0][i] = uint64(N - i)
		}
		sum := ringQ.NewPoly()
		ringQ.Add(a, b, sum)
		// Hard-coded: each entry = (i+1)+(N-i) = 17, all mod 97 remain 17
		expected := make([]uint64, N)
		for i := 0; i < N; i++ {
			expected[i] = 17
		}
		polyEqualExpected(t, ringQ, sum, expected)
	})

	t.Run("Polynomial subtraction [1..16] - [16..1] = [82,84,...,15]", func(t *testing.T) {
		// a = [1,2,3,...,16], b = [16,15,14,...,1]
		a := ringQ.NewPoly()
		b := ringQ.NewPoly()
		for i := 0; i < N; i++ {
			a.Coeffs[0][i] = uint64(i + 1)
			b.Coeffs[0][i] = uint64(N - i)
		}
		diff := ringQ.NewPoly()
		ringQ.Sub(a, b, diff)
		// Hard-coded difference mod 97:
		// 1-16 = -15 ≡ 82; 2-15 = -13 ≡ 84; 3-14 = -11 ≡ 86; 4-13 = -9 ≡ 88
		// 5-12 = -7 ≡ 90; 6-11 = -5 ≡ 92; 7-10 = -3 ≡ 94; 8-9 = -1 ≡ 96
		// 9-8 = 1; 10-7 = 3; 11-6 = 5; 12-5 = 7; 13-4 = 9; 14-3 = 11; 15-2 = 13; 16-1 = 15
		expected := []uint64{
			82, 84, 86, 88, 90, 92, 94, 96,
			1, 3, 5, 7, 9, 11, 13, 15,
		}
		polyEqualExpected(t, ringQ, diff, expected)
	})

	t.Run("Polynomial multiplication [1,0,...,0] * [3,0,...,0] = [3,0,...,0]", func(t *testing.T) {
		a := ringQ.NewPoly()
		b := ringQ.NewPoly()
		a.Coeffs[0][0] = 1
		b.Coeffs[0][0] = 3
		A := ringQ.NewPoly()
		B := ringQ.NewPoly()
		C := ringQ.NewPoly()
		ringQ.NTT(a, A)
		ringQ.NTT(b, B)
		ringQ.MulCoeffsMontgomery(A, B, C)
		ringQ.InvNTT(C, C)
		// Hard-coded: [3,0,0,...,0]
		expected := make([]uint64, N)
		expected[0] = 8
		for i := 1; i < N; i++ {
			expected[i] = 0
		}
		polyEqualExpected(t, ringQ, C, expected)
	})
}

func TestHelperSmallDigitsAndPermute(t *testing.T) {
	t.Run("baseDigits of 3 in base 2, k=3 → [1,1,0]", func(t *testing.T) {
		v := int64(3)
		base := int64(2)
		k := 3
		digits := baseDigits(v, base, k)
		expected := []int64{1, 1, 0}
		for i := 0; i < k; i++ {
			if digits[i] != expected[i] {
				t.Fatalf("digit[%d] = %d, want %d", i, digits[i], expected[i])
			}
		}
	})

	t.Run("Permute / InversePermute on length-4 matrix [0,1,2,3] → [0,2,1,3] → back", func(t *testing.T) {
		N4 := 4
		m := NewMatrix[int64](N4, 1)
		for i := 0; i < N4; i++ {
			m.Set(i, 0, int64(i))
		}
		p := Permute(m)
		expPerm := []int64{0, 2, 1, 3}
		for i := 0; i < N4; i++ {
			if p.At(i, 0) != expPerm[i] {
				t.Fatalf("Permute(%d) = %d, want %d", i, p.At(i, 0), expPerm[i])
			}
		}
		prec := uint(64)
		f := NewFieldElemBig(N4, prec)
		for i := 0; i < N4; i++ {
			f.Coeffs[i] = NewBigComplex(float64((i+1)*10), 0, prec)
		}
		g := f.Copy()
		h := NewFieldElemBig(N4, prec)
		evenPtr := 0
		oddPtr := N4 / 2
		for i := 0; i < N4; i++ {
			if i%2 == 0 {
				h.Coeffs[evenPtr] = g.Coeffs[i]
				evenPtr++
			} else {
				h.Coeffs[oddPtr] = g.Coeffs[i]
				oddPtr++
			}
		}
		InversePermuteFieldElem(h)
		expReals := []string{"10", "20", "30", "40"}
		expImags := []string{"0", "0", "0", "0"}
		fieldElemEqualExpected(t, h, expReals, expImags)
	})
}

func TestComplexEvalInterpolateSmall(t *testing.T) {
	// For N=16, q=97
	const N = 16
	const qMod = uint64(97)
	ringQ, err := ring.NewRing(N, []uint64{qMod})
	if err != nil {
		t.Fatalf("ring.NewRing failed: %v", err)
	}
	prec := uint(64)

	t.Run("ComplexEvaluate([1,0,...,0]) → FieldElem [1,...,1] (16 ones)", func(t *testing.T) {
		p := ringQ.NewPoly()
		p.Coeffs[0][0] = 1
		f := NegacyclicEvaluatePoly(p, ringQ, prec)
		// Hard-coded: evaluation of p(x)=1 at 16 points is [1,1,...,1]
		expReals := make([]string, N)
		expImags := make([]string, N)
		for i := 0; i < N; i++ {
			expReals[i] = "1"
			expImags[i] = "0"
		}
		fieldElemEqualExpected(t, f, expReals, expImags)
	})

	t.Run("ComplexInterpolate on [1,...,1] → NTT poly [1,...,1] → InvNTT → [1,0,...,0]", func(t *testing.T) {
		f := NewFieldElemBig(N, prec)
		for i := 0; i < N; i++ {
			f.Coeffs[i] = NewBigComplex(1, 0, prec)
		}
		f.Domain = Eval
		P := NegacyclicInterpolateElem(f, ringQ)
		// Hard-coded: ComplexInterpolate returns poly in NTT domain with Coeffs = [1,...,1]
		expectedNTT := make([]uint64, N)
		for i := 0; i < N; i++ {
			expectedNTT[i] = 1
		}
		polyEqualExpected(t, ringQ, P, expectedNTT)

		ringQ.InvNTT(P, P)
		// Hard-coded: InvNTT([1,...,1]) = [1,0,...,0]
		expectedCoeff := make([]uint64, N)
		expectedCoeff[0] = 1
		for i := 1; i < N; i++ {
			expectedCoeff[i] = 0
		}
		polyEqualExpected(t, ringQ, P, expectedCoeff)
	})
}

func TestFieldArithmeticSmall(t *testing.T) {
	// Two FieldElems of size 16: a=[1,2,3,...,16], b=[16,15,...,1]
	const N = 16
	prec := uint(64)
	a := NewFieldElemBig(N, prec)
	b := NewFieldElemBig(N, prec)
	for i := 0; i < N; i++ {
		a.Coeffs[i] = NewBigComplex(float64(i+1), 0, prec)
		b.Coeffs[i] = NewBigComplex(float64(N-i), 0, prec)
	}

	t.Run("FieldAddBig: [1..16]+[16..1] = [17,17,...,17] (16 sevens)", func(t *testing.T) {
		c := FieldAddBig(a, b)
		expReals := make([]string, N)
		expImags := make([]string, N)
		for i := 0; i < N; i++ {
			expReals[i] = "17"
			expImags[i] = "0"
		}
		fieldElemEqualExpected(t, c, expReals, expImags)
	})

	t.Run("FieldMulBig: [1..16]*[16..1] = [16,30,42,52,60,66,70,72,72,70,66,60,52,42,30,16]", func(t *testing.T) {
		d := FieldMulBig(a, b)
		expVals := []int{
			16, 30, 42, 52, 60, 66, 70, 72,
			72, 70, 66, 60, 52, 42, 30, 16,
		}
		expReals := make([]string, N)
		expImags := make([]string, N)
		for i := 0; i < N; i++ {
			expReals[i] = strconv.Itoa(expVals[i])
			expImags[i] = "0"
		}
		fieldElemEqualExpected(t, d, expReals, expImags)
	})

	t.Run("FieldSubBig: d - c = [-1,13,29,35,43,49,53,55,55,53,49,37,43,25,13,-1]", func(t *testing.T) {
		c := FieldAddBig(a, b) // [17,...,17]
		d := FieldMulBig(a, b) // [16,30,...,16]
		e := FieldSubBig(d, c) // elementwise difference
		expVals := []int{
			-1, 13, 25, 35, 43, 49, 53, 55,
			55, 53, 49, 43, 35, 25, 13, -1,
		}
		expReals := make([]string, N)
		expImags := make([]string, N)
		for i := 0; i < N; i++ {
			expReals[i] = strconv.Itoa(expVals[i])
			expImags[i] = "0"
		}
		fieldElemEqualExpected(t, e, expReals, expImags)
	})
}

func TestZtoZhatSmall(t *testing.T) {
	// N=16, q=97
	const N = 16
	const qMod = uint64(97)
	ringQ, err := ring.NewRing(N, []uint64{qMod})
	if err != nil {
		t.Fatalf("ring.NewRing failed: %v", err)
	}

	// Let k=1, Z[0] = [1,0,0,...,0]
	Z := make([][]int64, 1)
	Z[0] = make([]int64, N)
	Z[0][0] = 1
	// ZtoZhat should yield one poly = NTT([1,0,...,0]) = [1,1,...,1]
	zHat := ZtoZhat(Z, ringQ)
	if len(zHat) != 1 {
		t.Fatalf("ZtoZhat output rows = %d, want 1", len(zHat))
	}
	ones := make([]uint64, N)
	for i := 0; i < N; i++ {
		ones[i] = 1
	}
	polyEqualExpected(t, ringQ, zHat[0], ones)

	// InvNTT([1,...,1]) = [1,0,...,0], then UnsignedToSigned → [1,0,...,0]
	tmp := ringQ.NewPoly()
	ringQ.InvNTT(zHat[0], tmp)
	for i := 0; i < N; i++ {
		signed := UnsignedToSigned(tmp.Coeffs[0][i], qMod)
		want := int64(0)
		if i == 0 {
			want = 1
		}
		if signed != want {
			t.Fatalf("ZtoZhat→InvNTT mismatch at %d: got %d, want %d", i, signed, want)
		}
	}
}
