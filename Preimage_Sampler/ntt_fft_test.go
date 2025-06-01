// fft_ntt_test.go
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

// maxAbsDiff computes the maximum absolute difference between two real‐valued slices.
func maxAbsDiff(a, b []complex128) float64 {
	if len(a) != len(b) {
		panic("length mismatch")
	}
	max := 0.0
	for i := range a {
		d := math.Abs(real(a[i]) - real(b[i]))
		if d > max {
			max = d
		}
		d2 := math.Abs(imag(a[i]) - imag(b[i]))
		if d2 > max {
			max = d2
		}
	}
	return max
}

// TestFFTvsNTT tests, for several random integer polynomials, the exactness of
// (1) NTT ↔ InvNTT (should be exact in Zq) and (2) FFT ↔ IFFT (floating–point).
func TestFFTvsNTT(t *testing.T) {
	// Parameters for the test:
	const (
		N    = 512               // ring dimension (power of two)
		qMod = uint64(1<<30 + 1) // a prime ≡ 1 mod 2N; here 1073741825 = 2^30 +1 is not prime,
		// so pick a known NTT‐friendly prime, e.g.  2147483649 ≡ 1 mod 1024?
		// For safety, use one that Lattigo accepts. Let's use  2013265921 = 15·2^27+1,
		// which supports N up to 2^27. We'll just do N=512.
	)
	var (
		primeList = []uint64{2013265921} // 2013265921 ≡ 1 mod 2^27, so certainly ≡1 mod 1024
	)
	// Build ringQ
	ringQ, err := ring.NewRing(N, primeList)
	if err != nil {
		t.Fatalf("ring.NewRing(%d, %v) failed: %v", N, primeList, err)
	}

	// Number of random polynomials to test
	const trials = 10

	// Helper: for each trial, build a random integer polynomial p (coeffs in [0, q−1]),
	// then:
	//   1) check NTT↔InvNTT: exactly recover p
	//   2) check FFT↔IFFT:
	//       a) represent p in ℤ (by SignedToSigned), build a length‐N complex slice
	//          vals[i] = p_i (as float64), call FFT(vals,N) → eval, then IFFT(eval,N) → rec,
	//          measure max|rec[i]−p_i|. That should be O(1e−12 · N) or so in double precision.
	//       b) also compare rec (rounded) to exact integer p_i, count how many entries deviate by ≥0.5.

	for trial := 0; trial < trials; trial++ {
		// 1) generate a random "coeff domain" polynomial p
		p := ringQ.NewPoly()
		for i := 0; i < N; i++ {
			p.Coeffs[0][i] = rand.Uint64() % qMod
		}

		// 2) NTT ↔ InvNTT round‐trip check
		cpy := ringQ.NewPoly()
		// Copy p → cpy
		for i := 0; i < N; i++ {
			cpy.Coeffs[0][i] = p.Coeffs[0][i]
		}
		// Forward NTT, then InvNTT
		ringQ.NTT(cpy, cpy)
		ringQ.InvNTT(cpy, cpy)
		// Now cpy should match p exactly
		for i := 0; i < N; i++ {
			if cpy.Coeffs[0][i] != p.Coeffs[0][i] {
				t.Fatalf("NTT↔InvNTT failed at trial %d, index %d: got %d, want %d",
					trial, i, cpy.Coeffs[0][i], p.Coeffs[0][i])
			}
		}

		// 3) FFT ↔ IFFT check
		// Build a slice of length N of complex128: vals[i] = signed(p_i)
		vals := make([]complex128, N)
		for i := 0; i < N; i++ {
			// Convert to "signed" in [−q/2 … +q/2]
			s := UnsignedToSigned(p.Coeffs[0][i], qMod)
			vals[i] = complex(float64(s), 0)
		}
		// Compute FFT(vals) → eval
		eval := FFT(vals, N)
		// Compute IFFT(eval) → rec
		rec := IFFT(eval, N)

		// 3a) measure raw floating‐point error: max|real(rec[i]) − vals[i].Real|
		maxErr := maxAbsDiff(vals, rec)
		t.Logf("trial %d: FFT↔IFFT max absolute error = %g", trial, maxErr)

		// 3b) count how many coefficients would round "incorrectly" to the nearest integer
		wrongCount := 0
		for i := 0; i < N; i++ {
			orig := real(vals[i])
			rounded := math.Round(real(rec[i]))
			if math.Abs(rounded-orig) > 0.49 {
				wrongCount++
			}
		}
		if wrongCount > 0 {
			t.Logf("  → %d coefficients would round incorrectly (of %d) at trial %d",
				wrongCount, N, trial)
		} else {
			t.Logf("  → all %d coefficients would round correctly at trial %d", N, trial)
		}
	}
}

const (
	nTest     = 16
	qTest     = uint64(97) // ≡ 1 mod 32, so n=16 works
	maxTrials = 20
)

// modInverseUint64 returns a^{-1} mod m (if it exists), or 0 if not invertible.
func modInverseUint64(a, m uint64) uint64 {
	ai := new(big.Int).SetUint64(a)
	mi := new(big.Int).SetUint64(m)
	inv := new(big.Int).ModInverse(ai, mi)
	if inv == nil {
		return 0
	}
	return inv.Uint64()
}

// invertPolyNTT computes the exact inverse of A(X) mod (X^n+1, q) by NTT → pointwise → InvNTT.
// Returns the inverse in *NTT‐domain*.  If A is not invertible, returns an error.
func invertPolyNTT(Acoeff *ring.Poly, ringQ *ring.Ring) (*ring.Poly, error) {
	n := ringQ.N

	// 1) Compute A_ntt = NTT_q(Acoeff)
	AcoeffCopy := ringQ.NewPoly()
	copy(AcoeffCopy.Coeffs[0], Acoeff.Coeffs[0])
	ringQ.NTT(AcoeffCopy, AcoeffCopy)

	// 2) Form Ainv_ntt by inverting each slot in Z_q
	AinvNTT := ringQ.NewPoly()
	for i := 0; i < n; i++ {
		slot := AcoeffCopy.Coeffs[0][i]
		if slot == 0 {
			return nil, fmt.Errorf("zero encountered at NTT slot %d; polynomial not invertible", i)
		}
		inv := modInverseUint64(slot, ringQ.Modulus[0])
		if inv == 0 {
			return nil, fmt.Errorf("no inverse for %d mod %d", slot, ringQ.Modulus[0])
		}
		AinvNTT.Coeffs[0][i] = inv
	}

	// 3) Optionally verify that InvNTT(AinvNTT) is indeed the coefficient‐vector of A^{-1}:
	//    But here we only need the NTT‐domain form.  If needed, one can do:
	//    invCoe := ringQ.NewPoly(); ringQ.InvNTT(AinvNTT, invCoe); ringQ.NTT(invCoe, AinvNTT)

	return AinvNTT, nil
}

// fftInvert computes “FFT‐based” inversion: Apply ComplexEvaluate → invert complex slots → ComplexInterpolate → NTT.
func fftInvert(Acoeff *ring.Poly, ringQ *ring.Ring, prec uint) *ring.Poly {
	n := ringQ.N

	// (a) Compute A_coeff explicitly in coefficient domain (in case Acoeff is already NTT).
	AcoeffCopy := ringQ.NewPoly()
	copy(AcoeffCopy.Coeffs[0], Acoeff.Coeffs[0])
	// We want AcoeffCopy in coefficient form; assume input is coefficient form.

	// (b) ComplexEvaluate → returns a FieldElem in evaluation domain (complex double vector).
	fEval := ComplexEvaluate(AcoeffCopy, ringQ, prec)
	// Now fEval.Domain == Eval.

	// (c) Invert each complex slot: z → 1/z = (a - ib)/(a^2 + b^2).
	fInv := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		slot := fEval.Coeffs[i]
		a, _ := slot.Real.Float64() // Coeffs[i] is a complex number
		b, _ := slot.Imag.Float64() // Coeffs[i] is a complex number
		den := a*a + b*b
		// Avoid divide-by-zero; if den is very small, rounding will be huge.
		fInv.Coeffs[i] = NewBigComplex(a/den, -b/den, prec)
	}
	fInv.Domain = Eval

	// (d) ComplexInterpolate: returns a *NTT‐domain* polynomial under modulus qTest.
	//     That is, ComplexInterpolate(fInv) yields an NTT‐form ring.Poly so that InvNTT would give the coefficient form.
	PinvNTT := ComplexInterpolate(fInv, ringQ)

	return PinvNTT // already NTT‐domain under qTest
}

func TestFFTvsExactInverse(t *testing.T) {
	rand.Seed(time.Now().UnixNano())

	// Build ringQ mod qTest
	ringQ, err := ring.NewRing(nTest, []uint64{qTest})
	if err != nil {
		t.Fatalf("failed to create ringQ: %v", err)
	}

	prec := uint(64) // 64‐bit precision for BigFloat

	for trial := 0; trial < maxTrials; trial++ {
		// 1) Pick a random Acoeff in coefficient form mod 17.
		Acoeff := ringQ.NewPoly()
		for i := 0; i < nTest; i++ {
			Acoeff.Coeffs[0][i] = uint64(rand.Intn(int(qTest)))
		}

		// 2) Ensure A is invertible: compute exact inversion via NTT method.
		AinvExactNTT, err := invertPolyNTT(Acoeff, ringQ)
		if err != nil {
			// If non‐invertible, skip this trial.
			continue
		}

		// 3) Compute floating‐point FFT‐based inverse:
		AinvFFTNTT := fftInvert(Acoeff, ringQ, prec)

		// 4) Compare every slot of AinvExactNTT vs AinvFFTNTT.
		//    Because NTT‐domain values must match exactly in Z_q, any difference ≠0 indicates rounding error.
		for i := 0; i < nTest; i++ {
			exact := AinvExactNTT.Coeffs[0][i]
			fft := AinvFFTNTT.Coeffs[0][i]
			if exact != fft {
				t.Logf("trial %d, slot %d: exact inverse = %d, FFT‐inv = %d (difference = %d)",
					trial, i, exact, fft, (fft+qTest-exact)%qTest)
				t.Fail()
				break
			}
		}
	}
}
