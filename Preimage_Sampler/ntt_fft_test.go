// fft_ntt_test.go
package Preimage_Sampler

import (
	"math"
	"math/big"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
)

//--------------------------------------------------------------------------------------------------
// TestFFTvsNTT: NTT↔InvNTT and FFTBig↔IFFTBig (as before), plus additional “benchmark” operations.
//--------------------------------------------------------------------------------------------------

func TestFFTvsNTT(t *testing.T) {
	const (
		N    = 512
		qMod = uint64(2013265921) // prime ≡ 1 mod 1024
	)
	primeList := []uint64{qMod}

	// Build ringQ
	ringQ, err := ring.NewRing(N, primeList)
	if err != nil {
		t.Fatalf("ring.NewRing(%d, %v) failed: %v", N, primeList, err)
	}

	const trials = 10
	prec := uint(128) // 128-bit precision for BigFloat

	for trial := 0; trial < trials; trial++ {
		// 1) Generate a random polynomial p ∈ R_q[x]/(x^N+1), coefficients in [0, qMod)
		p := ringQ.NewPoly()
		for i := 0; i < N; i++ {
			p.Coeffs[0][i] = rand.Uint64() % qMod
		}

		// 2) NTT ↔ InvNTT round-trip (must be exact in ℤ/q):
		buf := ringQ.NewPoly()
		copy(buf.Coeffs[0], p.Coeffs[0])
		ringQ.NTT(buf, buf)
		ringQ.InvNTT(buf, buf)
		for i := 0; i < N; i++ {
			if buf.Coeffs[0][i] != p.Coeffs[0][i] {
				t.Fatalf("NTT↔InvNTT failed at trial %d, index %d: got %d, want %d",
					trial, i, buf.Coeffs[0][i], p.Coeffs[0][i])
			}
		}

		// 3) FFTBig ↔ IFFTBig check:
		//    (a) Build signed integer slice orig[i] = UnsignedToSigned(p_i, qMod).
		orig := make([]float64, N)
		for i := 0; i < N; i++ {
			orig[i] = float64(UnsignedToSigned(p.Coeffs[0][i], qMod))
		}

		//    (b) Lift orig to []*BigComplex and call FFTBig.
		bigIn := make([]*BigComplex, N)
		for i := 0; i < N; i++ {
			reBF := new(big.Float).SetPrec(prec).SetFloat64(orig[i])
			imBF := new(big.Float).SetPrec(prec).SetFloat64(0.0)
			bigIn[i] = NewBigComplexFromFloat(reBF, imBF)
		}
		evalBig := FFTBig(bigIn, prec)

		//    (c) Call IFFTBig(evalBig) → coeffBig.
		coeffBig := IFFTBig(evalBig, prec)

		//    (d) Convert each coeffBig[i].Real back to float64, round, compare to orig[i].
		recovered := make([]float64, N)
		maxErr := 0.0
		wrongCount := 0
		for i := 0; i < N; i++ {
			f64, _ := coeffBig[i].Real.Float64()
			recovered[i] = f64
			errVal := math.Abs(f64 - orig[i])
			if errVal > maxErr {
				maxErr = errVal
			}
			if math.Abs(math.Round(f64)-orig[i]) > 0.49 {
				wrongCount++
			}
		}
		t.Logf("trial %d: FFTBig↔IFFTBig max float error = %.5e; misrounds = %d/%d",
			trial, maxErr, wrongCount, N)
	}
}

//--------------------------------------------------------------------------------------------------
// TestFieldInverseDiagWithNorm: as before, verifies that aEval·invEval ≈ 1 slotwise.
//--------------------------------------------------------------------------------------------------

func TestFieldInverseDiagWithNorm(t *testing.T) {
	const n = 16
	const q = uint64(97) // 97 ≡ 1 mod 32, supports NTT size 16
	ringQ, err := ring.NewRing(n, []uint64{q})
	if err != nil {
		t.Fatalf("Failed to make ring: %v", err)
	}

	prec := uint(64)
	rand.Seed(time.Now().UnixNano())

	for trial := 0; trial < 20; trial++ {
		// (a) Pick a random nonzero polynomial aCoeff in coeff domain
		aCoeff := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			aCoeff.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}

		// (b) Embed into evaluation (FFT) domain:
		aEval := ComplexEvaluate(aCoeff, ringQ, prec)
		aEval.Domain = Eval

		// (c) Compute invConj,norms := FieldInverseDiagWithNorm(aEval)
		invConj, norms := FieldInverseDiagWithNorm(aEval)
		invEval := FieldScalarDiv(invConj, norms)
		invEval.Domain = Eval

		// (d) Multiply aEval * invEval slotwise → prod
		prod := FieldMulBig(aEval, invEval) // domain=Eval

		// Check that each slot of prod ≈ 1
		tol := new(big.Float).SetPrec(prec).SetFloat64(1e-30)
		for i := 0; i < n; i++ {
			re := prod.Coeffs[i].Real
			im := prod.Coeffs[i].Imag

			one := new(big.Float).SetPrec(prec).SetFloat64(1.0)
			zero := new(big.Float).SetPrec(prec).SetFloat64(0.0)

			diffRe := new(big.Float).Sub(re, one)
			diffIm := new(big.Float).Sub(im, zero)

			if diffRe.Abs(diffRe).Cmp(tol) > 0 || diffIm.Abs(diffIm).Cmp(tol) > 0 {
				r64, _ := re.Float64()
				i64, _ := im.Float64()
				t.Fatalf("trial %d, slot %d: aEval·invEval ≠ 1 (got % .6g + % .6gi)",
					trial, i, r64, i64)
			}
		}

		// (e) Finally, apply ComplexInterpolate(prod) → get onesPolyNTT, then InvNTT → onesCoeff
		onesPolyNTT := ComplexInterpolate(prod, ringQ)
		onesCoeff := ringQ.NewPoly()
		ringQ.InvNTT(onesPolyNTT, onesCoeff)
		for i := 0; i < n; i++ {
			if onesCoeff.Coeffs[0][i] != 1 {
				t.Fatalf("trial %d, slot %d: after Interpolate→InvNTT, expected 1 but got %d",
					trial, i, onesCoeff.Coeffs[0][i])
			}
		}
	}
}

//--------------------------------------------------------------------------------------------------
// TestFFTFieldOperationsAccuracy
//
// Now we add a suite of operations (multiply, invert, add) performed via FFT<Big> + field‐ops,
// and compare each result against the exact “NTT‐based” result.  Any discrepancy indicates
// bit-level inaccuracies in FFTBig/IFFTBig or in the chi-squared coeff→evaluation→coeff round-trip.
//--------------------------------------------------------------------------------------------------

func TestFFTFieldOperationsAccuracy(t *testing.T) {
	const (
		N    = 512
		qMod = uint64(2013265921)
	)
	primeList := []uint64{qMod}

	ringQ, err := ring.NewRing(N, primeList)
	if err != nil {
		t.Fatalf("ring.NewRing(%d, %v) failed: %v", N, primeList, err)
	}

	prec := uint(128)
	rand.Seed(time.Now().UnixNano())

	trials := 10
	for trial := 0; trial < trials; trial++ {
		//------------------------------------------------------------
		// 1) Pick two random polynomials p, q ∈ R_q[x]/(x^N+1)
		//------------------------------------------------------------
		p := ringQ.NewPoly()
		q := ringQ.NewPoly()
		for i := 0; i < N; i++ {
			p.Coeffs[0][i] = rand.Uint64() % qMod
			q.Coeffs[0][i] = rand.Uint64() % qMod
		}

		//------------------------------------------------------------
		// 2) EXACT polynomial multiplication via NTT:
		//    pNTT = NTT(p), qNTT = NTT(q)
		//    exactMulNTT = pNTT ∘ qNTT (pointwise)
		//    InvNTT(exactMulNTT) = exactMulCoefficients
		//------------------------------------------------------------
		pNTT := ringQ.NewPoly()
		qNTT := ringQ.NewPoly()
		copy(pNTT.Coeffs[0], p.Coeffs[0])
		copy(qNTT.Coeffs[0], q.Coeffs[0])
		ringQ.NTT(pNTT, pNTT)
		ringQ.NTT(qNTT, qNTT)

		exactMulNTT := ringQ.NewPoly()
		ringQ.MulCoeffsMontgomery(pNTT, qNTT, exactMulNTT)

		exactMulCoe := ringQ.NewPoly()
		ringQ.InvNTT(exactMulNTT, exactMulCoe)

		//------------------------------------------------------------
		// 3) “FFT<Big>” polynomial multiplication:
		//    aEval = ComplexEvaluate(p), bEval = ComplexEvaluate(q)
		//    prodEval = FieldMulBig(aEval, bEval)
		//    fftMulNTT = ComplexInterpolate(prodEval)
		//    InvNTT(fftMulNTT) = fftMulCoe
		//------------------------------------------------------------
		aEval := ComplexEvaluate(p, ringQ, prec)
		bEval := ComplexEvaluate(q, ringQ, prec)
		aEval.Domain = Eval
		bEval.Domain = Eval

		prodEval := FieldMulBig(aEval, bEval) // still in Eval
		fftMulNTT := ComplexInterpolate(prodEval, ringQ)

		fftMulCoe := ringQ.NewPoly()
		ringQ.InvNTT(fftMulNTT, fftMulCoe)

		//------------------------------------------------------------
		// 4) Compare exactMulCoe vs fftMulCoe
		//------------------------------------------------------------
		mismatchMul := 0
		for i := 0; i < N; i++ {
			if exactMulCoe.Coeffs[0][i] != fftMulCoe.Coeffs[0][i] {
				mismatchMul++
			}
		}
		if mismatchMul > 0 {
			t.Errorf("trial %d: FFT‐based multiply mismatches at %d coefficient(s)", trial, mismatchMul)
		} else {
			t.Logf("trial %d: FFT‐based multiply OK", trial)
		}

		//------------------------------------------------------------
		//    exactAddCoe = p + q (coefficient-wise mod q)
		//    fftAddNTT = ComplexInterpolate(FieldAddBig(aEval, bEval))
		//    fftAddCoe = InvNTT(fftAddNTT)
		//------------------------------------------------------------
		exactAddCoe := ringQ.NewPoly()
		ringQ.Add(p, q, exactAddCoe)

		sumEval := FieldAddBig(aEval, bEval)
		sumNTT := ComplexInterpolate(sumEval, ringQ)
		sumCoe := ringQ.NewPoly()
		ringQ.InvNTT(sumNTT, sumCoe)

		mismatchAdd := 0
		for i := 0; i < N; i++ {
			if exactAddCoe.Coeffs[0][i] != sumCoe.Coeffs[0][i] {
				mismatchAdd++
			}
		}
		if mismatchAdd > 0 {
			t.Errorf("trial %d: FFT‐based add mismatches at %d coefficient(s)", trial, mismatchAdd)
		} else {
			t.Logf("trial %d: FFT‐based add OK", trial)
		}
	}
}

// naiveMul computes (a·b) mod (x^n+1, q) by direct convolution:
//
//	c_k = sum_{i+j = k} a_i b_j  −  sum_{i+j = k+n} a_i b_j   (all mod q)
func naiveMul(a, b *ring.Poly, ringQ *ring.Ring) *ring.Poly {
	n := ringQ.N
	q := int64(ringQ.Modulus[0])
	c := ringQ.NewPoly()
	for k := 0; k < n; k++ {
		var sum int64 = 0
		for i := 0; i < n; i++ {
			j := k - i
			sign := int64(1)
			if j < 0 {
				j += n
				sign = -1
			}
			ai := int64(a.Coeffs[0][i])
			bj := int64(b.Coeffs[0][j])
			sum += sign * ai * bj
		}
		sum %= q
		if sum < 0 {
			sum += q
		}
		c.Coeffs[0][k] = uint64(sum)
	}
	return c
}

// TestFFTvsExactMulSmall checks that “naïveMul” and the FFT-based approach agree
// exactly in R = Z_97[x]/(x^16+1).  We must use ComplexEvaluateSub / ComplexInterpolateSub
// because our ring is modulo x^n + 1, not x^n − 1.
func TestFFTvsExactMulSmall(t *testing.T) {
	const (
		n16 = 16
		q97 = uint64(97)
	)

	// Build ℤ_97[x]/(x^16+1).
	ringR, err := ring.NewRing(n16, []uint64{q97})
	if err != nil {
		t.Fatalf("ring.NewRing(%d, [%d]) failed: %v", n16, q97, err)
	}

	prec := uint(64)
	rand.Seed(time.Now().UnixNano())

	const trials = 20
	for trial := 0; trial < trials; trial++ {
		// 1) pick two random coefficient‐polynomials a,b of length 16
		a := ringR.NewPoly()
		b := ringR.NewPoly()
		for i := 0; i < n16; i++ {
			a.Coeffs[0][i] = uint64(rand.Intn(int(q97)))
			b.Coeffs[0][i] = uint64(rand.Intn(int(q97)))
		}

		// 2) compute cExact = naiveMul(a,b)
		cExact := naiveMul(a, b, ringR)

		// 3) compute cFFT by going “(x^n+1)→ evaluation at 2n-th roots→ interp”
		//   3a) aEval = ComplexEvaluateSub(a,16,ringR,prec)
		aEval := ComplexEvaluateSub(a, n16, ringR, prec)
		bEval := ComplexEvaluateSub(b, n16, ringR, prec)
		aEval.Domain = Eval
		bEval.Domain = Eval

		//   3b) pointwise multiply in K_{2n}:
		prodEval := NewFieldElemBig(n16, prec)
		for i := 0; i < n16; i++ {
			prodEval.Coeffs[i] = aEval.Coeffs[i].Mul(bEval.Coeffs[i])
		}
		prodEval.Domain = Eval

		//   3c) ComplexInterpolateSub(prodeval,16,ringR) → a *ring.Poly* in coeff form
		cFFT := ComplexInterpolateSub(prodEval, n16, ringR)

		// 4) compare cExact vs cFFT, coefficient‐by‐coefficient
		var mismatches int
		for i := 0; i < n16; i++ {
			e := cExact.Coeffs[0][i]
			f := cFFT.Coeffs[0][i]
			if e != f {
				mismatches++
				if mismatches <= 5 {
					t.Logf("trial %d, idx %d: naive=%d, FFT=%d", trial, i, e, f)
				}
			}
		}
		if mismatches != 0 {
			t.Fatalf("trial %d: found %d mismatches in n=16, q=97", trial, mismatches)
		}
	}
}
