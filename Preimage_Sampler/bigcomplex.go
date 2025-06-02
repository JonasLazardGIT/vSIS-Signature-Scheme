package Preimage_Sampler

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"math/rand"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Domain indicates whether a CyclotomicFieldElem is in coefficient or evaluation domain.
type Domain int

const (
	Coeff Domain = iota
	Eval
)

// BigComplex represents a complex number with arbitrary-precision parts.
type BigComplex struct {
	Real *big.Float
	Imag *big.Float
}

// NewBigComplex creates a BigComplex with given real, imag and precision.
func NewBigComplex(real, imag float64, prec uint) *BigComplex {
	return &BigComplex{
		Real: new(big.Float).SetPrec(prec).SetFloat64(real),
		Imag: new(big.Float).SetPrec(prec).SetFloat64(imag),
	}
}

// NewBigComplexFromFloat copies floats into a BigComplex.
func NewBigComplexFromFloat(re, im *big.Float) *BigComplex {
	return &BigComplex{
		Real: new(big.Float).Copy(re),
		Imag: new(big.Float).Copy(im),
	}
}

// NewBigComplexZero returns zero BigComplex at precision.
func NewBigComplexZero(prec uint) *BigComplex {
	return NewBigComplex(0, 0, prec)
}

// Add returns z + w.
func (z *BigComplex) Add(w *BigComplex) *BigComplex {
	return &BigComplex{
		Real: new(big.Float).Add(z.Real, w.Real),
		Imag: new(big.Float).Add(z.Imag, w.Imag),
	}
}

// Sub returns z - w.
func (z *BigComplex) Sub(w *BigComplex) *BigComplex {
	return &BigComplex{
		Real: new(big.Float).Sub(z.Real, w.Real),
		Imag: new(big.Float).Sub(z.Imag, w.Imag),
	}
}

// Mul returns z * w.
func (z *BigComplex) Mul(w *BigComplex) *BigComplex {
	ac := new(big.Float).Mul(z.Real, w.Real)
	bd := new(big.Float).Mul(z.Imag, w.Imag)
	ad := new(big.Float).Mul(z.Real, w.Imag)
	bc := new(big.Float).Mul(z.Imag, w.Real)
	return &BigComplex{
		Real: new(big.Float).Sub(ac, bd),
		Imag: new(big.Float).Add(ad, bc),
	}
}

// Conj returns complex conjugate.
func (z *BigComplex) Conj() *BigComplex {
	return &BigComplex{
		Real: new(big.Float).Copy(z.Real),
		Imag: new(big.Float).Neg(z.Imag),
	}
}

// AbsSquared returns |z|^2.
func (z *BigComplex) AbsSquared() *big.Float {
	r2 := new(big.Float).Mul(z.Real, z.Real)
	i2 := new(big.Float).Mul(z.Imag, z.Imag)
	return new(big.Float).Add(r2, i2)
}

// Inv returns 1/z.
func (z *BigComplex) Inv() *BigComplex {
	conj := z.Conj()
	d := z.AbsSquared()
	return &BigComplex{
		Real: new(big.Float).Quo(conj.Real, d),
		Imag: new(big.Float).Quo(conj.Imag, d),
	}
}

// Copy returns a deep copy.
func (z *BigComplex) Copy() *BigComplex {
	return &BigComplex{
		Real: new(big.Float).Copy(z.Real),
		Imag: new(big.Float).Copy(z.Imag),
	}
}

// DivBy divides z by scalar.
func (z *BigComplex) DivBy(scalar *big.Float) *BigComplex {
	return &BigComplex{
		Real: new(big.Float).Quo(z.Real, scalar),
		Imag: new(big.Float).Quo(z.Imag, scalar),
	}
}

// CyclotomicFieldElem is an element in K_{2N}.
type CyclotomicFieldElem struct {
	N      int
	Coeffs []*BigComplex
	Domain Domain
}

// NewFieldElemBig allocates a zero field element in coeff domain.
func NewFieldElemBig(n int, prec uint) *CyclotomicFieldElem {
	coeffs := make([]*BigComplex, n)
	for i := range coeffs {
		coeffs[i] = NewBigComplexZero(prec)
	}
	return &CyclotomicFieldElem{N: n, Coeffs: coeffs, Domain: Coeff}
}

// DggSampler interface for Gaussian samplers.
type DggSampler interface {
	SampleKarney(mean, stddev float64) int64
}

// Matrix is a generic matrix.
type Matrix[T any] struct {
	Rows, Cols int
	Data       []T
}

// NewMatrix allocates a matrix of zero values.
func NewMatrix[T any](rows, cols int) *Matrix[T] {
	return &Matrix[T]{Rows: rows, Cols: cols, Data: make([]T, rows*cols)}
}

// At returns element (i,j).
func (m *Matrix[T]) At(i, j int) T {
	return m.Data[i*m.Cols+j]
}

// Set assigns element (i,j) = v.
func (m *Matrix[T]) Set(i, j int, v T) {
	m.Data[i*m.Cols+j] = v
}

// ----------------------------------------------------------------
// FFTBig: take “coeffs” (length n, a power of two) in coefficient form and
//
//	return the length-n slice of *BigComplex values evaluated at the
//	primitive 2n-th roots of unity e^{-2πi·k/size} (forward FFT).
//
//	- coeffs: slice of length n of *BigComplex (domain: “time” or “coeff”).
//	- prec: the precision (in bits) that should be used for every BigComplex
//	        in the output and intermediate twiddle factors.
//	- Returns a brand-new slice of length n of *BigComplex in the “evaluation” domain.
//
//	**IMPORTANT**: internally we compute the complex exponentials e^{-2πi/size}
//	by calling math.Cos, math.Sin in float64 and then lifting to big.Float.
//	If you need **full** big‐float accuracy on those sines/cosines, you must
//	replace the calls to math.Cos/math.Sin by a high‐precision big‐float routine.
//
//	NOTE: this is the standard in-place Cooley–Tukey; we allocate a new buffer
//	so that we do not destroy the input argument.
//
// ----------------------------------------------------------------
func FFTBig(coeffs []*BigComplex, prec uint) []*BigComplex {
	n := len(coeffs)
	if n == 0 || (n&(n-1)) != 0 {
		panic("FFTBig: length must be a nonzero power of 2")
	}

	// 1) Copy input into result array
	result := make([]*BigComplex, n)
	for i := 0; i < n; i++ {
		result[i] = coeffs[i].Copy()
	}

	// 2) Bit-reversal reordering
	logN := bits.Len(uint(n)) - 1
	for i := 0; i < n; i++ {
		j := bitReverseBig(i, logN)
		if i < j {
			result[i], result[j] = result[j], result[i]
		}
	}

	// 3) Iterative Cooley–Tukey
	//    “size” loops over 2, 4, 8, …, n
	for size := 2; size <= n; size <<= 1 {
		half := size >> 1

		// angle = −2π/size  (forward FFT => negative)
		// we compute cos(angle) and sin(angle) as float64, then lift to big.Float
		angleF := -2.0 * math.Pi / float64(size)
		cosF := big.NewFloat(0).SetPrec(prec).SetFloat64(math.Cos(angleF))
		sinF := big.NewFloat(0).SetPrec(prec).SetFloat64(math.Sin(angleF))

		// wn = e^{i * angle} in *BigComplex form
		wn := &BigComplex{
			Real: new(big.Float).Copy(cosF),
			Imag: new(big.Float).Copy(sinF),
		}

		for start := 0; start < n; start += size {
			// w = 1 + 0i (as *BigComplex)
			w := &BigComplex{
				Real: big.NewFloat(1).SetPrec(prec),
				Imag: big.NewFloat(0).SetPrec(prec),
			}

			for j := 0; j < half; j++ {
				idx1 := start + j
				idx2 := start + j + half

				// temp = w * result[idx2]
				temp := result[idx2].Mul(w)

				// result[idx2] = result[idx1] - temp
				result[idx2] = result[idx1].Sub(temp)

				// result[idx1] = result[idx1] + temp
				result[idx1] = result[idx1].Add(temp)

				// w = w * wn   (update twiddle)
				w = w.Mul(wn)
			}
		}
	}

	return result
}

// ----------------------------------------------------------------
// IFFTBig: take a slice of length n = power‐of‐two in “evaluation” domain
//
//	(i.e. values at the 2n-th roots e^{-2πik/size}), and return the inverse FFT
//	so that the output is back in “coefficient” form. We scale by 1/n at the end.
//
//	- evals:   slice of length n of *BigComplex in evaluation domain.
//	- prec:    desired precision (bits) for all intermediate big‐floats.
//	- Returns length-n slice of *BigComplex in coefficient domain.
//
//	Again, we do exactly the same bit-reversal + Cooley–Tukey, but we use
//	the opposite sign of the root (angle = +2π/size), and in the end divide
//	every coordinate by n (via big-float division).
//
// ----------------------------------------------------------------
func IFFTBig(evals []*BigComplex, prec uint) []*BigComplex {
	n := len(evals)
	if n == 0 || (n&(n-1)) != 0 {
		panic("IFFTBig: length must be a nonzero power of 2")
	}

	// 1) Copy evals into result
	result := make([]*BigComplex, n)
	for i := 0; i < n; i++ {
		result[i] = evals[i].Copy()
	}

	// 2) Bit-reversal
	logN := bits.Len(uint(n)) - 1
	for i := 0; i < n; i++ {
		j := bitReverseBig(i, logN)
		if i < j {
			result[i], result[j] = result[j], result[i]
		}
	}

	// 3) Inverse Cooley–Tukey (angle = +2π/size)
	for size := 2; size <= n; size <<= 1 {
		half := size >> 1

		angleF := 2.0 * math.Pi / float64(size) // **positive** for inverse FFT
		cosF := big.NewFloat(0).SetPrec(prec).SetFloat64(math.Cos(angleF))
		sinF := big.NewFloat(0).SetPrec(prec).SetFloat64(math.Sin(angleF))

		wn := &BigComplex{
			Real: new(big.Float).Copy(cosF),
			Imag: new(big.Float).Copy(sinF),
		}

		for start := 0; start < n; start += size {
			w := &BigComplex{
				Real: big.NewFloat(1).SetPrec(prec),
				Imag: big.NewFloat(0).SetPrec(prec),
			}
			for j := 0; j < half; j++ {
				idx1 := start + j
				idx2 := start + j + half

				temp := result[idx2].Mul(w)

				result[idx2] = result[idx1].Sub(temp)
				result[idx1] = result[idx1].Add(temp)

				w = w.Mul(wn)
			}
		}
	}

	// 4) Divide everything by n (scale by 1/n in big-float)
	bigN := big.NewFloat(0).SetPrec(prec).SetFloat64(float64(n))
	invN := new(big.Float).SetPrec(prec).Quo(big.NewFloat(1).SetPrec(prec), bigN)

	for i := 0; i < n; i++ {
		result[i].Real = result[i].Real.Mul(result[i].Real, invN)
		result[i].Imag = result[i].Imag.Mul(result[i].Imag, invN)
	}

	return result
}

// ----------------------------------------------------------------
// bitReverseBig: compute the bit-reversal of “i” in a log2(n)-bit space.
//
//	E.g. if n=8 (log2(n)=3), then bitReverseBig(3,3) = 6 because (011)_2 → (110)_2.
//
// ----------------------------------------------------------------
func bitReverseBig(i, logN int) int {
	var rev int
	for b := 0; b < logN; b++ {
		rev = (rev << 1) | ((i >> b) & 1)
	}
	return rev
}

// --------------------------------------------------------------------------------
//
//	ConvertFromPolyBig
//
//	Given a ring.Ring “r” and a *ring.Poly “p” (in coefficient form modulo q),
//	produce a *CyclotomicFieldElem of length n=r.N in the EVAL (FFT) domain.
//	We do exactly:
//
//	    1) Read p.Coeffs[0][i] as an integer mod q, convert to float64.
//	    2) Lift to a BigComplex with (real=float64, imag=0), to high precision.
//	    3) Call FFTBig(...) on that slice of length n.
//	    4) Store the resulting *BigComplex directly into CyclotomicFieldElem.Coeffs.
//
//	Note: we do not round to float64 until *after* the FFT is done in high precision.
//
// --------------------------------------------------------------------------------
func ConvertFromPolyBig(r *ring.Ring, p *ring.Poly, prec uint) *CyclotomicFieldElem {
	n := r.N

	// 1) Build a slice of length n of *BigComplex from p.Coeffs:
	bigSlice := make([]*BigComplex, n)
	for i := 0; i < n; i++ {
		// Convert p.Coeffs[0][i] (a uint64) into a float64 in [0, q).
		f64 := ModQToFloat64(p.Coeffs[0][i], r.Modulus[0])
		// Create a BigComplex(re=f64, im=0) at precision “prec”:
		re := new(big.Float).SetPrec(prec).SetFloat64(f64)
		im := new(big.Float).SetPrec(prec).SetFloat64(0.0)
		bigSlice[i] = NewBigComplexFromFloat(re, im)
	}

	// 2) Run the high‐precision FFT:
	fftResult := FFTBig(bigSlice, prec)

	// 3) Copy the result into a new CyclotomicFieldElem (in Eval domain):
	out := NewFieldElemBig(n, prec)
	out.Domain = Eval
	for i := 0; i < n; i++ {
		out.Coeffs[i] = fftResult[i].Copy()
	}
	return out
}

// --------------------------------------------------------------------------------
//
//	ConvertToPolyBig
//
//	Given a *CyclotomicFieldElem “f” (in Eval domain), interpolate it back into
//	a *ring.Poly (coeff domain), then reduce mod q.  Steps:
//
//	  1) Take f.Coeffs[i] (each is *BigComplex), build a slice []*BigComplex.
//	  2) Call IFFTBig(...) to recover the length‐n coefficient vector (high precision).
//	  3) For each output index i: take (Real part, a *big.Float), convert to float64
//	     (via .Float64), round to nearest int64, and reduce mod q.
//
// --------------------------------------------------------------------------------
func ConvertToPolyBig(f *CyclotomicFieldElem, r *ring.Ring) *ring.Poly {
	n := f.N

	// 1) Build a slice of length n of *BigComplex (the “evaluation” values).
	evalSlice := make([]*BigComplex, n)
	for i := 0; i < n; i++ {
		evalSlice[i] = f.Coeffs[i].Copy()
	}

	// 2) Call IFFTBig at the same precision
	prec := f.Coeffs[0].Real.Prec()
	coeffSlice := IFFTBig(evalSlice, prec)

	// 3) Create a new ring.Poly and round each coefficient
	P := r.NewPoly()
	for i := 0; i < n && i < r.N; i++ {
		// coeffSlice[i].Real is a high‐precision big.Float
		// Convert to float64 and round to nearest integer:
		realBF := coeffSlice[i].Real
		f64, _ := realBF.Float64()   // lose some bits here, but f64 should be within ±2^53
		ri := int64(math.Round(f64)) // nearest integer
		// reduce mod q
		qi := int64(r.Modulus[0])
		ri = ((ri % qi) + qi) % qi
		P.Coeffs[0][i] = uint64(ri)
	}
	return P
}

// --------------------------------------------------------------------------------
//
//	SampleGaussianFieldElemBig
//
//	“f” is an evaluation‐domain element of length n (FFT domain).  We wish to sample
//	a cyclotomic field element whose coordinates (in Eval domain) are i.i.d.  N(0,σ²).
//	We can do that coordinate‐wise by drawing two independent Normal()’s in float64,
//	then lifting each into a BigComplex at precision “prec.”  (That is exactly what
//	your old code was doing, so no change is strictly required here.)
//
// --------------------------------------------------------------------------------
func SampleGaussianFieldElemBig(n int, sigma float64, prec uint) *CyclotomicFieldElem {
	out := NewFieldElemBig(n, prec)
	out.Domain = Eval
	for i := 0; i < n; i++ {
		re64 := rand.NormFloat64() * sigma
		im64 := rand.NormFloat64() * sigma
		re := new(big.Float).SetPrec(prec).SetFloat64(re64)
		im := new(big.Float).SetPrec(prec).SetFloat64(im64)
		out.Coeffs[i] = NewBigComplexFromFloat(re, im)
	}
	return out
}

// ----------------------------------------------------------------
// NegacyclicEvaluatePoly
//
// Given:
//   - p      : *ring.Poly in “coefficient” form (length = m = ringQ.N, mod q).
//   - ringQ  : *ring.Ring which knows m and the modulus.
//   - prec   : desired bit-precision for all BigComplex arithmetic.
//
// Returns:
//   - *CyclotomicFieldElem of length m, containing p evaluated at the
//     2m-th roots of unity “ω^{2k+1}” (i.e. a negacyclic FFT).
//
// How it works:
//  1. Zero-pad to length 2m.
//  2. Call FFTBig on that 2m-length buffer.
//  3. Extract only the odd indices (2k+1), producing an m-length slice.
//  4. Copy into a *CyclotomicFieldElem with Domain=Eval.
//
// ----------------------------------------------------------------
func NegacyclicEvaluatePoly(p *ring.Poly, ringQ *ring.Ring, prec uint) *CyclotomicFieldElem {
	// m := ring dimension = ringQ.N
	m := ringQ.N
	twoM := 2 * m

	// 1) Build a length-2m slice of *BigComplex, zero-padding after index m-1
	//    We interpret each coefficient p.Coeffs[0][i] as a signed integer in [−q/2..+q/2).

	// allocate and fill the first m entries
	A := make([]*BigComplex, twoM)
	for i := 0; i < m; i++ {
		// Convert p.Coeffs[0][i] to signed int64 in [−q/2..+q/2)
		// signed := UnsignedToSigned(p.Coeffs[0][i], ringQ.Modulus[0])
		// lift to BigComplex
		reBF := new(big.Float).SetPrec(prec).SetFloat64(float64(p.Coeffs[0][i])) //!signed
		imBF := new(big.Float).SetPrec(prec).SetFloat64(0.0)
		A[i] = &BigComplex{Real: reBF, Imag: imBF}
	}
	// fill indices m..2m−1 with exact zero
	zeroBF := new(big.Float).SetPrec(prec).SetFloat64(0.0)
	zeroBC := &BigComplex{Real: zeroBF, Imag: zeroBF}
	for i := m; i < twoM; i++ {
		A[i] = zeroBC
	}

	// 2) perform length-2m forward FFTBig (angle = −2π/(2m))
	B := FFTBig(A, prec) // B has length = 2m

	// 3) extract only the odd indices B[2k+1], k=0..m−1
	evals := make([]*BigComplex, m)
	for k := 0; k < m; k++ {
		evals[k] = B[2*k+1].Copy()
	}

	// 4) pack into a CyclotomicFieldElem in “evaluation” domain
	out := NewFieldElemBig(m, prec)
	out.Domain = Eval
	for k := 0; k < m; k++ {
		out.Coeffs[k] = evals[k]
	}
	return out
}

// ----------------------------------------------------------------
// NegacyclicInterpolateElem
//
// Given:
//   - f      : *CyclotomicFieldElem of length m, Domain=Eval (i.e. f[k] = p(ω^{2k+1})),
//     where ω = e^{2πi/(2m)}.
//   - ringQ  : *ring.Ring which knows m and the modulus.
//   - prec   : bit-precision used when f.Coeffs[*] were constructed.
//
// Returns:
//   - *ring.Poly (length = m, in coefficient form) corresponding to p(x) mod (x^m + 1).
//
// How it works:
//  1. Form a length-2m slice by inserting f.Coeffs[k] into index 2k+1, zeros at evens.
//  2. Call IFFTBig on that length-2m buffer (angle = +2π/(2m), then divide by 2m).
//  3. Take the first m BigComplex results, round Real→int, reduce mod q.
//
// ----------------------------------------------------------------
func NegacyclicInterpolateElem(f *CyclotomicFieldElem, ringQ *ring.Ring) *ring.Poly {
	// m := ring dimension = ringQ.N
	m := ringQ.N
	twoM := 2 * m

	// 1) form length-2m slice A such that A[2k+1] = f.Coeffs[k], A[even]=0
	A := make([]*BigComplex, twoM)

	prec := f.Coeffs[0].Real.Prec()
	zeroBF := new(big.Float).SetPrec(prec).SetFloat64(0.0)
	zeroBC := &BigComplex{Real: zeroBF, Imag: zeroBF}
	// fill evens with zero
	for i := 0; i < twoM; i += 2 {
		A[i] = zeroBC
	}
	// copy f.Coeffs[k] into odd slots
	for k := 0; k < m; k++ {
		A[2*k+1] = f.Coeffs[k].Copy()
	}

	// 2) perform length-2m inverse FFTBig (angle = +2π/(2m), then scale by 1/(2m))
	inv := IFFTBig(A, prec) // inv has length = 2m

	P := ringQ.NewPoly()
	q := int64(ringQ.Modulus[0])
	for j := 0; j < m; j++ {
		realBF := inv[j].Real
		f64, _ := realBF.Float64()

		// Multiply by 2 to undo the “÷2” that IFFTBig implicitly introduced ★
		rInt := int64(math.Round(f64 * 2.0))

		// reduce mod q into [0..q−1]
		rInt = ((rInt % q) + q) % q
		P.Coeffs[0][j] = uint64(rInt)
	}

	return P
}

// ToEvalNegacyclic converts a CyclotomicFieldElem that is currently
// in COEFF domain into the corresponding negacyclic‐FFT “Eval” form.
//   - e:      input in COEFF domain (BigComplex entries represent integers mod q).
//   - ringQ:  the ring.Ring (defines n and Modulus).
//   - prec:   desired precision (in bits) for BigFloat twiddles.
//
// Returns a fresh *CyclotomicFieldElem in Eval domain.
// Panics if e.Domain != Coeff.
func ToEvalNegacyclic(e *CyclotomicFieldElem, ringQ *ring.Ring, prec uint) *CyclotomicFieldElem {
	if e.Domain != Coeff {
		panic("ToEvalNegacyclic: input must be in Coeff domain")
	}
	n := e.N
	// 1) Build a ring.Poly by rounding each BigComplex coefficient → uint64 mod q.
	P := ringQ.NewPoly()
	q := int64(ringQ.Modulus[0])
	for i := 0; i < n && i < ringQ.N; i++ {
		// e.Coeffs[i].Real is a *big.Float that should be very close to an integer.
		realBF := e.Coeffs[i].Real
		f64, _ := realBF.Float64()
		rInt := int64(math.Round(f64))
		// reduce mod q
		rInt = ((rInt % q) + q) % q
		P.Coeffs[0][i] = uint64(rInt)
	}
	// 2) Call the negacyclic evaluator (returns an element in Eval domain).
	out := NegacyclicEvaluatePoly(P, ringQ, prec)
	// out.Domain will be Eval already, and out.Coeffs[i] is BigComplex(evaluated at e^{-iπj/n}).
	return out
}

// ToCoeffNegacyclic converts a CyclotomicFieldElem that is currently
// in negacyclic‐Eval domain back into a CyclotomicFieldElem in COEFF domain.
//   - e:      input in Eval domain (BigComplex = values at 2n-roots of −1).
//   - ringQ:  the ring.Ring (defines n and Modulus).
//   - prec:   desired precision (in bits) for the output BigComplex coefficients.
//
// Returns a fresh *CyclotomicFieldElem in Coeff domain.
// Panics if e.Domain != Eval.
func ToCoeffNegacyclic(e *CyclotomicFieldElem, ringQ *ring.Ring, prec uint) *CyclotomicFieldElem {
	if e.Domain != Eval {
		panic("ToCoeffNegacyclic: input must be in Eval domain")
	}
	n := e.N
	// 1) Inverse‐transform via NegacyclicInterpolateElem → a ring.Poly in coefficient form.
	P := NegacyclicInterpolateElem(e, ringQ)
	// 2) Copy P’s uint64 coefficients into a new CyclotomicFieldElem (in Coeff domain).
	out := NewFieldElemBig(n, prec)
	for i := 0; i < n && i < ringQ.N; i++ {
		// lift integer P.Coeffs[0][i] → BigFloat
		reBF := new(big.Float).SetPrec(prec).SetFloat64(float64(P.Coeffs[0][i]))
		imBF := new(big.Float).SetPrec(prec).SetFloat64(0.0)
		out.Coeffs[i] = NewBigComplexFromFloat(reBF, imBF)
	}
	out.Domain = Coeff
	return out
}

// Field operations
func FieldAddBig(a, b *CyclotomicFieldElem) *CyclotomicFieldElem {
	if a.N != b.N {
		panic("FieldAddBig: dimension mismatch")
	}
	res := NewFieldElemBig(a.N, a.Coeffs[0].Real.Prec())
	for i := 0; i < a.N; i++ {
		res.Coeffs[i] = a.Coeffs[i].Add(b.Coeffs[i])
	}
	return res
}

func FieldSubBig(a, b *CyclotomicFieldElem) *CyclotomicFieldElem {
	if a.N != b.N {
		panic("FieldSubBig: dimension mismatch")
	}
	res := NewFieldElemBig(a.N, a.Coeffs[0].Real.Prec())
	for i := 0; i < a.N; i++ {
		res.Coeffs[i] = a.Coeffs[i].Sub(b.Coeffs[i])
	}
	return res
}

func FieldMulBig(a, b *CyclotomicFieldElem) *CyclotomicFieldElem {
	if a.N != b.N {
		panic("FieldMulBig: dimension mismatch")
	}
	res := NewFieldElemBig(a.N, a.Coeffs[0].Real.Prec())
	for i := 0; i < a.N; i++ {
		res.Coeffs[i] = a.Coeffs[i].Mul(b.Coeffs[i])
	}
	return res
}

// Field scalar multiplication and division
type FieldScalar struct{ Val *BigComplex }

func FieldScalarMulBig(a *CyclotomicFieldElem, s *BigComplex) *CyclotomicFieldElem {
	res := NewFieldElemBig(a.N, a.Coeffs[0].Real.Prec())
	for i := 0; i < a.N; i++ {
		res.Coeffs[i] = a.Coeffs[i].Mul(s)
	}
	return res
}

func FieldScalarDiv(a *CyclotomicFieldElem, norm []*big.Float) *CyclotomicFieldElem {
	if len(norm) != a.N {
		panic("FieldScalarDiv: length mismatch")
	}
	res := NewFieldElemBig(a.N, a.Coeffs[0].Real.Prec())
	for i := 0; i < a.N; i++ {
		if norm[i].Cmp(big.NewFloat(0).SetPrec(a.Coeffs[0].Real.Prec())) == 0 {
			panic(fmt.Sprintf("division by zero norm at %d", i))
		}
		res.Coeffs[i] = a.Coeffs[i].DivBy(norm[i])
	}
	return res
}

// HermitianTransposeFieldElem computes the “polynomial transpose” f ↦ fᵗ.
//   - If f.Domain == Coeff:  fᵗ(x) = f₀ − fₙ₋₁·x − fₙ₋₂·x² − … − f₁·xⁿ⁻¹
//   - If f.Domain == Eval:   fᵗ(ζ²ⁿᵏ⁺¹) = f(ζ²ⁿ⁻(²ⁿᵏ⁺¹)) = f(ζ²ⁿ⁻(²ᵏ⁺¹)) = f(ζ²(n−k−1)+1),
//     so out.Eval[k] = in.Eval[n−k−1].
func HermitianTransposeFieldElem(f *CyclotomicFieldElem) *CyclotomicFieldElem {
	n := f.N
	prec := f.Coeffs[0].Real.Prec()
	out := NewFieldElemBig(n, prec)

	switch f.Domain {
	case Coeff:
		// Coefficient‐space transpose:
		//   out[0] = f[0],
		//   out[i] = − f[n−i],   for i = 1..n−1
		for i := 0; i < n; i++ {
			if i == 0 {
				// copy f.Coeffs[0] exactly
				out.Coeffs[0] = f.Coeffs[0].Copy()
			} else {
				// take f.Coeffs[n−i], multiply by −1
				src := f.Coeffs[n-i]
				negReal := new(big.Float).SetPrec(prec).Neg(src.Real)
				negImag := new(big.Float).SetPrec(prec).Neg(src.Imag)
				out.Coeffs[i] = &BigComplex{
					Real: negReal,
					Imag: negImag,
				}
			}
		}
		out.Domain = Coeff

	case Eval:
		// Evaluation‐space transpose = “automorphism by 2n−1”:
		//   out.Eval[k] = f.Eval[n−k−1].  (no extra conjugation here)
		for k := 0; k < n; k++ {
			// f.Coeffs[n−k−1] is already a *BigComplex; copy it
			out.Coeffs[k] = f.Coeffs[n-k-1].Copy()
		}
		out.Domain = Eval

	default:
		panic("HermitianTransposeFieldElem: unknown domain")
	}
	return out
}

// Field inverse with norms.
func FieldInverseDiagWithNorm(d *CyclotomicFieldElem) (*CyclotomicFieldElem, []*big.Float) {
	n := d.N
	prec := d.Coeffs[0].Real.Prec()
	inv := NewFieldElemBig(n, prec)
	norms := make([]*big.Float, n)
	zero := new(big.Float).SetPrec(prec).SetFloat64(0)
	for i := 0; i < n; i++ {
		re := d.Coeffs[i].Real
		im := d.Coeffs[i].Imag
		re2 := new(big.Float).Mul(re, re)
		im2 := new(big.Float).Mul(im, im)
		n := new(big.Float).Add(re2, im2)
		if n.Cmp(zero) == 0 {
			panic(fmt.Sprintf("zero norm at %d", i))
		}
		conj := &BigComplex{Real: new(big.Float).Copy(re), Imag: new(big.Float).Neg(im)}
		inv.Coeffs[i] = conj
		norms[i] = n
	}
	return inv, norms
}

// PstrideBig splits even and odd coefficients.
func PstrideBig(c *CyclotomicFieldElem) (*CyclotomicFieldElem, *CyclotomicFieldElem) {
	n := c.N
	h := n / 2
	prec := c.Coeffs[0].Real.Prec()
	c0 := NewFieldElemBig(h, prec)
	c1 := NewFieldElemBig(h, prec)
	for i := 0; i < h; i++ {
		c0.Coeffs[i] = c.Coeffs[2*i]
		c1.Coeffs[i] = c.Coeffs[2*i+1]
	}
	return c0, c1
}

// SubScalar subtracts scalar from all coords.
func (f *CyclotomicFieldElem) SubScalar(c *BigComplex) {
	for i := range f.Coeffs {
		f.Coeffs[i] = f.Coeffs[i].Sub(c)
	}
}

// AddScalar adds scalar to all coords.
func (f *CyclotomicFieldElem) AddScalar(c *BigComplex) {
	for i := range f.Coeffs {
		f.Coeffs[i] = f.Coeffs[i].Add(c)
	}
}

// Copy returns a deep copy of the CyclotomicFieldElem, including all BigComplex entries
func (f *CyclotomicFieldElem) Copy() *CyclotomicFieldElem {
	// Determine precision from an existing coefficient
	prec := f.Coeffs[0].Real.Prec()

	// Allocate a new element with the same length and precision
	out := NewFieldElemBig(f.N, prec)
	out.Domain = f.Domain

	// Deep‐copy each BigComplex
	for i := 0; i < f.N; i++ {
		out.Coeffs[i] = f.Coeffs[i].Copy()
	}

	return out
}

// Conj returns the coordinate‐wise complex conjugate of a CyclotomicFieldElem.
func (f *CyclotomicFieldElem) Conj() *CyclotomicFieldElem {
	prec := f.Coeffs[0].Real.Prec()
	out := NewFieldElemBig(f.N, prec)
	out.Domain = f.Domain
	for i := 0; i < f.N; i++ {
		out.Coeffs[i] = f.Coeffs[i].Conj()
	}
	return out
}

// ToComplex converts BigComplex to built-in complex128 (losing precision).
func (b *BigComplex) ToComplex() complex128 {
	r, _ := b.Real.Float64()
	i, _ := b.Imag.Float64()
	return complex(r, i)
}

// SetCoeffs copies the coefficients from src into e and returns e.
// Panics if dimensions do not match.
func (e *CyclotomicFieldElem) SetCoeffs(src *CyclotomicFieldElem) *CyclotomicFieldElem {
	if e.N != src.N {
		panic("SetCoeffs: dimension mismatch")
	}
	for i := 0; i < e.N; i++ {
		e.Coeffs[i] = src.Coeffs[i]
	}
	return e
}

// ExtractEven returns a new CyclotomicFieldElem containing the even-indexed coefficients of e.
func (e *CyclotomicFieldElem) ExtractEven() *CyclotomicFieldElem {
	m := e.N
	half := m / 2
	out := NewFieldElemBig(half, e.Coeffs[0].Real.Prec())
	for i := 0; i < half; i++ {
		out.Coeffs[i] = e.Coeffs[2*i]
	}
	return out
}

// ExtractOdd returns a new CyclotomicFieldElem containing the odd-indexed coefficients of e.
func (e *CyclotomicFieldElem) ExtractOdd() *CyclotomicFieldElem {
	m := e.N
	half := m / 2
	out := NewFieldElemBig(half, e.Coeffs[0].Real.Prec())
	for i := 0; i < half; i++ {
		out.Coeffs[i] = e.Coeffs[2*i+1]
	}
	return out
}

func InversePermuteFieldElem(x *CyclotomicFieldElem) {
	n := x.N
	tmp := make([]*BigComplex, n)
	half := n / 2
	e, o := 0, half
	// even positions come from the first half, odd from the second
	for i := 0; e < half; i += 2 {
		tmp[i] = x.Coeffs[e]
		tmp[i+1] = x.Coeffs[o]
		e++
		o++
	}
	copy(x.Coeffs, tmp)
}

// RingFreeToCoeffNegacyclic treats f ∈ Kₙ in EVAL domain (length=n, BigComplex values
// at the 2n-th roots of −1).  It returns a *CyclotomicFieldElem in COEFF domain by
// performing the inverse negacyclic FFT, all without ever using ring.Ring.
//
//   - n   = f.N
//   - q   = modulus
//   - prec = big-float precision.
func RingFreeToCoeffNegacyclic(
	f *CyclotomicFieldElem,
	n int,
	q uint64,
	prec uint,
) *CyclotomicFieldElem {
	if f.Domain != Eval {
		panic("RingFreeToCoeffNegacyclic: input must be in Eval domain")
	}
	m := n
	twoM := 2 * m

	// 1) Build length-2m slice A such that A[2k+1] = f.Coeffs[k], A[even] = 0
	A := make([]*BigComplex, twoM)
	zeroBF := new(big.Float).SetPrec(prec).SetFloat64(0)
	zeroBC := &BigComplex{Real: zeroBF, Imag: zeroBF}
	for i := 0; i < twoM; i += 2 {
		A[i] = zeroBC
	}
	for k := 0; k < m; k++ {
		A[2*k+1] = f.Coeffs[k].Copy()
	}

	// 2) Perform length-2m inverse FFTBig (angle = +2π/(2m)), then scale by 1/(2m).
	inv := IFFTBig(A, prec) // length=2m

	// 3) The “negacyclic IFFT” says: take the first m slots inv[0..m−1], round Real→int,
	//    multiply by 2 to reverse the “divide by 2” inside IFFTBig, then reduce mod q.
	out := NewFieldElemBig(m, prec)
	bigQ := int64(q)
	for j := 0; j < m; j++ {
		f64, _ := inv[j].Real.Float64()
		// IFFTBig already divided each coordinate by (2m).  Because
		//   IFFTBig does:   (   “regular inverse DFT”   )  / n   (n=2m),
		// our “true” negacyclic inverse wants 1/(2m) instead of 1/n=1/(2m).
		// In code above, your original interpo did “f64 * 2” to re-scale
		// from 1/(2m) to 1/m.  But the classical negacyclic inverse over
		// x^m+1 requires only a 1/m factor, not 1/(2m).  Since IFFTBig
		// used a 1/n=1/(2m) scaling, we must multiply by 2 to achieve
		// a net 1/m.  Hence the “f64*2.0” below.
		rInt := int64(math.Round(f64 * 2.0))
		// reduce mod q
		rInt = ((rInt % bigQ) + bigQ) % bigQ
		// lift back into *BigComplex with imag=0
		reBF := new(big.Float).SetPrec(prec).SetFloat64(float64(rInt))
		imBF := new(big.Float).SetPrec(prec).SetFloat64(0)
		out.Coeffs[j] = &BigComplex{Real: reBF, Imag: imBF}
	}
	out.Domain = Coeff
	return out
}

// RingFreeToEvalNegacyclic treats e ∈ Kₙ in COEFF domain (BigComplex entries ≈ integers mod q),
// and returns its negacyclic‐FFT in EVAL domain, all without ever constructing a ring.Ring.
//
//   - n   = e.N  (the length of the polynomial ring / cyclotomic degree).
//   - q   = the prime modulus (≡ 1 mod 2n).
//   - prec = big-float precision in bits.
//
// The result is a new *CyclotomicFieldElem of length n with Domain=Eval, whose
// .Coeffs[k] ≈ e(ω^{2k+1}) in high precision.
func RingFreeToEvalNegacyclic(
	e *CyclotomicFieldElem,
	n int,
	q uint64,
	prec uint,
) *CyclotomicFieldElem {
	if e.Domain != Coeff {
		panic("RingFreeToEvalNegacyclic: input must be in Coeff domain")
	}
	// 1) Build a “length-n” slice of uint64 coefficients in [0..q-1].
	//    We round each BigFloat→float64→int64, then reduce mod q.
	polyCoeffs := make([]uint64, n)
	bigQ := int64(q)
	for i := 0; i < n; i++ {
		// e.Coeffs[i].Real should be very close to an integer.
		f64, _ := e.Coeffs[i].Real.Float64()
		rInt := int64(math.Round(f64))
		// reduce into [0..q−1]
		rInt = ((rInt % bigQ) + bigQ) % bigQ
		polyCoeffs[i] = uint64(rInt)
	}

	// 2) Call our home-grown NegacyclicEvaluatePoly, except that it expects
	//    a *ring.Poly.  Instead, we can write a tiny adapter that treats
	//    “polyCoeffs” as if it were a ring.Poly.Coeffs[0][⋯].
	//
	//    In other words, we create a dummy *ring.Poly‐like wrapper that
	//    simply returns polyCoeffs[i] whenever NegacyclicEvaluatePoly asks
	//    for p.Coeffs[0][i].  To avoid pulling in ring.Ring entirely, one
	//    can copy NegacyclicEvaluatePoly’s logic here verbatim but replace
	//    “p.Coeffs[0][i]” by “polyCoeffs[i].”  For maximum clarity, we’ll
	//    just inline that part:

	// ———————————— Inlined “NegacyclicEvaluatePoly” but reading from polyCoeffs[] ————————————
	m := n
	twoM := 2 * m

	// Build length-2m slice A of *BigComplex, zero-padding after index m−1.
	A := make([]*BigComplex, twoM)
	for i := 0; i < m; i++ {
		// Interpret polyCoeffs[i] as an *unsigned* in [0..q), and
		// treat it directly as a BigComplex(re=that integer, im=0).
		reBF := new(big.Float).SetPrec(prec).SetFloat64(float64(polyCoeffs[i]))
		imBF := new(big.Float).SetPrec(prec).SetFloat64(0)
		A[i] = &BigComplex{Real: reBF, Imag: imBF}
	}
	zeroBF := new(big.Float).SetPrec(prec).SetFloat64(0)
	zeroBC := &BigComplex{Real: zeroBF, Imag: zeroBF}
	for i := m; i < twoM; i++ {
		A[i] = zeroBC
	}

	// 2b) Perform length-2m forward FFTBig (angle = −2π/(2m)).
	B := FFTBig(A, prec) // length=2m

	// 3) Extract only the odd indices B[2k+1], k=0..m−1
	evals := make([]*BigComplex, m)
	for k := 0; k < m; k++ {
		evals[k] = B[2*k+1].Copy()
	}

	// 4) Pack into a *CyclotomicFieldElem with Domain=Eval
	out := NewFieldElemBig(m, prec)
	out.Domain = Eval
	for i := 0; i < m; i++ {
		out.Coeffs[i] = evals[i]
	}
	return out
}
