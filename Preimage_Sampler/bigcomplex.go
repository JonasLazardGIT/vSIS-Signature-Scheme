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

// FFTBig takes a slice of length n (power of two) of *BigComplex (coefficient form),
// and returns a new slice of length n (evaluation form) at the primitive 2n-th roots.
// “prec” is the bit‐precision for every BigComplex operation and twiddle‐factor.
//
// (Implementation shown in the previous answer.)
// --------------------------------------------------------------------------------

// IFFTBig takes a slice of length n (evaluation), and returns a slice of length n
// (coefficient) via the inverse FFT at the same precision.  See previous answer.
// --------------------------------------------------------------------------------

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

// --------------------------------------------------------------------------------
//
//	ComplexInterpolate
//
//	The old “ComplexInterpolate” did:
//
//	  1) Build []complex128 from f.Coeffs[*].Real/Imag via .Float64().
//	  2) Called IFFT(...) on []complex128 to get []complex128 of length n.
//	  3) Built a ring.Poly by rounding each real(...) to nearest integer mod q, and then NTT.
//
//	We now replace step (2) with IFFTBig on []*BigComplex, then round.
//
//	Inputs:
//	  - f:     *CyclotomicFieldElem (Domain should be Eval, but we only read .Coeffs).
//	  - ringQ: *ring.Ring used to create the new ring.Poly and do the final NTT.
//
//	Output:
//	  - a *ring.Poly (NTT‐domain).
//
// --------------------------------------------------------------------------------
func ComplexInterpolate(f *CyclotomicFieldElem, ringQ *ring.Ring) *ring.Poly {
	n := f.N

	// 1) Copy f.Coeffs into a new slice of *BigComplex
	prec := f.Coeffs[0].Real.Prec()
	evalSlice := make([]*BigComplex, n)
	for i := 0; i < n; i++ {
		evalSlice[i] = f.Coeffs[i].Copy()
	}

	// 2) Perform a high‐precision IFFT:
	coeffSlice := IFFTBig(evalSlice, prec)

	// 3) Round each big‐float to nearest integer mod q, place into a new ring.Poly
	P := ringQ.NewPoly()
	q := int64(ringQ.Modulus[0])
	for j := 0; j < n; j++ {
		realBF := coeffSlice[j].Real
		f64, _ := realBF.Float64()
		rInt := int64(math.Round(f64))
		rInt = ((rInt % q) + q) % q
		P.Coeffs[0][j] = uint64(rInt)
	}

	// 4) Finally, go to NTT domain:
	ringQ.NTT(P, P)
	return P
}

// --------------------------------------------------------------------------------
//
//	ComplexEvaluateSub
//
//	Old version did:
//	  1) Build evals[j] = s⋅e^{iθ} at float64 precision, for θ = π j / m.
//	  2) Called IFFT(evals, m) in float128, then put into *BigComplex.
//
//	We now do exactly the same, except we call IFFTBig instead of IFFT.
//
//	Inputs:
//	  - p:     *ring.Poly (coeff‐domain, length m ≤ ringQ.N).
//	  - m:     target transform length (power of two).
//	  - ringQ: *ring.Ring (only used to read modulus).
//	  - prec:  bit‐precision for all BigComplex arithmetic.
//
//	Output:
//	  - a new *CyclotomicFieldElem of length m, in domain Eval.
//
// --------------------------------------------------------------------------------
func ComplexEvaluateSub(p *ring.Poly, m int, ringQ *ring.Ring, prec uint) *CyclotomicFieldElem {
	mod := ringQ.Modulus[0]

	// 1) Build “twisted” BigComplex input of length m:
	bigSlice := make([]*BigComplex, m)
	for j := 0; j < m; j++ {
		sInt := UnsignedToSigned(p.Coeffs[0][j], mod)
		θ := math.Pi * float64(j) / float64(m)

		// ** MUST use e^(−i·θ )  =  cos(θ) − i·sin(θ) **
		re := float64(sInt) * math.Cos(θ)
		im := float64(sInt) * -math.Sin(θ) // ← negative sign here

		reBF := new(big.Float).SetPrec(prec).SetFloat64(re)
		imBF := new(big.Float).SetPrec(prec).SetFloat64(im)
		bigSlice[j] = &BigComplex{Real: reBF, Imag: imBF}
	}

	// 2) Call forward FFT (not IFFT):
	fftResult := FFTBig(bigSlice, prec)

	// 3) Copy result into a CyclotomicFieldElem, mark Domain=Eval
	out := NewFieldElemBig(m, prec)
	out.Domain = Eval
	for k := 0; k < m; k++ {
		out.Coeffs[k] = fftResult[k].Copy()
	}
	return out
}

// --------------------------------------------------------------------------------
//
//	ComplexInterpolateSub
//
//	Old version did:
//	  1) Build a []complex128 “evals” from f.Coeffs[*].Real/Imag.
//	  2) Called FFT(evals, m) in float128, then applied a “twist” cos/sin to each FFT[k].
//	  3) Rounded to ring.Poly.
//
//	We now simply call FFTBig(...) on []*BigComplex, then multiply each output by
//	the twist e^{i·(−π j / m)} (again as a BigComplex), then round.
//
//	Inputs:
//	  - f:     *CyclotomicFieldElem (length m, Domain=Eval).
//	  - m:     transform length.
//	  - ringQ: *ring.Ring used to know q.
//	Output:
//	  - a new *ring.Poly in coefficient form (no final NTT).
//
// --------------------------------------------------------------------------------
func ComplexInterpolateSub(f *CyclotomicFieldElem, m int, ringQ *ring.Ring) *ring.Poly {
	mod := int64(ringQ.Modulus[0])
	prec := f.Coeffs[0].Real.Prec()

	// 1) Copy f.Coeffs into a new slice of *BigComplex
	evalSlice := make([]*BigComplex, m)
	for k := 0; k < m; k++ {
		evalSlice[k] = f.Coeffs[k].Copy()
	}

	// 2) Compute high‐precision FFT:
	fftResult := FFTBig(evalSlice, prec)

	// 3) Apply the “twist” factor e^{−i π j / m} to each fftResult[j]:
	for j := 0; j < m; j++ {
		θ := -math.Pi * float64(j) / float64(m) // negative for this “twist”
		// Build a BigComplex representing e^{iθ} = cos(θ) + i sin(θ):
		cosF := big.NewFloat(0).SetPrec(prec).SetFloat64(math.Cos(θ))
		sinF := big.NewFloat(0).SetPrec(prec).SetFloat64(math.Sin(θ))
		twist := &BigComplex{Real: cosF, Imag: sinF}

		// Multiply fftResult[j] by “twist” in big‐float:
		fftResult[j] = fftResult[j].Mul(twist)
	}

	// 4) Now round each coordinate to nearest integer mod q:
	P := ringQ.NewPoly()
	for j := 0; j < m; j++ {
		realBF := fftResult[j].Real
		f64, _ := realBF.Float64()
		si := int64(math.Round(f64))
		si = ((si % mod) + mod) % mod
		P.Coeffs[0][j] = uint64(si)
	}

	return P
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

// HermitianTransposeFieldElem applies f^t(x).
func HermitianTransposeFieldElem(f *CyclotomicFieldElem) *CyclotomicFieldElem {
	n := f.N
	res := NewFieldElemBig(n, f.Coeffs[0].Real.Prec())
	for i := 0; i < n; i++ {
		rev := (n - i) % n
		res.Coeffs[i] = f.Coeffs[rev].Conj()
	}
	return res
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

// ComplexEvaluate maps p in NTT domain to evaluation in K_{2n}.
func ComplexEvaluate(p *ring.Poly, ringQ *ring.Ring, prec uint) *CyclotomicFieldElem {
	n := ringQ.N
	vals := make([]complex128, n)
	for j := 0; j < n; j++ {
		s := UnsignedToSigned(p.Coeffs[0][j], ringQ.Modulus[0])
		vals[j] = complex(float64(s), 0)
	}
	fft := FFT(vals, n)
	out := NewFieldElemBig(n, prec)
	for k := 0; k < n; k++ {
		out.Coeffs[k] = NewBigComplex(real(fft[k]), imag(fft[k]), prec)
	}
	out.Domain = Eval
	return out
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

// SwitchToCoeff ensures the element is in the coefficient domain.
func (e *CyclotomicFieldElem) SwitchToCoeff(ringQ *ring.Ring) {
	if e.Domain == Coeff {
		return
	}
	// 1) Interpolate from evaluation domain back to a ring.Poly
	poly := ConvertToPolyBig(e, ringQ)

	// 2) Overwrite e.Coeffs with the polynomial coefficients (as BigComplex)
	prec := e.Coeffs[0].Real.Prec()
	for i := 0; i < e.N && i < ringQ.N; i++ {
		// Convert uint64 coefficient to BigComplex
		realVal := float64(poly.Coeffs[0][i])
		re := new(big.Float).SetPrec(prec).SetFloat64(realVal)
		im := new(big.Float).SetPrec(prec).SetFloat64(0)
		e.Coeffs[i] = NewBigComplexFromFloat(re, im)
	}

	// 3) Mark domain
	e.Domain = Coeff
}

// SwitchToEval ensures the element is in the evaluation (NTT) domain.
func (e *CyclotomicFieldElem) SwitchToEval(ringQ *ring.Ring) {
	if e.Domain == Eval {
		return
	}
	// 1) Build a ring.Poly from the coefficient-domain BigComplex values
	P := ringQ.NewPoly()
	for i := 0; i < e.N && i < ringQ.N; i++ {
		realF, _ := e.Coeffs[i].Real.Float64()
		// Round and reduce mod q
		q0 := float64(ringQ.Modulus[0])
		ri := int64(math.Round(realF)) % int64(q0)
		if ri < 0 {
			ri += int64(q0)
		}
		P.Coeffs[0][i] = uint64(ri)
	}

	// 2) Convert polynomial to evaluation-domain CyclotomicFieldElem
	evalElem := ConvertFromPolyBig(ringQ, P, e.Coeffs[0].Real.Prec())

	// 3) Overwrite e.Coeffs with the FFT output
	e.Coeffs = evalElem.Coeffs

	// 4) Mark domain
	e.Domain = Eval
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
