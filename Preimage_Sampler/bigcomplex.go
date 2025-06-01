package Preimage_Sampler

import (
	"fmt"
	"math"
	"math/big"
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

// ConvertFromPolyBig embeds p into evaluation domain of K_{2N}.
func ConvertFromPolyBig(r *ring.Ring, p *ring.Poly, prec uint) *CyclotomicFieldElem {
	n := r.N
	vals := make([]complex128, n)
	for i := 0; i < n; i++ {
		vals[i] = complex(ModQToFloat64(p.Coeffs[0][i], r.Modulus[0]), 0)
	}
	fft := FFT(vals, n)
	res := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		re := big.NewFloat(real(fft[i])).SetPrec(prec)
		im := big.NewFloat(imag(fft[i])).SetPrec(prec)
		res.Coeffs[i] = NewBigComplexFromFloat(re, im)
	}
	return res
}

// ConvertToPolyBig interpolates back to ring.Poly in R_q.
func ConvertToPolyBig(f *CyclotomicFieldElem, r *ring.Ring) *ring.Poly {
	n := f.N
	vals := make([]complex128, n)
	for i := 0; i < n; i++ {
		real, _ := f.Coeffs[i].Real.Float64()
		im, _ := f.Coeffs[i].Imag.Float64()
		vals[i] = complex(real, im)
	}
	td := IFFT(vals, n)
	P := r.NewPoly()
	qf := float64(r.Modulus[0])
	for i := 0; i < n && i < r.N; i++ {
		v := real(td[i])
		ri := int64(math.Round(v))
		ri = (ri%int64(qf) + int64(qf)) % int64(qf)
		P.Coeffs[0][i] = uint64(ri)
	}
	return P
}

// SampleGaussianFieldElemBig samples Gaussian noise in evaluation form.
func SampleGaussianFieldElemBig(n int, sigma float64, prec uint) *CyclotomicFieldElem {
	res := NewFieldElemBig(n, prec)
	for i := 0; i < n; i++ {
		re := rand.NormFloat64() * sigma
		im := rand.NormFloat64() * sigma
		res.Coeffs[i] = NewBigComplexFromFloat(
			big.NewFloat(re).SetPrec(prec),
			big.NewFloat(im).SetPrec(prec),
		)
	}
	return res
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

// ComplexInterpolate back to Poly via IFFT and NTT.
func ComplexInterpolate(f *CyclotomicFieldElem, ringQ *ring.Ring) *ring.Poly {
	n := f.N
	vals := make([]complex128, n)
	for i := 0; i < n; i++ {
		real, _ := f.Coeffs[i].Real.Float64()
		imag, _ := f.Coeffs[i].Imag.Float64()
		vals[i] = complex(real, imag)
	}
	coeffs := IFFT(vals, n)
	P := ringQ.NewPoly()
	for j := 0; j < n; j++ {
		ri := int64(math.Round(real(coeffs[j])))
		ri = (ri%int64(ringQ.Modulus[0]) + int64(ringQ.Modulus[0])) % int64(ringQ.Modulus[0])
		P.Coeffs[0][j] = uint64(ri)
	}
	ringQ.NTT(P, P)
	return P
}

// ComplexEvaluateSub embeds with twist and IFFT.
func ComplexEvaluateSub(p *ring.Poly, m int, ringQ *ring.Ring, prec uint) *CyclotomicFieldElem {
	mod := ringQ.Modulus[0]
	evals := make([]complex128, m)
	for j := 0; j < m; j++ {
		s := UnsignedToSigned(p.Coeffs[0][j], mod)
		theta := math.Pi * float64(j) / float64(m)
		evals[j] = complex(float64(s)*math.Cos(theta), float64(s)*math.Sin(theta))
	}
	fft := IFFT(evals, m)
	out := NewFieldElemBig(m, prec)
	for k := 0; k < m; k++ {
		out.Coeffs[k] = NewBigComplex(real(fft[k]), imag(fft[k]), prec)
	}
	out.Domain = Eval
	return out
}

// ComplexInterpolateSub inverts embedding with FFT and twist.
func ComplexInterpolateSub(f *CyclotomicFieldElem, m int, ringQ *ring.Ring) *ring.Poly {
	mod := ringQ.Modulus[0]
	evals := make([]complex128, m)
	for k := 0; k < m; k++ {
		r, _ := f.Coeffs[k].Real.Float64()
		i, _ := f.Coeffs[k].Imag.Float64()
		evals[k] = complex(r, i)
	}
	b := FFT(evals, m)
	P := ringQ.NewPoly()
	for j := 0; j < m; j++ {
		theta := -math.Pi * float64(j) / float64(m)
		aj := b[j] * complex(math.Cos(theta), math.Sin(theta))
		si := int64(math.Round(real(aj)))
		si = (si%int64(mod) + int64(mod)) % int64(mod)
		P.Coeffs[0][j] = uint64(si)
	}
	return P
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
