package Preimage_Sampler

import (
	"math/rand"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
)

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

// referenceTranspose implements the polynomial transpose used in Field_Elements_test.go
func referenceTranspose(r *ring.Ring, p *ring.Poly) *ring.Poly {
	N := r.N
	out := r.NewPoly()
	for lvl, q := range r.Modulus {
		out.Coeffs[lvl][0] = p.Coeffs[lvl][0] % q
		for i := 1; i < N; i++ {
			coeff := p.Coeffs[lvl][N-i] % q
			if coeff != 0 {
				out.Coeffs[lvl][i] = (q - coeff) % q
			} else {
				out.Coeffs[lvl][i] = 0
			}
		}
	}
	return out
}

func TestAutomorphismTransposeRoundTrip(t *testing.T) {
	const (
		n = 16
		q = uint64(97)
	)
	r := makeSmallRing(n, q)

	for trial := 0; trial < 5; trial++ {
		p := r.NewPoly()
		for i := 0; i < n; i++ {
			p.Coeffs[0][i] = uint64(rand.Intn(int(q)))
		}
		t1 := AutomorphismTranspose(r, p)
		t2 := AutomorphismTranspose(r, t1)
		for i := 0; i < n; i++ {
			if p.Coeffs[0][i]%q != t2.Coeffs[0][i]%q {
				t.Fatalf("round-trip mismatch at %d: got %d want %d", i, t2.Coeffs[0][i], p.Coeffs[0][i])
			}
		}
	}
}

func TestAutomorphismTransposeReference(t *testing.T) {
	const (
		n = 16
	)
	moduli := []uint64{97, 193}
	r, err := ring.NewRing(n, moduli)
	if err != nil {
		t.Fatalf("failed creating ring: %v", err)
	}

	for trial := 0; trial < 5; trial++ {
		p := r.NewPoly()
		for lvl, q := range moduli {
			for i := 0; i < n; i++ {
				p.Coeffs[lvl][i] = uint64(rand.Intn(int(q)))
			}
		}

		want := referenceTranspose(r, p)
		got := AutomorphismTranspose(r, p)

		for lvl, q := range moduli {
			for i := 0; i < n; i++ {
				if want.Coeffs[lvl][i]%q != got.Coeffs[lvl][i]%q {
					t.Fatalf("lvl %d index %d: got %d want %d", lvl, i, got.Coeffs[lvl][i], want.Coeffs[lvl][i])
				}
			}
		}
	}
}
