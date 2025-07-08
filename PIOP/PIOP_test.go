package PIOP

import (
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
)

func TestInterpolate(t *testing.T) {
	q := uint64(97)
	xs := []uint64{3, 5, 7}
	ys := []uint64{10, 20, 30}
	poly := Interpolate(xs, ys, q)
	for i, x := range xs {
		if got := EvalPoly(poly, x, q); got != ys[i] {
			t.Fatalf("interpolation failed: P(%d)=%d want %d", x, got, ys[i])
		}
	}
}

func TestBuildRowPolynomial(t *testing.T) {
	N := 16
	q := uint64(97)
	ringQ, _ := ring.NewRing(N, []uint64{q})

	row := []uint64{11, 22, 33}
	omega := []uint64{2, 4, 6}
	ell := 2

	polyNTT, rPts, rVals, err := BuildRowPolynomial(ringQ, row, omega, ell)
	if err != nil {
		t.Fatal(err)
	}

	// back to coeff domain
	coeff := polyNTT.CopyNew()
	ringQ.InvNTT(coeff, coeff)

	// check omega evaluations
	for i, w := range omega {
		got := EvalPoly(coeff.Coeffs[0], w, q)
		if got != row[i] {
			t.Fatalf("P(ω_%d) mismatch: got %d want %d", i, got, row[i])
		}
	}
	// check random evals
	for i, r := range rPts {
		got := EvalPoly(coeff.Coeffs[0], r, q)
		if got != rVals[i] {
			t.Fatalf("P(r_%d) mismatch: got %d want %d", i, got, rVals[i])
		}
	}
	// degree bound
	maxDeg := len(row) + ell - 1
	// trim trailing zeros
	deg := len(coeff.Coeffs[0]) - 1
	for deg > 0 && coeff.Coeffs[0][deg] == 0 {
		deg--
	}
	if deg > maxDeg {
		t.Fatalf("degree too large: got %d, bound %d", deg, maxDeg)
	}
}

func TestBuildMaskPolynomials(t *testing.T) {
	N := 16
	q := uint64(97)
	ringQ, _ := ring.NewRing(N, []uint64{q})

	omega := []uint64{2, 4, 6}
	rho := 3
	dQ := 5

	masks := BuildMaskPolynomials(ringQ, rho, dQ, omega)

	for i, mNTT := range masks {
		coeff := mNTT.CopyNew()
		ringQ.InvNTT(coeff, coeff)
		// sum over omega
		sum := uint64(0)
		for _, w := range omega {
			sum = (sum + EvalPoly(coeff.Coeffs[0], w, q)) % q
		}
		if sum != 0 {
			t.Fatalf("mask %d constraint failed: Σ M_i(ω) = %d ≠ 0", i, sum)
		}
		// degree bound
		deg := len(coeff.Coeffs[0]) - 1
		for deg > 0 && coeff.Coeffs[0][deg] == 0 {
			deg--
		}
		if deg > dQ {
			t.Fatalf("mask %d degree %d exceeds dQ=%d", i, deg, dQ)
		}
	}
}
