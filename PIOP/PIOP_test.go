package PIOP

import (
	"math/big"
	mrand "math/rand"
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

func TestMaskCancellation(t *testing.T) {
	N := 16
	q := uint64(97)
	ringQ, _ := ring.NewRing(N, []uint64{q})
	omega := []uint64{2, 4, 6}
	rho := 2
	dQ := 5

	// simple Fpar/Fagg
	Fpar := []*ring.Poly{ringQ.NewPoly()}
	Fagg := []*ring.Poly{}
	sumFpar := sumPolyList(ringQ, Fpar, omega)
	Gamma := sampleFSMatrix(rho, len(Fpar), q, newFSRNG("g"))
	gamma := sampleFSMatrix(rho, len(Fagg), q, newFSRNG("h"))
	M := BuildMaskPolynomials(ringQ, rho, dQ, omega, Gamma, gamma, sumFpar, []uint64{})
	Q := BuildQ(ringQ, M, Fpar, Fagg, Gamma, gamma)
	if !VerifyQ(ringQ, Q, omega) {
		t.Fatalf("VerifyQ failed")
	}
}

func TestFieldOps(t *testing.T) {
	q := uint64(1<<61 - 1)
	rnd := mrand.New(mrand.NewSource(0))
	for i := 0; i < 100; i++ {
		a := rnd.Uint64() % q
		b := rnd.Uint64() % q
		bigQ := new(big.Int).SetUint64(q)
		bigA := new(big.Int).SetUint64(a)
		bigB := new(big.Int).SetUint64(b)
		if modAdd(a, b, q) != new(big.Int).Mod(new(big.Int).Add(bigA, bigB), bigQ).Uint64() {
			t.Fatal("modAdd mismatch")
		}
		if modSub(a, b, q) != new(big.Int).Mod(new(big.Int).Sub(bigA, bigB), bigQ).Uint64() {
			t.Fatal("modSub mismatch")
		}
		if modMul(a, b, q) != new(big.Int).Mod(new(big.Int).Mul(bigA, bigB), bigQ).Uint64() {
			t.Fatal("modMul mismatch")
		}
	}
}

func TestEvalPolyRandom(t *testing.T) {
	q := uint64(97)
	rnd := mrand.New(mrand.NewSource(1))
	coeffs := make([]uint64, 5)
	for i := range coeffs {
		coeffs[i] = rnd.Uint64() % q
	}
	x := rnd.Uint64() % q
	// Horner result
	got := EvalPoly(coeffs, x, q)
	// naive evaluation
	naive := uint64(0)
	for i := len(coeffs) - 1; i >= 0; i-- {
		naive = modMul(naive, x, q)
		naive = modAdd(naive, coeffs[i], q)
		if i == 0 {
			break
		}
	}
	if got != naive {
		t.Fatalf("EvalPoly mismatch: got %d want %d", got, naive)
	}
}

func TestOmegaHygiene(t *testing.T) {
	N := 16
	q := uint64(97)
	ringQ, _ := ring.NewRing(N, []uint64{q})
	omega := []uint64{1, 1, 2}
	defer func() {
		if r := recover(); r == nil {
			t.Fatalf("expected panic on duplicate Ω")
		}
	}()
	BuildMaskPolynomials(ringQ, 1, 3, omega, [][]uint64{{1}}, [][]uint64{}, []uint64{0}, []uint64{})
}

func TestBuildThetaPrimeSetMultiRow(t *testing.T) {
	N := 16
	q := uint64(97)
	ringQ, _ := ring.NewRing(N, []uint64{q})
	zero := ringQ.NewPoly()
	A := [][]*ring.Poly{{zero, zero}, {zero, zero}}
	b1 := []*ring.Poly{zero, zero}
	B0Const := []*ring.Poly{zero, zero}
	B0Msg := [][]*ring.Poly{{zero, zero}}
	B0Rnd := [][]*ring.Poly{{zero, zero}}
	omega := []uint64{1, 2}
	tp := BuildThetaPrimeSet(ringQ, A, b1, B0Const, B0Msg, B0Rnd, omega)
	if len(tp.ARows) != 2 || len(tp.ARows[0]) != 2 {
		t.Fatalf("unexpected ARows shape")
	}
}
