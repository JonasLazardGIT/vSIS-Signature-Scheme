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
        ringQ.NTT(Fpar[0], Fpar[0])
        Fagg := []*ring.Poly{}
        Gamma := sampleFSPolys(ringQ, rho, len(Fpar), len(omega), newFSRNG("g"))
        gamma := sampleFSMatrix(rho, len(Fagg), q, newFSRNG("h"))
        M := BuildMaskPolynomials(ringQ, dQ, omega, Fpar, Fagg, Gamma, gamma)
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
        Fpar := []*ring.Poly{}
        Fagg := []*ring.Poly{}
        Gamma := [][]*ring.Poly{{}}
        gamma := [][]uint64{{}}
        BuildMaskPolynomials(ringQ, 3, omega, Fpar, Fagg, Gamma, gamma)
}

func TestOmegaHygieneMultiple(t *testing.T) {
        N := 16
        q := uint64(13)
        ringQ, _ := ring.NewRing(N, []uint64{q})
        omega := make([]uint64, 13)
        for i := range omega {
                omega[i] = uint64(i)
        }
        defer func() {
                if r := recover(); r == nil {
                        t.Fatalf("expected panic when q divisible by |Ω|")
                }
        }()
        Fpar := []*ring.Poly{}
        Fagg := []*ring.Poly{}
        Gamma := [][]*ring.Poly{{}}
        gamma := [][]uint64{{}}
        BuildMaskPolynomials(ringQ, 3, omega, Fpar, Fagg, Gamma, gamma)
}

func TestConstPolyNTT(t *testing.T) {
        N := 16
        q := uint64(97)
        ringQ, _ := ring.NewRing(N, []uint64{q})
        tval := uint64(42)
        p := constPolyNTT(ringQ, tval)
        coeff := ringQ.NewPoly()
        ringQ.InvNTT(p, coeff)
        if coeff.Coeffs[0][0]%q != tval%q {
                t.Fatalf("a0 mismatch: got %d want %d", coeff.Coeffs[0][0]%q, tval%q)
        }
        for i := 1; i < ringQ.N; i++ {
                if coeff.Coeffs[0][i]%q != 0 {
                        t.Fatalf("non-zero coeff at %d", i)
                }
        }
}

func TestColumnsToRowsEval(t *testing.T) {
        N := 16
        q := uint64(97)
        ringQ, _ := ring.NewRing(N, []uint64{q})
        w1 := []*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}
        w2 := ringQ.NewPoly()
        w3 := []*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()}
        for _, p := range append(append([]*ring.Poly{}, w1...), append([]*ring.Poly{w2}, w3...)...) {
                ringQ.NTT(p, p)
        }
        ell := 1
        px := ringQ.NewPoly(); px.Coeffs[0][1] = 1
        pts := ringQ.NewPoly(); ringQ.NTT(px, pts)
        omega := pts.Coeffs[0][:ringQ.N-ell]
        rows := columnsToRows(ringQ, w1, w2, w3, ell, omega)
        qmod := ringQ.Modulus[0]
        coeff := ringQ.NewPoly()
        for i := 0; i < len(w1); i++ {
                ringQ.InvNTT(w1[i], coeff)
                for j, w := range omega {
                        if rows[i][j] != EvalPoly(coeff.Coeffs[0], w%qmod, qmod) {
                                t.Fatalf("row %d col %d mismatch", i, j)
                        }
                }
        }
}

func TestVerifyQ(t *testing.T) {
        N := 16
        q := uint64(97)
        ringQ, _ := ring.NewRing(N, []uint64{q})
        omega := []uint64{2, 4, 6}
        Fpar := []*ring.Poly{ringQ.NewPoly()}
        ringQ.NTT(Fpar[0], Fpar[0])
        Fagg := []*ring.Poly{}
        Gamma := sampleFSPolys(ringQ, 1, len(Fpar), len(omega), newFSRNG("g"))
        gamma := sampleFSMatrix(1, len(Fagg), q, newFSRNG("h"))
        M := BuildMaskPolynomials(ringQ, 3, omega, Fpar, Fagg, Gamma, gamma)
        Q := BuildQ(ringQ, M, Fpar, Fagg, Gamma, gamma)
        if !VerifyQ(ringQ, Q, omega) {
                t.Fatalf("VerifyQ failed on clean input")
        }
        coeff := ringQ.NewPoly()
        ringQ.InvNTT(Q[0], coeff)
        coeff.Coeffs[0][0] ^= 1
        ringQ.NTT(coeff, coeff)
        Q[0] = coeff
        if VerifyQ(ringQ, Q, omega) {
                t.Fatalf("VerifyQ should fail after tamper")
        }
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
