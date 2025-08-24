package PIOP

import (
	"math/big"
	mrand "math/rand"
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
)

func TestInterpolate(t *testing.T) {
	q := uint64(12289)
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
	q := uint64(12289)
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
	q := uint64(12289)
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
	q := uint64(12289)
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
	q := uint64(12289)
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
	q := uint64(12289)
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

func evalAt(r *ring.Ring, p *ring.Poly, x uint64) uint64 {
	coeff := r.NewPoly()
	r.InvNTT(p, coeff)
	return EvalPoly(coeff.Coeffs[0], x%r.Modulus[0], r.Modulus[0])
}

func TestIntegerL2Gadget(t *testing.T) {
	N := 16
	q := uint64(12289)
	ringQ, _ := ring.NewRing(N, []uint64{q})
	ell := 1
	beta := uint64(20)
	wbits := 3
	spec := NewBoundSpec(q, beta, wbits, ell)

	s := 5
	omega := []uint64{1, 2, 3, 4, 5}
	if err := checkOmega(omega, q); err != nil {
		t.Fatal(err)
	}
	S0 := uint64(len(omega))
	S0inv := modInv(S0, q)

	mSig := 2
	w1 := make([]*ring.Poly, mSig)
	for tIdx := 0; tIdx < mSig; tIdx++ {
		vals := make([]uint64, s)
		for j := 0; j < s; j++ {
			vals[j] = uint64(tIdx + j + 1)
		}
		w1[tIdx] = buildValueRow(ringQ, vals, omega, ell)
	}

	cols := appendDecompositionColumns(ringQ, spec.LS, spec.W)
	Wc := 1
	for (1 << (Wc - 1)) < len(omega) {
		Wc++
	}
	glob := appendGlobalCarrys(ringQ, spec.LS, Wc)
	slack := appendGlobalSlack(ringQ, spec.LS, spec.W)

	Sqs, err := ProverFillIntegerL2(ringQ, w1, mSig, spec, cols, glob, slack, omega, ell, S0, S0inv)
	if err != nil {
		t.Fatal(err)
	}

	// Per-column checks
	for j, w := range omega {
		t0 := evalAt(ringQ, cols.T[0], w)
		sVal := evalAt(ringQ, Sqs, w)
		if t0 != sVal {
			t.Fatalf("T0!=Sqs at col %d", j)
		}
		prev := t0
		for l := 0; l < spec.LS; l++ {
			dval := evalAt(ringQ, cols.D[l], w)
			next := evalAt(ringQ, cols.T[l+1], w)
			lhs := modSub(prev, dval, q)
			rhs := modMul(next, spec.R%q, q)
			if lhs != rhs {
				t.Fatalf("remainder chain fail at col %d limb %d", j, l)
			}
			// digits from bits
			var fromBits uint64
			for u := 0; u < spec.W; u++ {
				bit := evalAt(ringQ, cols.Bit[l][u], w)
				if modSub(modMul(bit, bit, q), bit, q) != 0 {
					t.Fatalf("bitness fail at l %d u %d", l, u)
				}
				fromBits = modAdd(fromBits, (bit<<uint(u))%q, q)
			}
			if fromBits != dval {
				t.Fatalf("digit mismatch at col %d limb %d", j, l)
			}
			prev = next
		}
		if prev != 0 {
			t.Fatalf("final remainder nonzero at col %d", j)
		}
	}

	// Aggregated limb identity
	for l := 0; l < spec.LS; l++ {
		sumD := uint64(0)
		for _, w := range omega {
			sumD = modAdd(sumD, evalAt(ringQ, cols.D[l], w), q)
		}
		cVal := evalAt(ringQ, glob.C[l], 0)
		dVal := evalAt(ringQ, slack.D[l], 0)
		bScaled := modMul(spec.Beta2[l]%q, S0inv%q, q)
		cNext := evalAt(ringQ, glob.C[l+1], 0)
		lhs := modAdd(modAdd(modMul(sumD, S0inv, q), cVal, q), dVal, q)
		lhs = modSub(lhs, bScaled, q)
		rhs := modMul(cNext, spec.R%q, q)
		if lhs != rhs {
			t.Fatalf("aggregated limb fail at %d", l)
		}
	}

	// Telescoping equality
	sumSqs := uint64(0)
	for _, w := range omega {
		sumSqs = modAdd(sumSqs, evalAt(ringQ, Sqs, w), q)
	}
	totalDelta := uint64(0)
	powR := uint64(1)
	for l := 0; l < spec.LS; l++ {
		dVal := evalAt(ringQ, slack.D[l], 0) // Δ_l * S0inv
		totalDelta = modAdd(totalDelta, modMul(dVal, powR, q), q)
		powR = modMul(powR, spec.R%q, q)
	}
	lhs := modAdd(sumSqs, modMul(S0, totalDelta, q), q)
	beta2 := (beta * beta) % q
	if lhs != beta2 {
		t.Fatalf("telescoping fail: got %d want %d", lhs, beta2)
	}
}

func TestIntegerL2GadgetTamper(t *testing.T) {
	N := 16
	q := uint64(12289)
	ringQ, _ := ring.NewRing(N, []uint64{q})
	ell := 1
	beta := uint64(20)
	wbits := 3
	spec := NewBoundSpec(q, beta, wbits, ell)
	s := 5
	omega := []uint64{1, 2, 3, 4, 5}
	if err := checkOmega(omega, q); err != nil {
		t.Fatal(err)
	}
	S0 := uint64(len(omega))
	S0inv := modInv(S0, q)
	mSig := 1
	w1 := make([]*ring.Poly, mSig)
	p := ringQ.NewPoly()
	for j := 0; j < s; j++ {
		p.Coeffs[0][j] = uint64(j + 1)
	}
	ringQ.NTT(p, p)
	w1[0] = p
	cols := appendDecompositionColumns(ringQ, spec.LS, spec.W)
	Wc := 1
	for (1 << (Wc - 1)) < len(omega) {
		Wc++
	}
	glob := appendGlobalCarrys(ringQ, spec.LS, Wc)
	slack := appendGlobalSlack(ringQ, spec.LS, spec.W)
	_, err := ProverFillIntegerL2(ringQ, w1, mSig, spec, cols, glob, slack, omega, ell, S0, S0inv)
	if err != nil {
		t.Fatal(err)
	}
	// tamper one bit
	coeff := ringQ.NewPoly()
	ringQ.InvNTT(cols.Bit[0][0], coeff)
	coeff.Coeffs[0][0] ^= 1
	ringQ.NTT(coeff, coeff)
	cols.Bit[0][0] = coeff

	w := omega[0]
	dval := evalAt(ringQ, cols.D[0], w)
	var fromBits uint64
	for u := 0; u < spec.W; u++ {
		bit := evalAt(ringQ, cols.Bit[0][u], w)
		fromBits = modAdd(fromBits, (bit<<uint(u))%q, q)
	}
	if fromBits == dval {
		t.Fatalf("tamper undetected")
	}
}
