// Package PIOP – prover‑side helpers for Section 6 PACS construction.
//
// This file contains utilities that turn a *row* of the witness matrix
// (w_{i,1},…,w_{i,s}) into the blinding polynomial P_i ∈ F_q[X] described in
// Protocol 6 of the SmallWood paper.  Each polynomial is defined by
//   - its s anchor points Ω = {ω₁,…,ω_s},
//   - ℓ extra random points r₁,…,r_ℓ (r_j ∉ Ω),
//   - uniformly random evaluations P_i(r_j) ∈ F_q.
//
// Degree ≤ s+ℓ−1 and perfect zero‑knowledge follow immediately.
//
// The code is written for tiny s,ℓ (≤64) as is common in practice;  O(n³)
// barycentric Lagrange interpolation is therefore acceptable.
package PIOP

import (
	"crypto/rand"
	"encoding/binary"
	"errors"
	"log"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
	prof "vSIS-Signature/prof"
)

// -----------------------------------------------------------------------------
//  field helpers (mod q)
// -----------------------------------------------------------------------------

func randUint64Mod(q uint64) uint64 {
	var bound uint64 = ^uint64(0) - (^uint64(0) % q)
	for {
		var buf [8]byte
		if _, err := rand.Read(buf[:]); err != nil {
			panic("randUint64Mod: entropy read failed: " + err.Error())
		}
		v := binary.LittleEndian.Uint64(buf[:])
		if v < bound {
			return v % q
		}
	}
}

// modAdd returns (a+b) mod q.
func modAdd(a, b, q uint64) uint64 {
	s := a + b
	if s >= q || s < a { // handle wrap‑around
		s -= q
	}
	return s
}

// modSub returns (a‑b) mod q.
func modSub(a, b, q uint64) uint64 {
	if a >= b {
		return a - b
	}
	return a + q - b
}

// modMul returns (a·b) mod q –  constant‑time big‑uint multiplication.
func modMul(a, b, q uint64) uint64 {
	return (uint64(new(big.Int).Mul(big.NewInt(int64(a)), big.NewInt(int64(b))).Mod(
		new(big.Int).Mul(big.NewInt(int64(a)), big.NewInt(int64(b))), big.NewInt(int64(q))).Uint64()))
}

// modInv returns a^{‑1} mod q (q must be prime).
func modInv(a, q uint64) uint64 {
	return ring.ModExp(a, q-2, q) // Fermat since q is prime in all params used.
}

// randFieldElem draws a uniform element in [0,q) not in the forbidden set.
func randFieldElem(q uint64, forbid map[uint64]struct{}) (uint64, error) {
	if q == 0 {
		return 0, errors.New("q=0")
	}
	bound := ^uint64(0) - (^uint64(0) % q)
	for {
		var buf [8]byte
		if _, err := rand.Read(buf[:]); err != nil {
			return 0, err
		}
		v := uint64(buf[0]) | uint64(buf[1])<<8 | uint64(buf[2])<<16 | uint64(buf[3])<<24 | uint64(buf[4])<<32 | uint64(buf[5])<<40 | uint64(buf[6])<<48 | uint64(buf[7])<<56
		if v >= bound {
			continue
		}
		v %= q
		if _, bad := forbid[v]; !bad {
			return v, nil
		}
	}
}

// -----------------------------------------------------------------------------
//  polynomial helpers
// -----------------------------------------------------------------------------

// polyMul naive O(n²)  – sufficient for n ≤ 64.
func polyMul(a, b []uint64, q uint64) []uint64 {
	out := make([]uint64, len(a)+len(b)-1)
	for i, av := range a {
		for j, bv := range b {
			out[i+j] = modAdd(out[i+j], modMul(av, bv, q), q)
		}
	}
	return out
}

// scalePoly returns c·p  (mod q).
func scalePoly(p []uint64, c, q uint64) []uint64 {
	out := make([]uint64, len(p))
	for i, v := range p {
		out[i] = modMul(v, c, q)
	}
	return out
}

// addInto in‑place: dst += src mod q (resize dst if needed).
func addIntoUint(dst *[]uint64, src []uint64, q uint64) {
	if len(src) > len(*dst) {
		newDst := make([]uint64, len(src))
		copy(newDst, *dst)
		*dst = newDst
	}
	for i, v := range src {
		(*dst)[i] = modAdd((*dst)[i], v, q)
	}
}

// lagrangeBasisNumerator returns Π_{j≠i} (X - x_j) as a coefficient slice.
func lagrangeBasisNumerator(xs []uint64, i int, q uint64) []uint64 {
	num := []uint64{1}
	for j, xj := range xs {
		if j == i {
			continue
		}
		num = polyMul(num, []uint64{modSub(0, xj, q), 1}, q) // (X - xj)
	}
	return num
}

// Interpolate returns the coefficients of the unique poly of degree <len(xs)
// that satisfies P(xs[k]) = ys[k].  xs must be distinct.
func Interpolate(xs, ys []uint64, q uint64) []uint64 {
	n := len(xs)
	res := make([]uint64, 1) // zero‑poly
	for i := 0; i < n; i++ {
		num := lagrangeBasisNumerator(xs, i, q)
		// denom = Π_{j≠i} (xs[i]-xs[j])
		denom := uint64(1)
		for j, xj := range xs {
			if j == i {
				continue
			}
			denom = modMul(denom, modSub(xs[i], xj, q), q)
		}
		coeff := modMul(ys[i], modInv(denom, q), q)
		term := scalePoly(num, coeff, q)
		addIntoUint(&res, term, q)
	}
	// trim trailing zeros
	for len(res) > 1 && res[len(res)-1] == 0 {
		res = res[:len(res)-1]
	}
	return res
}

// -----------------------------------------------------------------------------
//  Public helper – build P_i(X) with random blinding
// -----------------------------------------------------------------------------

// BuildRowPolynomial takes a witness row (length s), the corresponding
// domain Ω, and a blinding parameter ℓ.  It returns
//   - *ring.Poly in NTT form (degree ≤ s+ℓ-1),
//   - the extra points r[0:ℓ],
//   - their evaluations y[0:ℓ].
//
// Pre‑conditions:  len(row)==len(omega)==s,   ℓ≥1,   xs are all distinct.
func BuildRowPolynomial(ringQ *ring.Ring, row, omega []uint64, ell int) (poly *ring.Poly, rPoints, rEvals []uint64, err error) {
	defer prof.Track(time.Now(), "BuildRowPolynomial")
	if len(row) != len(omega) {
		return nil, nil, nil, errors.New("row and omega length mismatch")
	}
	if ell <= 0 {
		return nil, nil, nil, errors.New("ell must be ≥1")
	}
	q := ringQ.Modulus[0]

	// 1. choose ℓ random points outside Ω
	forbid := make(map[uint64]struct{}, len(omega))
	for _, w := range omega {
		forbid[w] = struct{}{}
	}
	rPoints = make([]uint64, ell)
	for i := 0; i < ell; i++ {
		rp, e := randFieldElem(q, forbid)
		if e != nil {
			return nil, nil, nil, e
		}
		forbid[rp] = struct{}{}
		rPoints[i] = rp
	}

	// 2. choose ℓ random evaluations y_i
	rEvals = make([]uint64, ell)
	for i := 0; i < ell; i++ {
		y, e := randFieldElem(q, nil)
		if e != nil {
			return nil, nil, nil, e
		}
		rEvals[i] = y
	}

	// 3. interpolate over xs = Ω ∪ rPoints, ys = row ∪ rEvals
	xs := append(append([]uint64{}, omega...), rPoints...)
	ys := append(append([]uint64{}, row...), rEvals...)
	coeffs := Interpolate(xs, ys, q) // coeff domain

	// 4. lift to NTT and wrap in *ring.Poly
	poly = ringQ.NewPoly()
	copy(poly.Coeffs[0], coeffs)
	ringQ.NTT(poly, poly)
	return poly, rPoints, rEvals, nil
}

// -----------------------------------------------------------------------------
//  BuildMaskPolynomials
// -----------------------------------------------------------------------------
/*
BuildMaskPolynomials returns ρ random polynomials M₁…Mρ of degree ≤ dQ
in NTT form such that      Σ_{ω∈Ω} M_i(ω) = 0     for every i.

Args:
  ringQ   – context (modulus q must not divide len(Ω))
  rho     – number of polynomials
  dQ      – max degree   (from Eq.(3) in the paper, typically dQ = s+ℓ−1 or 2d)
  omega   – evaluation set Ω, given as raw field elements (coeff-domain)

The construction picks random coefficients a₁…a_{dQ} and solves for a₀
so the linear constraint is met:

    a₀ = - ( Σ_{k=1}^{dQ} a_k S_k ) / S₀    with  S_k = Σ_{ω∈Ω} ω^k .

It panics if q | |Ω|.
*/
func BuildMaskPolynomials(ringQ *ring.Ring, rho, dQ int, omega []uint64) []*ring.Poly {
	defer prof.Track(time.Now(), "BuildMaskPolynomials")

	q := ringQ.Modulus[0]
	s := uint64(len(omega))

	if s == 0 {
		log.Fatal("Ω must be non-empty")
	}
	if q%s == 0 {
		log.Fatalf("BuildMaskPolynomials: q (= %d) is multiple of |Ω| (= %d) – constraint not solvable", q, s)
	}

	// -- pre-compute S_k = Σ ω^k  for  k=0…dQ modulo q -----------------------
	S := make([]uint64, dQ+1) // S[0]=|Ω|
	S[0] = s % q

	// powersTmp[j] will hold ω_j^k
	powersTmp := make([]uint64, len(omega))
	for k := 1; k <= dQ; k++ {
		sum := uint64(0)
		for j, w := range omega {
			if k == 1 {
				powersTmp[j] = w % q
			} else {
				powersTmp[j] = (powersTmp[j] * w) % q // ω^k
			}
			sum += powersTmp[j]
			if sum >= q {
				sum -= q
			}
		}
		S[k] = sum
	}

	// -- modular inverse of S₀ -------------------------------------------------
	// q is prime in lattice settings, so inverse exists.
	invS0 := ring.ModExp(S[0], q-2, q)

	// -- allocate output -------------------------------------------------------
	M := make([]*ring.Poly, rho)

	// -- build each polynomial -------------------------------------------------
	for i := 0; i < rho; i++ {
		coeffs := make([]uint64, ringQ.N)
		// 1) random a_k  for k=1…dQ
		for k := 1; k <= dQ; k++ {
			coeffs[k] = randUint64Mod(q)
		}
		// 2) solve for a₀
		var tmp uint64
		for k := 1; k <= dQ; k++ {
			tmp = (tmp + coeffs[k]*S[k]) % q
		}
		coeffs[0] = (q - tmp) % q           // -Σ a_k S_k
		coeffs[0] = (coeffs[0] * invS0) % q // * (S₀)^{-1}

		// 3) build, lift to NTT
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs[:])
		ringQ.NTT(p, p)
		M[i] = p
	}
	return M
}
