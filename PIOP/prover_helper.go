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
	"crypto/sha256"
	"encoding/binary"
	"errors"
	"fmt"
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
	prof "vSIS-Signature/prof"
)

var (
	bigA, bigB, bigQ, bigAB big.Int
)

const DEBUG_SUMS = false

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

// modMul returns (a·b) mod q using big.Int temporaries without allocations.
func modMul(a, b, q uint64) uint64 {
	bigA.SetUint64(a)
	bigB.SetUint64(b)
	bigQ.SetUint64(q)
	bigAB.Mul(&bigA, &bigB)
	bigAB.Mod(&bigAB, &bigQ)
	return bigAB.Uint64()
}

// modInv returns a^{‑1} mod q (q must be prime).
func modInv(a, q uint64) uint64 {
	return ring.ModExp(a, q-2, q) // Fermat since q is prime in all params used.
}

// -----------------------------------------------------------------------------
// Fiat–Shamir: tiny deterministic PRF stream (SHA‑256(counter || seed))
// -----------------------------------------------------------------------------

type fsRNG struct {
	seed [32]byte
	ctr  uint64
}

// newFSRNG derives a PRF seed from a label and arbitrary transcript material.
func newFSRNG(label string, material ...[]byte) *fsRNG {
	h := sha256.New()
	h.Write([]byte(label))
	for _, m := range material {
		h.Write(m)
	}
	var s [32]byte
	copy(s[:], h.Sum(nil))
	return &fsRNG{seed: s}
}

func (r *fsRNG) nextU64() uint64 {
	var in [40]byte
	copy(in[:32], r.seed[:])
	binary.LittleEndian.PutUint64(in[32:], r.ctr)
	sum := sha256.Sum256(in[:])
	r.ctr++
	return binary.LittleEndian.Uint64(sum[:])
}

// Helpers to serialize inputs for FS binding.
func bytesU64Vec(v []uint64) []byte {
	out := make([]byte, 8*len(v))
	for i, x := range v {
		binary.LittleEndian.PutUint64(out[8*i:], x)
	}
	return out
}

func bytesU64Mat(M [][]uint64) []byte {
	var out []byte
	for i := range M {
		out = append(out, bytesU64Vec(M[i])...)
	}
	return out
}

func bytesPolys(pp []*ring.Poly) []byte {
	var out []byte
	for _, p := range pp {
		for _, c := range p.Coeffs[0] {
			var b [8]byte
			binary.LittleEndian.PutUint64(b[:], c)
			out = append(out, b[:]...)
		}
	}
	return out
}

// sampleFSPolys(rows × cols) of degree < s, lifted to NTT.
func sampleFSPolys(r *ring.Ring, rows, cols, s int, rng *fsRNG) [][]*ring.Poly {
	out := make([][]*ring.Poly, rows)
	q := r.Modulus[0]
	for i := 0; i < rows; i++ {
		out[i] = make([]*ring.Poly, cols)
		for j := 0; j < cols; j++ {
			p := r.NewPoly()
			for k := 0; k < s; k++ {
				p.Coeffs[0][k] = rng.nextU64() % q
			}
			r.NTT(p, p)
			out[i][j] = p
		}
	}
	return out
}

// sampleFSMatrix(rows × cols) with entries in F_q.
func sampleFSMatrix(rows, cols int, q uint64, rng *fsRNG) [][]uint64 {
	M := make([][]uint64, rows)
	for i := 0; i < rows; i++ {
		M[i] = make([]uint64, cols)
		for j := 0; j < cols; j++ {
			M[i][j] = rng.nextU64() % q
		}
	}
	return M
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

// resetPoly sets all coefficients of p to zero.
func resetPoly(p *ring.Poly) {
	v := p.Coeffs[0]
	for i := range v {
		v[i] = 0
	}
}

// isZeroPoly returns true if all coefficients of p are zero.
func isZeroPoly(p *ring.Poly) bool {
	for _, v := range p.Coeffs[0] {
		if v != 0 {
			return false
		}
	}
	return true
}

// sumEvals returns Σ_{ω∈Ω} P(ω) mod q. scratch must be a *ring.Poly reused by caller.
func sumEvals(r *ring.Ring, P *ring.Poly, omega []uint64, scratch *ring.Poly) uint64 {
	q := r.Modulus[0]
	r.InvNTT(P, scratch)
	coeffs := scratch.Coeffs[0]
	sum := uint64(0)
	for _, w := range omega {
		sum = modAdd(sum, EvalPoly(coeffs, w%q, q), q)
	}
	return sum
}

// sumPolyList computes ΣΩ for each polynomial in list.
func sumPolyList(r *ring.Ring, polys []*ring.Poly, omega []uint64) []uint64 {
	out := make([]uint64, len(polys))
	scratch := r.NewPoly()
	for i, p := range polys {
		out[i] = sumEvals(r, p, omega, scratch)
	}
	return out
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
BuildMaskPolynomials returns ρ random polynomials M₁…Mρ of degree ≤ dQ in
NTT form whose constant term cancels the Ω‑sum of the final
Qᵢ(X) = Mᵢ(X) + Σ_t Γ′_{i,t}F_par,t(X) + Σ_u γ′_{i,u}F_agg,u(X).

Args:
  ringQ       – context (modulus q must not divide len(Ω))
  dQ          – max degree (typically s+ℓ−1)
  omega       – evaluation set Ω
  Fpar        – slice of parallel constraint polynomials
  Fagg        – slice of aggregated constraint polynomials
  GammaPrime  – ρ×|F_par| polynomials Γ′
  gammaPrime  – ρ×|F_agg| scalars γ′

For each i, random coefficients a₁…a_{dQ} are chosen, then a₀ is set so that
ΣΩ Qᵢ(ω) = 0. It panics if q divides |Ω| or Ω has duplicates.
*/
func BuildMaskPolynomials(ringQ *ring.Ring, dQ int, omega []uint64, Fpar []*ring.Poly, Fagg []*ring.Poly, GammaPrime [][]*ring.Poly, gammaPrime [][]uint64) []*ring.Poly {
       defer prof.Track(time.Now(), "BuildMaskPolynomials")

       q := ringQ.Modulus[0]
       s := uint64(len(omega))

       if s == 0 {
               panic("Ω must be non-empty")
       }
       seen := make(map[uint64]struct{}, len(omega))
       for _, w := range omega {
               wm := w % q
               if _, ok := seen[wm]; ok {
                       panic(fmt.Sprintf("BuildMaskPolynomials: Ω contains duplicate element %d (mod q)", wm))
               }
               seen[wm] = struct{}{}
       }
       if q%s == 0 {
               panic(fmt.Sprintf("BuildMaskPolynomials: q (= %d) is multiple of |Ω| (= %d) – constraint not solvable", q, s))
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

       // Precompute evaluations of Fpar at Ω
       fpEval := make([][]uint64, len(Fpar))
       tmpPoly := ringQ.NewPoly()
       for t := range Fpar {
               ringQ.InvNTT(Fpar[t], tmpPoly)
               fpEval[t] = make([]uint64, len(omega))
               for j, w := range omega {
                       fpEval[t][j] = EvalPoly(tmpPoly.Coeffs[0], w%q, q)
               }
       }

       // Precompute ΣΩ Fagg_u
       sumFagg := make([]uint64, len(Fagg))
       for u := range Fagg {
               ringQ.InvNTT(Fagg[u], tmpPoly)
               var sum uint64
               for _, w := range omega {
                       sum = modAdd(sum, EvalPoly(tmpPoly.Coeffs[0], w%q, q), q)
               }
               sumFagg[u] = sum
       }

       // -- modular inverse of S₀ -------------------------------------------------
       // q is prime in lattice settings, so inverse exists.
       invS0 := ring.ModExp(S[0], q-2, q)

       // -- allocate output -------------------------------------------------------
       rho := len(GammaPrime)
       M := make([]*ring.Poly, rho)

       gpCoeff := ringQ.NewPoly()
       for i := 0; i < rho; i++ {
               coeffs := make([]uint64, ringQ.N)
               for k := 1; k <= dQ; k++ {
                       coeffs[k] = randUint64Mod(q)
               }
               var tmp uint64
               for k := 1; k <= dQ; k++ {
                       tmp = modAdd(tmp, modMul(coeffs[k], S[k], q), q)
               }
               for t := range Fpar {
                       ringQ.InvNTT(GammaPrime[i][t], gpCoeff)
                       for j, w := range omega {
                               g := EvalPoly(gpCoeff.Coeffs[0], w%q, q)
                               tmp = modAdd(tmp, modMul(g, fpEval[t][j], q), q)
                       }
               }
               for u := range Fagg {
                       tmp = modAdd(tmp, modMul(gammaPrime[i][u], sumFagg[u], q), q)
               }
               coeffs[0] = modMul(modSub(0, tmp%q, q), invS0, q)

               p := ringQ.NewPoly()
               copy(p.Coeffs[0], coeffs[:])
               ringQ.NTT(p, p)
               if DEBUG_SUMS {
                       coeff := ringQ.NewPoly()
                       ringQ.InvNTT(p, coeff)
                       sum := uint64(0)
                       for _, w := range omega {
                               sum = modAdd(sum, EvalPoly(coeff.Coeffs[0], w%q, q), q)
                       }
                       fmt.Printf("[mask %d] ΣΩ M_i = %d\n", i, sum)
               }
               M[i] = p
       }
       return M
}
