// quadraticgate/quadratic_gate.go
package quadraticgate

import (
	"errors"
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// -----------------------------------------------------------------------------
// Convenience helpers
// -----------------------------------------------------------------------------

// zeroPoly allocates an all-zero polynomial in NTT form.
func zeroPoly(r *ring.Ring) *ring.Poly { p := r.NewPoly(); return p }

// copyPoly returns a deep copy of p.
func copyPoly(r *ring.Ring, p *ring.Poly) *ring.Poly {
	out := r.NewPoly()
	ring.Copy(p, out)
	return out
}

// mulConst  :  q · x   where q = const mod Modulus[0]
func mulConst(r *ring.Ring, p *ring.Poly, c uint64) *ring.Poly {
	out := r.NewPoly()
	r.MulScalar(p, c, out)
	return out
}

// addInto  :  dst += src
func addInto(r *ring.Ring, dst, src *ring.Poly) { r.Add(dst, src, dst) }

// -----------------------------------------------------------------------------
// Public / private containers
// -----------------------------------------------------------------------------

type GatePublic struct {
	R2 [][]*ring.Poly // sparse symmetric (size k×k)
	R1 []*ring.Poly   // length-k vector
	R0 *ring.Poly
}

type GatePrivate struct {
	X []*ring.Poly // witness vector  (s ‖ x1 ‖ u ‖ x0)
}

// -----------------------------------------------------------------------------
//  BuildGate  —— main entry
// -----------------------------------------------------------------------------
//
// ringQ          :  R_q   (single-prime ring used everywhere)
// A              :  [n][m] matrix – each entry an NTT polynomial
// b1             :  length-n vector (polynomials)
// B0Const        :  length-n vector (constant column of B0)
// B0Msg, B0Rnd   :  slices of columns of B0 (message & randomness)
// rho            :  length-n compression vector (native uint64 coeffs) *
//
// s, x1, u, x0   :  secret vectors (polynomials) that should already
//                  satisfy the proof-friendly equation.
// -----------------------------------------------------------------------------

func BuildGate(
	ringQ *ring.Ring,
	A [][]*ring.Poly,
	b1 []*ring.Poly,
	B0Const []*ring.Poly,
	B0Msg [][]*ring.Poly,
	B0Rnd [][]*ring.Poly,
	rho []uint64, // assumed in coefficient domain, centred in [0,q)
	/* private */
	s []*ring.Poly,
	x1 *ring.Poly,
	u []*ring.Poly,
	x0 []*ring.Poly,
) (pub GatePublic, priv GatePrivate, err error) {

	n := len(A)           // #rows in A
	m := len(s)           // len(signature vector)
	lu := len(u)          // len(message block
	lx0 := len(x0)        // len(mask block
	k := m + 1 + lu + lx0 // witness length

	// -------------------------------------------------------------------------
	// 0)  Sanity-check dimensions
	// -------------------------------------------------------------------------
	if len(b1) != n || len(B0Const) != n {
		return pub, priv, errors.New("dimension mismatch in public vectors")
	}
	for _, row := range A {
		if len(row) != m {
			return pub, priv, errors.New("wrong A row length ≠ m")
		}
	}

	// -------------------------------------------------------------------------
	// 1)  Verify the proof-friendly equation
	// -------------------------------------------------------------------------
	// (b1 ⊙ A)·s
	left1 := make([]*ring.Poly, n)
	for j := 0; j < n; j++ {
		left1[j] = zeroPoly(ringQ)
		for t := 0; t < m; t++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(b1[j], A[j][t], tmp) // b₁ⱼ * Aⱼ,t
			ringQ.MulCoeffs(tmp, s[t], tmp)
			addInto(ringQ, left1[j], tmp)
		}
	}

	// (A·s) * x1
	left2 := make([]*ring.Poly, n)
	for j := 0; j < n; j++ {
		left2[j] = zeroPoly(ringQ)
		for t := 0; t < m; t++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(A[j][t], s[t], tmp) // Aⱼ,t * s_t
			ringQ.MulCoeffs(tmp, x1, tmp)
			addInto(ringQ, left2[j], tmp)
		}
	}

	// B0(1;u;x0)
	right := make([]*ring.Poly, n)
	one := zeroPoly(ringQ)
	one.Coeffs[0][0] = 1 // constant 1
	for j := 0; j < n; j++ {
		right[j] = copyPoly(ringQ, B0Const[j]) // 1 · B0const

		// + message part
		for i := 0; i < lu; i++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(B0Msg[i][j], u[i], tmp)
			addInto(ringQ, right[j], tmp)
		}
		// + randomness part
		for i := 0; i < lx0; i++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(B0Rnd[i][j], x0[i], tmp)
			addInto(ringQ, right[j], tmp)
		}
	}

	// Check equality row by row
	for j := 0; j < n; j++ {
		tmp := ringQ.NewPoly()
		ringQ.Sub(left1[j], left2[j], tmp) // (b⊙A)s − (A s)x1
		ringQ.Sub(tmp, right[j], tmp)      // − B0(...)
		if !ringQ.Equal(tmp, ringQ.NewPoly()) {
			fmt.Printf("Want 0 got %d\n", tmp.Coeffs[0][j])
			return pub, priv, fmt.Errorf("proof-friendly eq. fails on row %d", j)
		}
	}

	// -------------------------------------------------------------------------
	// 2)  Compute the compressed constants  ρᵀA, ρᵀ(b⊙A), ρᵀB0
	// -------------------------------------------------------------------------
	hatA := make([]*ring.Poly, m)
	hatb := make([]*ring.Poly, m)
	for t := 0; t < m; t++ {
		hatA[t] = zeroPoly(ringQ)
		hatb[t] = zeroPoly(ringQ)
		for j := 0; j < n; j++ {
			rho_j := rho[j] % ringQ.Modulus[0]
			tmp := mulConst(ringQ, A[j][t], rho_j)
			addInto(ringQ, hatA[t], tmp)

			// b1_j * A_jt first
			tmp2 := ringQ.NewPoly()
			ringQ.MulCoeffs(b1[j], A[j][t], tmp2)
			tmp2 = mulConst(ringQ, tmp2, rho_j)
			addInto(ringQ, hatb[t], tmp2)
		}
	}

	//  ρᵀ.B0 (split in three pieces)
	hatCconst := zeroPoly(ringQ)
	for j, col := range B0Const { // ← vector
		addInto(ringQ, hatCconst, mulConst(ringQ, col, rho[j]))
	}
	hatCmsg := make([]*ring.Poly, lu)
	for i := range hatCmsg {
		hatCmsg[i] = zeroPoly(ringQ)
		for j := 0; j < n; j++ {
			addInto(ringQ, hatCmsg[i], mulConst(ringQ, B0Msg[i][j], rho[j]))
		}
	}
	hatCrnd := make([]*ring.Poly, lx0)
	for i := range hatCrnd {
		hatCrnd[i] = zeroPoly(ringQ)
		for j := 0; j < n; j++ {
			addInto(ringQ, hatCrnd[i], mulConst(ringQ, B0Rnd[i][j], rho[j]))
		}
	}

	// -------------------------------------------------------------------------
	// 3)  Build  R2   (only off-diag block with −½ · hatA_t)
	// -------------------------------------------------------------------------
	pub.R2 = make([][]*ring.Poly, k)
	for i := 0; i < k; i++ {
		pub.R2[i] = make([]*ring.Poly, k)
		for j := 0; j < k; j++ {
			pub.R2[i][j] = zeroPoly(ringQ)
		}
	}
	q := ringQ.Modulus[0]
	inv2 := (q + 1) >> 1 // 1/2 mod q
	negHalf := q - inv2  // −½ mod q
	for t := 0; t < m; t++ {
		val := mulConst(ringQ, hatA[t], negHalf)
		pub.R2[t][m] = val
		pub.R2[m][t] = copyPoly(ringQ, val)
	}

	// -------------------------------------------------------------------------
	// 4)  Build  r1   (length-k)
	// -------------------------------------------------------------------------
	pub.R1 = make([]*ring.Poly, k)
	//  s-part
	for t := 0; t < m; t++ {
		pub.R1[t] = copyPoly(ringQ, hatb[t])
	}
	//  x1 entry
	pub.R1[m] = zeroPoly(ringQ)
	//  u block
	for i := 0; i < lu; i++ {
		pub.R1[m+1+i] = copyPoly(ringQ, hatCmsg[i])
		ringQ.Neg(pub.R1[m+1+i], pub.R1[m+1+i])
	}
	//  x0 block
	for i := 0; i < lx0; i++ {
		idx := m + 1 + lu + i
		pub.R1[idx] = copyPoly(ringQ, hatCrnd[i])
		ringQ.Neg(pub.R1[idx], pub.R1[idx])
	}

	// -------------------------------------------------------------------------
	// 5)  r0  = − hatCconst
	// -------------------------------------------------------------------------
	pub.R0 = copyPoly(ringQ, hatCconst)
	ringQ.Neg(pub.R0, pub.R0)

	// -------------------------------------------------------------------------
	// 6)  Pack witness vector
	// -------------------------------------------------------------------------
	priv.X = make([]*ring.Poly, 0, k)
	priv.X = append(priv.X, s...)
	priv.X = append(priv.X, x1)
	priv.X = append(priv.X, u...)
	priv.X = append(priv.X, x0...)

	return pub, priv, nil
}

// -----------------------------------------------------------------------------
//
//	VerifyGate – recompute  sᵀR₂s + r₁ᵀs + r₀  and check it is zero
//
// -----------------------------------------------------------------------------
func VerifyGate(r *ring.Ring, pub GatePublic, priv GatePrivate) error {

	k := len(priv.X)

	// --- 1. quadratic part  xᵀ R2 x -----------------------------------------
	quad := r.NewPoly()

	for i := 0; i < k; i++ {
		for j := 0; j < k; j++ {
			if r.Equal(pub.R2[i][j], r.NewPoly()) { // skip explicit zeros – sparse matrix
				continue
			}
			tmp := r.NewPoly()
			r.MulCoeffs(pub.R2[i][j], priv.X[i], tmp) // R2_ij * x_i
			r.MulCoeffs(tmp, priv.X[j], tmp)          // * x_j
			r.Add(quad, tmp, quad)
		}
	}

	// --- 2. linear part  r1ᵀ x ----------------------------------------------
	lin := r.NewPoly()
	for i := 0; i < k; i++ {
		if r.Equal(pub.R1[i], r.NewPoly()) {
			continue
		}
		tmp := r.NewPoly()
		r.MulCoeffs(pub.R1[i], priv.X[i], tmp)
		r.Add(lin, tmp, lin)
	}

	// --- 3. total = quad + lin + r0 ------------------------------------------
	total := r.NewPoly()
	r.Add(quad, lin, total)
	r.Add(total, pub.R0, total)

	// --- 4. success ? ---------------------------------------------------------
	if r.Equal(total, r.NewPoly()) {
		return nil
	}
	return fmt.Errorf("quadratic-gate check failed (first slot = %d)", total.Coeffs[0][0])
}
