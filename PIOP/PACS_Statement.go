// quadratic_gate_pacs.go
package PIOP

import (
	"encoding/binary"
	"fmt"
	"log"
	"time"
	signer "vSIS-Signature/Signer"
	prof "vSIS-Signature/prof"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// -----------------------------------------------------------------------------
// Θ   – all-zero constants for the parallel constraint f
// -----------------------------------------------------------------------------

// BuildTheta returns a slice of *ring.Poly – one per witness column – that
// represent θ_{j,i,k}.  In our quadratic gate **all θ are zero**, so each
// polynomial is simply the zero-poly in NTT domain.
//
// The caller supplies s = len(w1) so the slice matches the number of columns.
func BuildTheta(ringQ *ring.Ring, s int) []*ring.Poly {
	return BuildThetaZeros(ringQ, s) // helper defined earlier
}

// -----------------------------------------------------------------------------
//  Small utility – evaluate a coefficient-domain polynomial  mod q
// -----------------------------------------------------------------------------

// evalPoly returns   P(x)  where  P(X)=Σ_i coeffs[i]·X^i  and all arithmetic
// is modulo q.  The slice is in ascending degree order (coeffs[0] = a₀).
func EvalPoly(coeffs []uint64, x, q uint64) uint64 {
	if len(coeffs) == 0 {
		return 0
	}
	// Horner scheme: (((a_d)·x + a_{d-1})·x + … + a₀)
	res := coeffs[len(coeffs)-1] % q
	for i := len(coeffs) - 2; i >= 0; i-- {
		res = modMul(res, x, q)
		res = modAdd(res, coeffs[i]%q, q)
	}
	return res
}

// -----------------------------------------------------------------------------
// Θ – all–zero parallel-constraint constants
// -----------------------------------------------------------------------------

// BuildThetaZeros returns s constant-zero polys (degree 0) in NTT form.
// They materialise θ_{j,i,k}=0 for every i,k in our quadratic gate.
func BuildThetaZeros(ringQ *ring.Ring, s int) []*ring.Poly {
	zeros := make([]*ring.Poly, s)
	z := ringQ.NewPoly() // already 0 in coeff & NTT
	for i := 0; i < s; i++ {
		zeros[i] = z.CopyNew()
	}
	return zeros
}

// -----------------------------------------------------------------------------
// Θ′ – interpolating polys for public coefficients in f′
// -----------------------------------------------------------------------------

// BuildThetaPrime constructs a polynomial for every public coefficient that
// depends on the column index k:  A_{j,k},  b1_j,  B0msg_{i,k},  B0rnd_{i,k}.
//
//	values[k] must be the table of that coefficient for k=0..s-1.
//	omega     is the evaluation set Ω of size s.
//
// The result is an NTT poly whose degree ≤ s-1 and which satisfies
//
//	P(omega[k]) = values[k]  for all k.
func BuildThetaPrime(ringQ *ring.Ring, values, omega []uint64) *ring.Poly {
	if len(values) != len(omega) {
		panic("BuildThetaPrime: length mismatch")
	}
	q := ringQ.Modulus[0]
	coeffs := Interpolate(omega, values, q)
	p := ringQ.NewPoly()
	copy(p.Coeffs[0], coeffs)
	ringQ.NTT(p, p)
	return p
}

// -----------------------------------------------------------------------------
// Pure evaluators for  f  and  f′
// -----------------------------------------------------------------------------

// EvalParallel returns   w3_k - w1_k*w2   ∈ F_q
func EvalParallel(ringQ *ring.Ring, w1k, w2, w3k uint64) uint64 {
	q := ringQ.Modulus[0]
	return modSub(w3k, modMul(w1k, w2, q), q)
}

// EvalAggregated returns
// (b1·A)s - (A·s)x1 - B0(1;u;x0)   for one *row* j   (mod q).
//
// Inputs are already the field sums for that row.
func EvalAggregated(ringQ *ring.Ring, term1, term2, B0 uint64) uint64 {
	q := ringQ.Modulus[0]
	return modSub(modSub(term1, term2, q), B0, q)
}

// -----------------------------------------------------------------------------
// Θ′  – interpolating polys for every public coefficient appearing in f′
// -----------------------------------------------------------------------------

type ThetaPrime struct {
	ARows   [][]*ring.Poly // [row][col]
	B1Rows  []*ring.Poly   // one polynomial per row  j (|Ω| evaluations)
	B0Const []*ring.Poly   // idem
	B0Msg   [][]*ring.Poly // [msgIdx][row]
	B0Rnd   [][]*ring.Poly // [rndIdx][row]
}

// BuildThetaPrimeSet builds Θ′ for *all* public tables.  omega is the
// evaluation set Ω = {ω₁,…,ω_s} (length s) used when you interpolated witness
// rows into P_i(X).
func BuildThetaPrimeSet(
	ringQ *ring.Ring,
	A [][]*ring.Poly, // [row][col]
	b1 []*ring.Poly, // [row]
	B0Const []*ring.Poly, // [row]
	B0Msg, B0Rnd [][]*ring.Poly, // [msgIdx][row]  /  [rndIdx][row]
	omega []uint64, // |Ω| = s
) *ThetaPrime {

	q := ringQ.Modulus[0]
	s := len(omega)

	// -- A rows/cols ----------------------------------------------------------
	aRows := make([][]*ring.Poly, len(A))
	for i := range A {
		aRows[i] = make([]*ring.Poly, len(A[i]))
		for k := range A[i] {
			coeff := ringQ.NewPoly()
			ringQ.InvNTT(A[i][k], coeff)
			vals := make([]uint64, s)
			for j := 0; j < s; j++ {
				vals[j] = EvalPoly(coeff.Coeffs[0], omega[j]%q, q)
			}
			aRows[i][k] = BuildThetaPrime(ringQ, vals, omega)
		}
	}

	// -- helper for row-wise constants ----------------------------------------
	buildRowPolys := func(src []*ring.Poly) []*ring.Poly {
		out := make([]*ring.Poly, len(src))
		for j, pj := range src {
			coeff := ringQ.NewPoly()
			ringQ.InvNTT(pj, coeff)
			vals := make([]uint64, s)
			for t := 0; t < s; t++ {
				vals[t] = EvalPoly(coeff.Coeffs[0], omega[t]%q, q)
			}
			out[j] = BuildThetaPrime(ringQ, vals, omega)
		}
		return out
	}

	b1Rows := buildRowPolys(b1)
	b0cRows := buildRowPolys(B0Const)

	build2D := func(src [][]*ring.Poly) [][]*ring.Poly {
		out := make([][]*ring.Poly, len(src))
		for idx := range src {
			out[idx] = buildRowPolys(src[idx])
		}
		return out
	}

	return &ThetaPrime{
		ARows:   aRows,
		B1Rows:  b1Rows,
		B0Const: b0cRows,
		B0Msg:   build2D(B0Msg),
		B0Rnd:   build2D(B0Rnd),
	}
}

// q_polys.go  (same PIOP package)

// BuildQ constructs the vector of polynomials Q_i(X) as in Eq.(4).
//
// inputs:
//
//	M          – slice of *ring.Poly masking rows  (len = rho)
//	Fpar       – slice of *ring.Poly for every parallel constraint F_j(X)
//	Fagg       – slice of *ring.Poly for every aggregated constraint F'_j(X)
//	GammaPrime – [][]*ring.Poly   len(rho) × m1    (random deg≤s-1)
//	gammaPrime – [][]uint64       len(rho) × m2    (random scalars)
//
// output:  []*ring.Poly   len = rho
func BuildQ(
	ringQ *ring.Ring,
	M []*ring.Poly,
	Fpar []*ring.Poly,
	Fagg []*ring.Poly,
	GammaPrime [][]uint64,
	gammaPrime [][]uint64,
) []*ring.Poly {
	defer prof.Track(time.Now(), "BuildQ")

	rho := len(M)
	m1 := len(Fpar)
	m2 := len(Fagg)

	Q := make([]*ring.Poly, rho)
	tmp := ringQ.NewPoly()
	for i := 0; i < rho; i++ {
		Qi := M[i].CopyNew() // start with M_i(X)

		// Σ_j Γ'_{i,j} * F_j(X)
		for j := 0; j < m1; j++ {
			mulScalarNTT(ringQ, Fpar[j], GammaPrime[i][j], tmp)
			addInto(ringQ, Qi, tmp)
		}
		// Σ_j γ'_{i,j} * F'_j(X)
		for j := 0; j < m2; j++ {
			mulScalarNTT(ringQ, Fagg[j], gammaPrime[i][j], tmp)
			addInto(ringQ, Qi, tmp)
		}
		Q[i] = Qi
	}
	return Q
}

// verify_q.go

// VerifyQ checks, for every i∈[ρ], that
//
//	Σ_{ω∈Ω} Q_i(ω) = 0   (Eq.(7))
//
// The caller has already performed the Merkle-consistency check.
func VerifyQ(
	ringQ *ring.Ring,
	Q []*ring.Poly,
	omega []uint64,
) bool {
	defer prof.Track(time.Now(), "VerifyQ")
	coeff := ringQ.NewPoly()
	q := ringQ.Modulus[0]

	seen := make(map[uint64]struct{}, len(omega))
	for _, w := range omega {
		wm := w % q
		if _, ok := seen[wm]; ok {
			log.Fatalf("VerifyQ: Ω has duplicate element %d", wm)
		}
		seen[wm] = struct{}{}
	}
	if q%uint64(len(omega)) == 0 {
		log.Fatalf("VerifyQ: q (= %d) is multiple of |Ω| (= %d)", q, len(omega))
	}

	for i, Qi := range Q {
		ringQ.InvNTT(Qi, coeff)
		sum := uint64(0)
		for _, w := range omega {
			sum = modAdd(sum, EvalPoly(coeff.Coeffs[0], w, q), q)
		}
		if DEBUG_SUMS {
			fmt.Printf("[Q %d] ΣΩ Q_i = %d\n", i, sum)
		}
		if sum != 0 {
			return false
		}
	}
	return true
}

// -----------------------------------------------------------------------------
//  Full PACS verification  (f and f′ exactly as in the paper)
// -----------------------------------------------------------------------------

// VerifyFullPACS checks the two sets of constraints verbatim:
//
//   - ∀k :  f ( w1_k , w2 , w3_k ) = 0                    (parallel)
//   - ∀j :  Σ_k f′( column_k , θ′_j,k ) = 0              (aggregated)
//
// It returns true iff both families hold.
func VerifyFullPACS(
	ringQ *ring.Ring,
	w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly,
	A [][]*ring.Poly, b1 []*ring.Poly,
	B0Const []*ring.Poly, B0Msg, B0Rnd [][]*ring.Poly,
) bool {

	s := len(w1) // number of columns

	// ------------------------------------------------------------------ f ----
	tmp := ringQ.NewPoly()
	zero := ringQ.NewPoly()
	for k := 0; k < s; k++ {
		ringQ.MulCoeffs(w1[k], w2, tmp)
		ringQ.Sub(w3[k], tmp, tmp)
		if !ringQ.Equal(tmp, zero) {
			return false // f constraint breaks for this k
		}
	}

	// -------------------------------------------------------------- f′ summed
	mSig := len(w1) - len(B0Msg) - len(B0Rnd) // length of signature vector s
	if mSig < 0 {
		log.Fatalf("VerifyFullPACS: negative mSig (len(w1)=%d)", len(w1))
	}
	nRows := len(A)

	left1 := ringQ.NewPoly()
	left2 := ringQ.NewPoly()
	right := ringQ.NewPoly()

	for j := 0; j < nRows; j++ {
		// reset accumulators
		for _, p := range []*ring.Poly{left1, left2, right} {
			resetPoly(p)
		}

		// Σ_k (b1⊙A)_j·s
		for t := 0; t < mSig; t++ {
			ringQ.MulCoeffs(b1[j], A[j][t], tmp)
			ringQ.MulCoeffs(tmp, w1[t], tmp)
			addInto(ringQ, left1, tmp)
		}
		// Σ_k (A s) * x1
		for t := 0; t < mSig; t++ {
			ringQ.MulCoeffs(A[j][t], w1[t], tmp)
			ringQ.MulCoeffs(tmp, w2, tmp)
			addInto(ringQ, left2, tmp)
		}
		// Σ_k B0(1;u;x0)
		addInto(ringQ, right, B0Const[j])
		for i := range B0Msg {
			ringQ.MulCoeffs(B0Msg[i][j], w1[mSig+i], tmp)
			addInto(ringQ, right, tmp)
		}
		offset := mSig + len(B0Msg)
		for i := range B0Rnd {
			ringQ.MulCoeffs(B0Rnd[i][j], w1[offset+i], tmp)
			addInto(ringQ, right, tmp)
		}

		// f′_sum = left1 - left2 - right  must be zero
		ringQ.Sub(left1, left2, tmp)
		ringQ.Sub(tmp, right, tmp)
		if !ringQ.Equal(tmp, zero) {
			return false // some aggregated constraint failed
		}
	}
	return true
}

// BuildQFromDisk loads the same JSON fixtures used by VerifyGHFromDisk,
// samples fresh Fiat-Shamir randomness Γ' and γ', computes
//
//	Q = (Q₁,…,Q_ρ)
//
// as in Eq.(4) of the SmallWood paper, and returns it.
//
// The function is **self-contained**: you can call it from a unit-test and
// check the Ω-sum condition with VerifyQ (see q_polys.go).
func BuildQFromDisk() (Q []*ring.Poly, omega []uint64, ringQ *ring.Ring) {

	//-------------------------------------------------------------------[0] ring
	par, err := signer.LoadParams("../Parameters/Parameters.json")
	if err != nil {
		log.Fatalf("Parameters.json: %v", err)
	}
	ringQ, _ = ring.NewRing(par.N, []uint64{par.Q})
	toNTT := func(p *ring.Poly) { ringQ.NTT(p, p) }

	//-------------------------------------------------------------------[1] A,pk
	pk, err := signer.LoadPublicKey("../public_key/public_key.json")
	if err != nil {
		log.Fatalf("public_key.json: %v", err)
	}
	A := [][]*ring.Poly{make([]*ring.Poly, len(pk.A))}
	for i, c := range pk.A {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], c)
		A[0][i] = p
	}

	//-------------------------------------------------------------------[2] B-mat
	Bcoeffs, _ := loadBmatrixCoeffs("../Parameters/Bmatrix.json")
	B0Const := []*ring.Poly{toNTTwrap(ringQ, Bcoeffs[0], toNTT)}
	B0Msg := [][]*ring.Poly{{toNTTwrap(ringQ, Bcoeffs[1], toNTT)}}
	B0Rnd := [][]*ring.Poly{{toNTTwrap(ringQ, Bcoeffs[2], toNTT)}}
	b1 := []*ring.Poly{toNTTwrap(ringQ, Bcoeffs[3], toNTT)}

	//-------------------------------------------------------------------[3] sign
	sig, _ := loadSignature("../Signature/Signature.json")
	m := toNTTwrap(ringQ, sig.Message, toNTT)
	x0 := toNTTwrap(ringQ, sig.X0, toNTT)
	x1 := toNTTwrap(ringQ, sig.X1, toNTT)
	s := makeSigPolys(ringQ, sig.Signature)

	//-------------------------------------------------------------------[4] witness
	w1, w2, w3 := BuildWitness(ringQ,
		A, b1, B0Const, B0Msg, B0Rnd,
		/*private*/ s, x1, []*ring.Poly{m}, []*ring.Poly{x0})

	//-------------------------------------------------------------------[5] Ω  = first s points of the NTT evaluation grid
	sCols := len(w1)
	px := ringQ.NewPoly() // P(X)=X
	px.Coeffs[0][1] = 1
	pts := ringQ.NewPoly()
	ringQ.NTT(px, pts) // NTT(X) enumerates the grid (order consistent with ringQ)
	omega = make([]uint64, sCols)
	copy(omega, pts.Coeffs[0][:sCols])

	//-------------------------------------------------------------------[6] build F_j(X) & F'_j(X)
	Fpar := buildFpar(ringQ, w1, w2, w3) // len = s
	Fagg := buildFagg(ringQ, w1, w2,
		A, b1, B0Const, B0Msg, B0Rnd) // len = nRows(=1)

	//-------------------------------------------------------------------[7] sample Γ' , γ'
	rho := 1 // one masking row is enough here
	ell := 1 // expose 1 random point per row
	dQ := sCols + ell - 1

	// Bind to public inputs via FS: include Ω and public tables
	q := ringQ.Modulus[0]
	concatPolys := func(pp []*ring.Poly) []byte {
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
	fsOmg := bytesU64Vec(omega)
	fsA := concatPolys(A[0])
	fsB1 := concatPolys(b1)
	fsGamma := newFSRNG("GammaPrime:offline", fsOmg, fsA, fsB1)
	fsGammaSc := newFSRNG("gammaPrime:offline", fsOmg, fsA, fsB1)
	GammaPrime := sampleFSMatrix(rho, len(Fpar), q, fsGamma)
	gammaPrime := sampleFSMatrix(rho, len(Fagg), q, fsGammaSc)

	// precompute Ω-sums of Fpar and Fagg
	sumFpar := sumPolyList(ringQ, Fpar, omega)
	sumFagg := sumPolyList(ringQ, Fagg, omega)

	M := BuildMaskPolynomials(ringQ, rho, dQ, omega, GammaPrime, gammaPrime, sumFpar, sumFagg)

	//-------------------------------------------------------------------[8] build Q
	Q = BuildQ(ringQ, M, Fpar, Fagg, GammaPrime, gammaPrime)

	fmt.Println("[BuildQFromDisk]  built", len(Q), "polys  (deg ≤", dQ, ")")
	return Q, omega, ringQ
}

// ---------- tiny helpers -----------------------------------------------------

func toNTTwrap(r *ring.Ring, coeffs []uint64, lift func(*ring.Poly)) *ring.Poly {
	p := r.NewPoly()
	copy(p.Coeffs[0], coeffs)
	lift(p)
	return p
}
func makeSigPolys(r *ring.Ring, rows [][]uint64) []*ring.Poly {
	out := make([]*ring.Poly, len(rows))
	for i, c := range rows {
		p := r.NewPoly()
		copy(p.Coeffs[0], c)
		out[i] = p
	}
	return out
}

// Fpar_k(X)  = w3_k - w1_k·w2
func buildFpar(r *ring.Ring, w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly) []*ring.Poly {
	defer prof.Track(time.Now(), "buildFpar")
	out := make([]*ring.Poly, len(w1))
	for k := range w1 {
		out[k] = makeProductConstraint(r, w1[k], w2, w3[k])
	}
	return out
}

// Fagg_j(X)  = (b1⊙A)s − (A·s)x1 − B0(1;u;x0)
func buildFagg(r *ring.Ring,
	w1 []*ring.Poly, w2 *ring.Poly,
	A [][]*ring.Poly, b1 []*ring.Poly,
	B0Const []*ring.Poly, B0Msg, B0Rnd [][]*ring.Poly) []*ring.Poly {
	defer prof.Track(time.Now(), "buildFagg")

	mSig := len(w1) - len(B0Msg) - len(B0Rnd)
	out := make([]*ring.Poly, len(A))
	tmp := r.NewPoly()
	left1 := r.NewPoly()
	left2 := r.NewPoly()
	right := r.NewPoly()

	for j := range A {
		clear := func(p *ring.Poly) {
			for i := range p.Coeffs[0] {
				p.Coeffs[0][i] = 0
			}
		}
		clear(left1)
		clear(left2)
		clear(right)

		for t := 0; t < mSig; t++ {
			r.MulCoeffs(b1[j], A[j][t], tmp)
			r.MulCoeffs(tmp, w1[t], tmp)
			addInto(r, left1, tmp)

			r.MulCoeffs(A[j][t], w1[t], tmp)
			r.MulCoeffs(tmp, w2, tmp)
			addInto(r, left2, tmp)
		}
		addInto(r, right, B0Const[j])
		for i := range B0Msg {
			r.MulCoeffs(B0Msg[i][j], w1[mSig+i], tmp)
			addInto(r, right, tmp)
		}
		offset := mSig + len(B0Msg)
		for i := range B0Rnd {
			r.MulCoeffs(B0Rnd[i][j], w1[offset+i], tmp)
			addInto(r, right, tmp)
		}

		r.Sub(left1, left2, tmp)
		r.Sub(tmp, right, tmp)
		out[j] = tmp.CopyNew()
	}
	return out
}

// buildFaggOnOmega – BBS aggregated rows (Eq. 3) using Θ′(X) on Ω.
//
// Fagg_j(X) = (b1⊙A)_j(X) · S(X)  −  (A_j(X)·S(X))·X1(X)  −  B0_j(X)
func buildFaggOnOmega(
	r *ring.Ring,
	w1 []*ring.Poly, w2 *ring.Poly,
	theta *ThetaPrime, // built with BuildThetaPrimeSet(..., omega)
	mSig int,
) []*ring.Poly {
	defer prof.Track(time.Now(), "buildFaggOnOmega")

	out := make([]*ring.Poly, len(theta.ARows))
	tmp := r.NewPoly()
	left1 := r.NewPoly()
	left2 := r.NewPoly()
	right := r.NewPoly()

	for j := range theta.ARows {
		// clear accumulators
		resetPoly(left1)
		resetPoly(left2)
		resetPoly(right)

		// Σ_t  (b1⊙A)_j,t  *  s_t
		for t := 0; t < mSig; t++ {
			r.MulCoeffs(theta.B1Rows[j], theta.ARows[j][t], tmp) // b1_j * A_j,t
			r.MulCoeffs(tmp, w1[t], tmp)
			addInto(r, left1, tmp)
		}
		// Σ_t  (A_j,t * s_t) * x1
		for t := 0; t < mSig; t++ {
			r.MulCoeffs(theta.ARows[j][t], w1[t], tmp)
			r.MulCoeffs(tmp, w2, tmp)
			addInto(r, left2, tmp)
		}
		// B0 · (1; u; x0)  : constant + message + randomness blocks
		addInto(r, right, theta.B0Const[j])
		// message block
		for i := range theta.B0Msg {
			r.MulCoeffs(theta.B0Msg[i][j], w1[mSig+i], tmp)
			addInto(r, right, tmp)
		}
		// randomness block
		off := mSig + len(theta.B0Msg)
		for i := range theta.B0Rnd {
			r.MulCoeffs(theta.B0Rnd[i][j], w1[off+i], tmp)
			addInto(r, right, tmp)
		}

		// F'_j(X) = left1 - left2 - right
		r.Sub(left1, left2, tmp)
		r.Sub(tmp, right, tmp)
		out[j] = tmp.CopyNew()
	}
	return out
}

// random deg ≤ s-1 polys
func sampleRandPolys(r *ring.Ring, rows, cols, s int) [][]*ring.Poly {
	defer prof.Track(time.Now(), "sampleRandPolys")
	out := make([][]*ring.Poly, rows)
	for i := 0; i < rows; i++ {
		out[i] = make([]*ring.Poly, cols)
		for j := 0; j < cols; j++ {
			p := r.NewPoly()
			for k := 0; k < s; k++ {
				p.Coeffs[0][k] = randUint64Mod(r.Modulus[0])
			}
			r.NTT(p, p)
			out[i][j] = p
		}
	}
	return out
}

// random scalar matrix  rows×cols  in F_q
func sampleRandMatrix(rows, cols int, q uint64) [][]uint64 {
	defer prof.Track(time.Now(), "sampleRandMatrix")
	M := make([][]uint64, rows)
	for i := 0; i < rows; i++ {
		M[i] = make([]uint64, cols)
		for j := 0; j < cols; j++ {
			M[i][j] = randUint64Mod(q)
		}
	}
	return M
}
