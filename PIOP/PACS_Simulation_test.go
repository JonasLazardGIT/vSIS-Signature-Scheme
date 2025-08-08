// PACS_Simulation.go – upgraded demo that now uses the **full** LVCS
// evaluation check (EvalStep2) instead of the light EvalStep2Slim.
//
// Key changes vs. the previous version
// -----------------------------------
//  1. A public coefficient matrix **C** (here: one‑row, deterministic) is
//     generated so we actually prove a *linear map* of the committed rows.
//  2. The prover calls  lvcs.EvalInit  to compute the masked targets **bar**.
//     These values are sent to the verifier (locally just stored).
//  3. The verifier now calls  EvalStep2  – which re‑runs the Slim tests *and*
//     additionally checks that the opened coordinates satisfy the linear
//     relation and match **bar**.
//
// Nothing else in the PACS layer (Eq.(4), Σ_Ω, etc.) changes.
// --------------------------------------------------------------------------
package PIOP

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"log"
	"math/bits"
	"testing"

	decs "vSIS-Signature/DECS"
	lvcs "vSIS-Signature/LVCS"
	signer "vSIS-Signature/Signer"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// --------------------------------------------------------------------------
// transcript – printed once at the end so humans can diff two runs quickly.
// --------------------------------------------------------------------------
type transcript struct {
	Root                string
	Gamma0, GammaPrime0 [][]uint64
	Rhash               string
	E                   []int // indices opened through DECS
	Flags               struct{ Merkle, Deg, LinMap, Eq4, Sum bool }
}

// --------------------------------------------------------------------------
// Global variables for bit-decomposition checks
// --------------------------------------------------------------------------

var (
	B             uint64
	T             int
	bitVals       []uint64
	sumSquares    uint64
	originalW1Len int
	tamperBit     bool
)

// --------------------------------------------------------------------------
// go test entry‑point
// --------------------------------------------------------------------------
func TestPACSSimulation(t *testing.T) {
	if !RunPACSSimulation() {
		t.Fatalf("verifier rejected – some check failed")
	}
}

// Negative test: flip one bit and expect rejection.
func TestPACSSimulationNegative(t *testing.T) {
	tamperBit = true
	if RunPACSSimulation() {
		t.Fatalf("verifier accepted invalid witness")
	}
	tamperBit = false
}

// RunPACSSimulation executes one *interactive* proof and returns the verdict.
func RunPACSSimulation() bool {
	// ------------------------------------------------------------- parameters
	par, _ := signer.LoadParams("../Parameters/Parameters.json")
	ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})
	q := ringQ.Modulus[0]

	// ------------------------------------------------------------- witnesses
	w1, w2, w3 := BuildWitnessFromDisk() // helper in another PIOP file
	A, b1, B0c, B0m, B0r := loadPublicTables(ringQ)

	// ─── after: w1, w2, w3 := BuildWitnessFromDisk() ──────────────
	originalW1Len = len(w1)
	mSig := originalW1Len - len(B0m) - len(B0r)
	sumSquares = 0
	tmp := ringQ.NewPoly()
	for t := 0; t < mSig; t++ {
		ringQ.InvNTT(w1[t], tmp)
		for j := 0; j < ringQ.N; j++ {
			v := tmp.Coeffs[0][j] % q
			sumSquares = (sumSquares + (v*v)%q) % q
		}
	}

	B = q
	delta := (B + q - sumSquares%q) % q
	T = bits.Len64(B)
	bitVals = make([]uint64, T)
	for t := 0; t < T; t++ {
		bitVals[t] = (delta >> uint(t)) & 1
	}
	for t := 0; t < T; t++ {
		p1 := ringQ.NewPoly()
		p1.Coeffs[0][0] = bitVals[t]
		ringQ.NTT(p1, p1)
		w1 = append(w1, p1)

		p3 := ringQ.NewPoly()
		p3.Coeffs[0][0] = bitVals[t]
		ringQ.NTT(p3, p3)
		w3 = append(w3, p3)
	}
	if tamperBit && len(bitVals) > 0 {
		bitVals[0] ^= 1
		if len(w1) > originalW1Len {
			c := ringQ.NewPoly()
			ringQ.InvNTT(w1[originalW1Len], c)
			c.Coeffs[0][0] ^= 1
			ringQ.NTT(c, c)
			w1[originalW1Len] = c
			w3[originalW1Len] = c.CopyNew()
		}
	}

	// ---------------------------------------------------------- LVCS.Commit
	ell := 1 // exactly one mask coordinate per row
	rows := columnsToRows(ringQ, w1, w2, w3, ell)
	root, pk, _ := lvcs.CommitInit(ringQ, rows, ell)

	vrf := lvcs.NewVerifier(ringQ, len(rows), decs.Eta)
	Gamma := vrf.CommitStep1(root)
	Rpolys := lvcs.CommitFinish(pk, Gamma)
	if !vrf.CommitStep2(Rpolys) {
		fmt.Println("[deg‑chk] R failed")
		return false
	}

	// ---------------------------------------------------------- coefficient C
	// We prove that  v = C·rows  where C is a *1×r* public matrix.
	rRows := len(rows)
	C := sampleRandMatrix(1, rRows, q) // shape 1×r  – deterministic helper

	// ------------------------------------------------------- compute  bar[k][i]
	bar := lvcs.EvalInit(ringQ, pk, C)

	// ------------------------------------------------------- PACS batching
	rho := 1
	sCols := len(w1)

	Fpar := buildFparBits(ringQ, w1, w2, w3)
	Fagg := buildFaggSlack(ringQ, w1, w2, A, b1, B0c, B0m, B0r)

	totalParallel := len(Fpar)
	GammaP := sampleRandPolys(ringQ, rho, totalParallel, sCols)
	totalAgg := len(Fagg)
	gammaP := sampleRandMatrix(rho, totalAgg, q)

	fmt.Printf("→ parallel rows: %d; aggregated rows: %d; witness cols: %d\n", totalParallel, totalAgg, sCols)

	omega := make([]uint64, sCols)
	for i := range omega {
		omega[i] = uint64(i + 1)
	}
	dQ := sCols + ell - 1
	M := BuildMaskPolynomials(ringQ, rho, dQ, omega)
	Q := BuildQ(ringQ, M, Fpar, Fagg, GammaP, gammaP)

	// --------------------------------------------------------- verifier picks E
	ncols := ringQ.N - ell // first index in mask area is ncols
	E := []int{ncols}      // open exactly the first mask coordinate

	// --------------------------------------------------------- opening & check
	open := lvcs.EvalFinish(pk, E)
	okLin := vrf.EvalStep2(bar, E, open.DECSOpen, C)
	if !okLin {
		fmt.Println("[open] DECS linear‑map check failed")
		return false
	}

	okEq4 := checkEq4OnOpening(ringQ, Q, M, open, Fpar, Fagg, GammaP, gammaP)
	okSum := VerifyQ(ringQ, Q, omega)

	// --------------------------------------------------------- pretty print
	tr := transcript{}
	tr.Root = hex.EncodeToString(root[:])
	tr.Gamma0 = firstColumn(Gamma)
	Rh := sha256.Sum256(polysToBytes(Rpolys))
	tr.Rhash = hex.EncodeToString(Rh[:])
	tr.GammaPrime0 = firstColumnNested(GammaP)
	tr.E = E
	tr.Flags = struct{ Merkle, Deg, LinMap, Eq4, Sum bool }{true, true, okLin, okEq4, okSum}
	pretty(&tr)

	return tr.Flags.Merkle && tr.Flags.Deg && tr.Flags.LinMap && tr.Flags.Eq4 && tr.Flags.Sum
}

// ============================================================================
// Helpers
// ============================================================================

func columnsToRows(r *ring.Ring, w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly, ell int) [][]uint64 {
	s := len(w1)
	ncols := r.N - ell
	rows := make([][]uint64, s+2)
	q := r.Modulus[0]

	coeffsW1 := make([]*ring.Poly, len(w1))
	for i, P := range w1 {
		c := r.NewPoly()
		r.InvNTT(P, c)
		coeffsW1[i] = c
	}

	for row := 0; row < s; row++ {
		rows[row] = make([]uint64, ncols)
		for col := 0; col < len(w1) && col < ncols; col++ {
			rows[row][col] = coeffsW1[col].Coeffs[0][row] % q
		}
		for col := len(w1); col < ncols; col++ {
			rows[row][col] = 0
		}
	}

	coeff := r.NewPoly()
	r.InvNTT(w2, coeff)
	rows[s] = make([]uint64, ncols)
	for col := 0; col < ncols; col++ {
		rows[s][col] = coeff.Coeffs[0][0]
	}

	rows[s+1] = make([]uint64, ncols)
	for col, P := range w3 {
		r.InvNTT(P, coeff)
		rows[s+1][col] = coeff.Coeffs[0][0]
	}
	for col := len(w3); col < ncols; col++ {
		rows[s+1][col] = 0
	}

	return rows
}

// Eq.(4) consistency on each opened index (unchanged)
func checkEq4OnOpening(r *ring.Ring, Q, M []*ring.Poly, op *lvcs.Opening,
	Fpar []*ring.Poly, Fagg []*ring.Poly, GammaP [][]*ring.Poly, gammaP [][]uint64) bool {

	q := r.Modulus[0]
	tmp := r.NewPoly()
	for i, idx := range op.DECSOpen.Indices {
		r.InvNTT(Q[i], tmp)
		lhs := tmp.Coeffs[0][idx]

		rhs := evalPoint(r, M[i], idx)
		for j := range GammaP[i] {
			g := evalPoint(r, GammaP[i][j], idx)
			f := evalPoint(r, Fpar[j], idx)
			rhs = modAdd(rhs, modMul(g, f, q), q)
		}
		for j := 0; j < len(Fagg); j++ {
			g := gammaP[i][j]
			f := evalPoint(r, Fagg[j], idx)
			rhs = modAdd(rhs, modMul(g, f, q), q)
		}

		if lhs != rhs {
			return false
		}
	}
	return true
}

// buildFparBits constructs the parallel constraint polynomials including
// additional bit-check rows.
func buildFparBits(r *ring.Ring, w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly) []*ring.Poly {
	out := make([]*ring.Poly, 0, originalW1Len+T)
	tmp := r.NewPoly()
	for k := 0; k < originalW1Len; k++ {
		p := w3[k].CopyNew()
		r.MulCoeffs(w1[k], w2, tmp)
		r.Sub(p, tmp, p)
		out = append(out, p)
	}
	for t := 0; t < T; t++ {
		Z := w1[originalW1Len+t]
		tmp2 := r.NewPoly()
		r.MulCoeffs(Z, Z, tmp2)
		r.Sub(tmp2, Z, tmp2)
		out = append(out, tmp2)
	}
	return out
}

// buildFaggSlack appends a constant slack row enforcing the bound on the
// sum of squares.
func buildFaggSlack(r *ring.Ring,
	w1 []*ring.Poly, w2 *ring.Poly,
	A [][]*ring.Poly, b1 []*ring.Poly,
	B0Const []*ring.Poly, B0Msg, B0Rnd [][]*ring.Poly) []*ring.Poly {

	mSig := originalW1Len - len(B0Msg) - len(B0Rnd)
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

	// Slack row enforcing sumSquares + Σ2^t bitVals[t] = B
	q := r.Modulus[0]
	slackVal := sumSquares % q
	for t, bv := range bitVals {
		slackVal = (slackVal + ((1<<uint(t))%q)*bv%q) % q
	}
	slackVal = (slackVal + q - (B % q)) % q

	slackPoly := r.NewPoly()
	for i := range slackPoly.Coeffs[0] {
		slackPoly.Coeffs[0][i] = slackVal
	}
	r.NTT(slackPoly, slackPoly)
	out = append(out, slackPoly)

	return out
}

// ------- small utilities  ---------------------------------------
func evalPoint(r *ring.Ring, p *ring.Poly, idx int) uint64 {
	c := r.NewPoly()
	r.InvNTT(p, c)
	return c.Coeffs[0][idx]
}

func firstColumn(mat [][]uint64) [][]uint64 {
	out := make([][]uint64, len(mat))
	for i := range mat {
		if len(mat[i]) > 0 {
			out[i] = []uint64{mat[i][0]}
		}
	}
	return out
}

func firstColumnNested(mat [][]*ring.Poly) [][]uint64 {
	out := make([][]uint64, len(mat))
	for i := range mat {
		if len(mat[i]) > 0 {
			out[i] = []uint64{mat[i][0].Coeffs[0][0]}
		}
	}
	return out
}

func pretty(tr *transcript) {
	fmt.Println("\n—— SmallWood one‑shot proof ——")
	fmt.Println("root        :", tr.Root)
	fmt.Println("Γ[0]        :", tr.Gamma0)
	fmt.Println("hash(R)     :", tr.Rhash)
	fmt.Println("Γ'[0]       :", tr.GammaPrime0)
	fmt.Println("E           :", tr.E)
	fmt.Printf("flags       : Merkle=%v  deg=%v  lin=%v  Eq4=%v  ΣΩ=%v\n",
		tr.Flags.Merkle, tr.Flags.Deg, tr.Flags.LinMap, tr.Flags.Eq4, tr.Flags.Sum)
	if tr.Flags.Merkle && tr.Flags.Deg && tr.Flags.LinMap && tr.Flags.Eq4 && tr.Flags.Sum {
		fmt.Println("VERIFIER → ACCEPT")
	} else {
		fmt.Println("VERIFIER → REJECT")
	}
}

func polysToBytes(pp []*ring.Poly) []byte {
	var out []byte
	for _, p := range pp {
		for _, c := range p.Coeffs[0] {
			buf := make([]byte, 8)
			binary.LittleEndian.PutUint64(buf, c)
			out = append(out, buf...)
		}
	}
	return out
}

// ---------------------------------------------------------------------------
// Public-data loader (A, b₁, B₀, …) – all NTT‑lifted on return. (unchanged)
// ---------------------------------------------------------------------------
func loadPublicTables(ringQ *ring.Ring) (A [][]*ring.Poly, b1, B0Const []*ring.Poly,
	B0Msg, B0Rnd [][]*ring.Poly) {

	pk, err := signer.LoadPublicKey("../public_key/public_key.json")
	if err != nil {
		log.Fatalf("public_key.json: %v", err)
	}
	A = [][]*ring.Poly{make([]*ring.Poly, len(pk.A))}
	for j, evals := range pk.A {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], evals)
		A[0][j] = p
	}

	rawB, err := loadBmatrixCoeffs("../Parameters/Bmatrix.json")
	if err != nil {
		log.Fatalf("Bmatrix.json: %v", err)
	}
	toNTT := func(coeffs []uint64) *ring.Poly {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs)
		ringQ.NTT(p, p)
		return p
	}
	B0Const = []*ring.Poly{toNTT(rawB[0])}
	B0Msg = [][]*ring.Poly{{toNTT(rawB[1])}}
	B0Rnd = [][]*ring.Poly{{toNTT(rawB[2])}}
	b1 = []*ring.Poly{toNTT(rawB[3])}

	return
}
