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
// go test entry‑point
// --------------------------------------------------------------------------
func TestPACSSimulation(t *testing.T) {
	if !RunPACSSimulation() {
		t.Fatalf("verifier rejected – some check failed")
	}
}

// RunPACSSimulation executes one *interactive* proof and returns the verdict.
func RunPACSSimulation() bool {
	// ------------------------------------------------------------- parameters
	par, _ := signer.LoadParams("../Parameters/Parameters.json")
	ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})
	q := ringQ.Modulus[0]

	// ------------------------------------------------------------- witnesses
	w1, w2, w3 := BuildWitnessFromDisk() // helper in another PIOP file

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
	GammaP := sampleRandPolys(ringQ, rho, sCols, sCols) // Γ′ (deg ≤ s−1)
	gammaP := sampleRandMatrix(rho, 1, q)               // γ′ (single row)

	Fpar := buildFpar(ringQ, w1, w2, w3)
	A, b1, B0c, B0m, B0r := loadPublicTables(ringQ)
	Fagg := buildFagg(ringQ, w1, w2, A, b1, B0c, B0m, B0r)

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
// Helpers (unchanged or lightly patched)
// ============================================================================

func columnsToRows(r *ring.Ring, w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly, ell int) [][]uint64 {
	s := len(w1)
	ncols := r.N - ell
	rows := make([][]uint64, s+2) // s rows of w1, then x1 and w3
	q := r.Modulus[0]
	coeff := r.NewPoly()

	for row := 0; row < s; row++ {
		rows[row] = make([]uint64, ncols)
		for col, P := range w1 {
			r.InvNTT(P, coeff)
			rows[row][col] = coeff.Coeffs[0][row] % q
		}
		for col := s; col < ncols; col++ {
			rows[row][col] = 0
		}
	}

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
	for col := s; col < ncols; col++ {
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
		g := gammaP[i][0]
		f := evalPoint(r, Fagg[0], idx)
		rhs = modAdd(rhs, modMul(g, f, q), q)

		if lhs != rhs {
			return false
		}
	}
	return true
}

// ------- small utilities (unchanged) ---------------------------------------
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
