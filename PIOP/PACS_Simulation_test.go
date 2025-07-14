// PACS_Simulation.go – **reference implementation** of the full SmallWood
// proof‑system stack on the quadratic‑gate demo instance.
// ---------------------------------------------------------------------------
// Layer diagram (each arrow is an *actual* function call in this file):
//
//	prover                                verifier
//	─────────────────────────────────────────────────────────────────────────
//	  rows ─────► lvcs.CommitInit ───► Merkle root (DECS)
//	                    │                        │
//	                    │  Γ (degree test)       │
//	                    ├────── CommitStep1 ◄────┤   derives Γ from root
//	                    │                        │
//	R‑polys ◄────────────┼── CommitFinish ◄──────┤   degree‑check on R
//	                    │                        │
//	           PACS batching (Γ', γ')
//	                    │                        │
//	build Q             ├──────────── Q ─────────►┤
//	                    │                        │
//	          choose E' │                        │ choose E', open P,M
//
// lvcs.EvalFinish ──────┼────────── open ───────►┤ Merkle+deg on opening
//
//	     │                        │ Eq.(4) on E'  +  ΣΩQ=0
//	     ▼                        ▼
//	ACCEPT / REJECT            ACCEPT / REJECT
//
// ---------------------------------------------------------------------------
// Package expectations
// --------------------
//   - **decs** – Merkle‑tree verification + exported `Degree` constant.
//   - **lvcs** – API:
//     CommitInit / CommitFinish / EvalFinish, NewVerifier, …
//     `ProverKey` must expose:  `MaskPolys` and `RowPolys` (both []*ring.Poly).
//     `Opening`   must embed:   `DECSOpen *decs.DECSOpening`.
//
// ---------------------------------------------------------------------------
//
//	go test -run TestPACSSimulation   # CI / regression
//	go run ./PIOP                    # manual transcript
//
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
//   Minimal transcript object – printed at the end for visual inspection
// ---------------------------------------------------------------------------

type transcript struct {
	Root                string
	Gamma0, GammaPrime0 [][]uint64 // only first column shown
	Rhash               string
	Eprime, Edecs       []int
	Flags               struct{ Merkle, Deg, Eq4, Sum bool }
}

// ---------------------------------------------------------------------------
//   Public entry‑point for `go test`
// ---------------------------------------------------------------------------

func TestPACSSimulation(t *testing.T) {
	if !RunPACSSimulation() {
		t.Fatalf("verifier rejected – at least one check failed")
	}
}

// RunPACSSimulation executes one end‑to‑end proof and returns the verdict.
func RunPACSSimulation() bool {
	// 0) ring / parameters --------------------------------------------------
	par, _ := signer.LoadParams("../Parameters/Parameters.json")
	ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})

	// 1) witness columns (quadratic‑gate helper) ---------------------------
	w1, w2, w3 := BuildWitnessFromDisk()

	// 2) convert to rows and run LVCS.Commit -------------------------------
	rows := columnsToRows(ringQ, w1, w2, w3)
	ell := 1
	root, pk, _ := lvcs.CommitInit(ringQ, rows, ell)

	vrf := lvcs.NewVerifier(ringQ, len(rows), decs.Eta)
	Gamma := vrf.CommitStep1(root)
	Rpolys := lvcs.CommitFinish(pk, Gamma)
	if !vrf.CommitStep2(Rpolys) {
		fmt.Println("[deg‑chk] R failed")
		return false
	}

	// 3) PACS batching randomness -----------------------------------------
	rho := 1
	sCols := len(w1)
	GammaP := sampleRandPolys(ringQ, rho, sCols, sCols)  // Γ′
	gammaP := sampleRandMatrix(rho, 1, ringQ.Modulus[0]) // γ′, m₂=1

	// 4) build Fₚₐᵣ and Fₐgg from witness columns
	Fpar := buildFpar(ringQ, w1, w2, w3)

	// Toy instance ⇒ one aggregated row j=0
	A, b1, B0c, B0m, B0r := loadPublicTables(ringQ)
	Fagg := buildFagg(ringQ, w1, w2, A, b1, B0c, B0m, B0r)
	// 5) Build Qᵢ(X) -------------------------------------------------------
	omega := make([]uint64, sCols)
	for i := range omega {
		omega[i] = uint64(i + 1)
	}

	dQ := sCols + ell - 1
	M := BuildMaskPolynomials(ringQ, rho, dQ, omega)
	Q := BuildQ(ringQ, M, Fpar, Fagg, GammaP, gammaP)
	// 6) choose E′, open P & M --------------------------------------------
	Eprime := []int{sCols + 2}
	open := lvcs.EvalFinish(pk, Eprime)
	if !vrf.EvalStep2Slim(open.DECSOpen) { // Merkle + deg on opening
		fmt.Println("[open] Merkle/deg failed")
		return false
	}

	okEq4 := checkEq4OnOpening(ringQ, Q, M, open, Fpar, Fagg, GammaP, gammaP)
	okSum := VerifyQ(ringQ, Q, omega)

	// 7) print transcript ---------------------------------------------------
	tr := transcript{}
	tr.Root = hex.EncodeToString(root[:])
	tr.Gamma0 = firstColumn(Gamma)
	Rh := sha256.Sum256(polysToBytes(Rpolys))
	tr.Rhash = hex.EncodeToString(Rh[:])
	tr.GammaPrime0 = firstColumnNested(GammaP)
	tr.Eprime, tr.Edecs = Eprime, open.DECSOpen.Indices
	tr.Flags = struct{ Merkle, Deg, Eq4, Sum bool }{true, true, okEq4, okSum}
	pretty(&tr)

	return tr.Flags.Merkle && tr.Flags.Deg && tr.Flags.Eq4 && tr.Flags.Sum
}

// ----------------------------------------------------------------------------
// columnsToRows – helper: transform witness columns into LVCS row matrix
// ----------------------------------------------------------------------------

func columnsToRows(r *ring.Ring, w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly) [][]uint64 {
	s := len(w1)
	rows := make([][]uint64, s+2) // s rows of w1, then x1 and w3
	q := r.Modulus[0]
	coeff := r.NewPoly()

	for row := 0; row < s; row++ {
		rows[row] = make([]uint64, s)
		for col, P := range w1 {
			r.InvNTT(P, coeff)
			rows[row][col] = coeff.Coeffs[0][row] % q
		}
	}
	r.InvNTT(w2, coeff)
	rows[s] = make([]uint64, s)
	for col := range rows[s] {
		rows[s][col] = coeff.Coeffs[0][0]
	}

	rows[s+1] = make([]uint64, s)
	for col, P := range w3 {
		r.InvNTT(P, coeff)
		rows[s+1][col] = coeff.Coeffs[0][0]
	}
	return rows
}

// ----------------------------------------------------------------------------
// Eq.(4) on opened evaluations (E′)
// ----------------------------------------------------------------------------

func checkEq4OnOpening(r *ring.Ring, Q, M []*ring.Poly, op *lvcs.Opening,
	Fpar []*ring.Poly, Fagg []*ring.Poly,
	GammaP [][]*ring.Poly, gammaP [][]uint64) bool {

	q := r.Modulus[0]
	coeff := r.NewPoly()
	for i, idx := range op.DECSOpen.Indices {
		r.InvNTT(Q[i], coeff)
		lhs := coeff.Coeffs[0][idx]
		rhs := evalPoint(r, M[i], idx)
		for j := range GammaP[i] {
			gp := evalPoint(r, GammaP[i][j], idx)
			fp := evalPoint(r, Fpar[j], idx)
			rhs = modAdd(rhs, modMul(gp, fp, q), q)
		}
		g := gammaP[i][0]
		fp := evalPoint(r, Fagg[0], idx)
		rhs = modAdd(rhs, modMul(g, fp, q), q)

		if lhs != rhs {
			return false
		}
	}
	return true
}

// ----------------------------------------------------------------------------
// tiny helpers (evaluation / pretty print / misc.)
// ----------------------------------------------------------------------------

func evalPoint(r *ring.Ring, p *ring.Poly, idx int) uint64 {
	coeff := r.NewPoly()
	r.InvNTT(p, coeff)
	return coeff.Coeffs[0][idx]
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
	fmt.Println("E' / E      :", tr.Eprime, "/", tr.Edecs)
	fmt.Printf("flags       : Merkle=%v  deg=%v  Eq4=%v  ΣΩ=%v\n",
		tr.Flags.Merkle, tr.Flags.Deg, tr.Flags.Eq4, tr.Flags.Sum)
	if tr.Flags.Merkle && tr.Flags.Deg && tr.Flags.Eq4 && tr.Flags.Sum {
		fmt.Println("VERIFIER → ACCEPT")
	} else {
		fmt.Println("VERIFIER → REJECT")
	}
}

func polysToBytes(pp []*ring.Poly) []byte {
	var out []byte
	for _, p := range pp {
		for _, c := range p.Coeffs[0] {
			b := make([]byte, 8)
			binary.LittleEndian.PutUint64(b, c)
			out = append(out, b...)
		}
	}
	return out
}

// ----------------------------------------------------------------------------
// External helpers reused from other PIOP files (just prototypes here)
// ----------------------------------------------------------------------------

// loadPublicTables reads the JSON fixtures shipping with the demo and
// returns **NTT-lifted** public data needed by the quadratic-gate checks.
//
// A        : one-row slice  [row][col]   (NTT domain)
// b1       : []*ring.Poly   (len = nRows)
// B0Const  : []*ring.Poly   ― constant column of B₀
// B0Msg    : [][]*ring.Poly ― message columns  (outer idx = msg block)
// B0Rnd    : [][]*ring.Poly ― randomness cols  (outer idx = rnd block)
//
// Panics on I/O or shape errors – fine for test-only helpers.
func loadPublicTables(ringQ *ring.Ring) (A [][]*ring.Poly,
	b1, B0Const []*ring.Poly,
	B0Msg, B0Rnd [][]*ring.Poly) {

	//----------------------------------------------------------------- A matrix
	pk, err := signer.LoadPublicKey("../public_key/public_key.json")
	if err != nil {
		log.Fatalf("public_key.json: %v", err)
	}
	A = [][]*ring.Poly{make([]*ring.Poly, len(pk.A))}
	for j, coeffs := range pk.A {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs) // already NTT in the fixture
		A[0][j] = p
	}

	//-------------------------------------------------------------- B-matrix set
	rawB, err := loadBmatrixCoeffs("../Parameters/Bmatrix.json")
	if err != nil {
		log.Fatalf("Bmatrix.json: %v", err)
	}
	toNTT := func(raw []uint64) *ring.Poly {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], raw)
		ringQ.NTT(p, p) // lift once
		return p
	}

	B0Const = []*ring.Poly{toNTT(rawB[0])}   // B₀,const
	B0Msg = [][]*ring.Poly{{toNTT(rawB[1])}} // one u-column
	B0Rnd = [][]*ring.Poly{{toNTT(rawB[2])}} // one x₀-column
	b1 = []*ring.Poly{toNTT(rawB[3])}        // b₁ row

	return
}
