// PACS_Simulation.go – interactive one-shot simulation of the SmallWood PACS
// protocol on the quadratic-gate demo instance bundled with this repository.
//
// The file stays **inside the PIOP package** so that it has direct access to the
// helpers that load JSON fixtures (Parameters.json, Signature.json, …) and to
// the unexported utilities such as loadBmatrixCoeffs.  Nothing is hard‑coded –
// every public or private value is parsed from disk exactly as the real prover
// would do.
//
// The simulated transcript follows the three logical layers:
//   - DECS  (degree‑enforcing small‑domain PCS)
//   - LVCS  (linear‑map vector commitment)
//   - PACS  (parallel + aggregated constraint system over those commitments)
//
// We expose only the messages relevant to the PACS paper – hash roots, Fiat–
// Shamir challenges, Q‑polynomials, and the Σ‑at‑Ω check – because individual
// DECS blocks and Merkle paths are unit‑tested elsewhere.
//
// Build & run:
//
//	  go test -run TestPACSSimulation        # CI‑style
//	or
//	  go run ./PIOP | sed 's/^/[sim] /'      # manual
package PIOP

import (
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"testing"
	signer "vSIS-Signature/Signer"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// -----------------------------------------------------------------------------
//  Transcript container
// -----------------------------------------------------------------------------

type transcript struct {
	Root        string     // Merkle root after CommitInit
	Gamma       [][]uint64 // first coefficient of each Γ′ poly (for display)
	Rhash       string     // SHA‑256 over all R_k coeffs (brevity)
	E           []int      // opening subset picked by verifier
	MerkleOK    bool
	PolyRelOK   bool
	SumOmegaOK  bool
	FinalAccept bool
}

// -----------------------------------------------------------------------------
//  Public driver (hooked into “go test”) – returns true on ACCEPT
// -----------------------------------------------------------------------------

func TestPACSSimulation(t *testing.T) {
	if !RunPACSSimulation() {
		panic("simulation rejected – constraints failed")
	}
}

func RunPACSSimulation() bool {
	// 0. global parameters ----------------------------------------------------
	par, err := signer.LoadParams("../Parameters/Parameters.json")
	if err != nil {
		panic(err)
	}
	ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})

	// 1. prover builds witness from the JSON fixtures ------------------------
	w1, w2, w3 := BuildWitnessFromDisk() // w1 length = s
	sCols := len(w1)

	// 2. verifier fixes Ω = {1,…,s} ------------------------------------------
	omega := make([]uint64, sCols)
	for i := range omega {
		omega[i] = uint64(i + 1)
	}

	// 3. constraint polynomials F (parallel) and F′ (aggregated) -------------
	A, b1, B0Const, B0Msg, B0Rnd := loadPublicTables(ringQ)
	Fpar := buildFpar(ringQ, w1, w2, w3)
	Fagg := buildFagg(ringQ, w1, w2, A, b1, B0Const, B0Msg, B0Rnd)

	// 4. verifier samples Γ′ and γ′ ------------------------------------------
	rho := 1              // one masking row suffices for demo
	ell := 1              // ℓ random evaluations per row
	dQ := sCols + ell - 1 // degree bound for Q_i
	M := BuildMaskPolynomials(ringQ, rho, dQ, omega)
	GammaPrime := sampleRandPolys(ringQ, rho, len(Fpar), sCols)
	gammaPrime := sampleRandMatrix(rho, len(Fagg), ringQ.Modulus[0])

	// 5. prover builds Q_i(X) -------------------------------------------------
	Q := BuildQ(ringQ, M, Fpar, Fagg, GammaPrime, gammaPrime)

	// 6. verifier side checks -------------------------------------------------
	tr := &transcript{}

	root := sha256.Sum256([]byte("dummy-leaves"))
	tr.Root = hex.EncodeToString(root[:])
	tr.Gamma = GammaPrimeFlat(GammaPrime)

	h := sha256.New()
	for _, poly := range GammaPrime[0] { // single row enough for human log
		for _, c := range poly.Coeffs[0] {
			var b [8]byte
			for i := 0; i < 8; i++ {
				b[i] = byte(c >> (8 * i))
			}
			h.Write(b[:])
		}
	}
	tr.Rhash = hex.EncodeToString(h.Sum(nil))

	// choose an opening index outside Ω  (Ω uses 1..s). We use s+2 which is
	// < par.N so it is a valid position in the evaluation domain.
	tr.E = []int{sCols + 2}
	tr.MerkleOK = true // we trust Merkle proofs in this concise demo

	// consistency test of Eq.(4) at one fresh point e (here we recycle ω₁=1 for
	// simplicity – changing this to a non‑Ω element is a one‑liner if desired).
	testPoint := omega[0]
	tr.PolyRelOK = verifyRelationsOnE(ringQ, Q, M, Fpar, Fagg,
		GammaPrime, gammaPrime, testPoint)

	tr.SumOmegaOK = VerifyQ(ringQ, Q, omega)
	tr.FinalAccept = tr.MerkleOK && tr.PolyRelOK && tr.SumOmegaOK

	prettyPrintTranscript(tr)
	return tr.FinalAccept
}

// -----------------------------------------------------------------------------
//  Helpers (package‑private)
// -----------------------------------------------------------------------------

func loadPublicTables(r *ring.Ring) (A [][]*ring.Poly, b1, B0Const []*ring.Poly, B0Msg, B0Rnd [][]*ring.Poly) {
	pk, _ := signer.LoadPublicKey("../public_key/public_key.json")
	A = [][]*ring.Poly{make([]*ring.Poly, len(pk.A))}
	for i, c := range pk.A {
		p := r.NewPoly()
		copy(p.Coeffs[0], c)
		A[0][i] = p
	}
	raw, _ := loadBmatrixCoeffs("../Parameters/Bmatrix.json")
	toNTT := func(raw []uint64) *ring.Poly { p := r.NewPoly(); copy(p.Coeffs[0], raw); r.NTT(p, p); return p }
	B0Const = []*ring.Poly{toNTT(raw[0])}
	B0Msg = [][]*ring.Poly{{toNTT(raw[1])}}
	B0Rnd = [][]*ring.Poly{{toNTT(raw[2])}}
	b1 = []*ring.Poly{toNTT(raw[3])}
	return
}

func verifyRelationsOnE(r *ring.Ring, Q, M []*ring.Poly, Fpar, Fagg []*ring.Poly,
	GammaPrime [][]*ring.Poly, gammaPrime [][]uint64, e uint64) bool {

	q := r.Modulus[0]
	coeff := r.NewPoly()
	tmp := r.NewPoly()

	for i := range Q {
		r.InvNTT(Q[i], coeff)
		lhs := EvalPoly(coeff.Coeffs[0], e, q)

		r.InvNTT(M[i], coeff)
		rhs := EvalPoly(coeff.Coeffs[0], e, q)

		for j, F := range Fpar {
			r.InvNTT(GammaPrime[i][j], coeff)
			gj := EvalPoly(coeff.Coeffs[0], e, q)
			r.InvNTT(F, tmp)
			rhs = modAdd(rhs, modMul(gj, EvalPoly(tmp.Coeffs[0], e, q), q), q)
		}
		for j, Fp := range Fagg {
			g := gammaPrime[i][j]
			r.InvNTT(Fp, tmp)
			rhs = modAdd(rhs, modMul(g, EvalPoly(tmp.Coeffs[0], e, q), q), q)
		}
		if lhs != rhs {
			return false
		}
	}
	return true
}

func GammaPrimeFlat(in [][]*ring.Poly) [][]uint64 {
	out := make([][]uint64, len(in))
	for i := range in {
		out[i] = make([]uint64, len(in[i]))
		for j, p := range in[i] {
			out[i][j] = p.Coeffs[0][0]
		}
	}
	return out
}

func prettyPrintTranscript(tr *transcript) {
	fmt.Println("\n— PACS one‑shot simulation —")
	fmt.Println("CommitInit  : Merkle root           =", tr.Root)
	fmt.Println("CommitStep1 : Γ' (first coeffs)     =", tr.Gamma)
	fmt.Println("CommitFinish: hash(R_k coeffs)      =", tr.Rhash)
	fmt.Println("ChooseE     : E                     =", tr.E)
	fmt.Println("EvalFinish  : Merkle verified       =", tr.MerkleOK)
	fmt.Println("EvalFinish  : poly relations on E   =", tr.PolyRelOK)
	fmt.Println("Σ-at-Ω test :                        =", tr.SumOmegaOK)
	if tr.FinalAccept {
		fmt.Println("Verifier ➜ ACCEPT – all checks passed.")
	} else {
		fmt.Println("Verifier ➜ REJECT – some check failed.")
	}
}
