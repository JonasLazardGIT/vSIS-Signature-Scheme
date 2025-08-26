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
	"fmt"
	"log"
	"math"
	"testing"
	"time"

	decs "vSIS-Signature/DECS"
	lvcs "vSIS-Signature/LVCS"
	signer "vSIS-Signature/Signer"
	prof "vSIS-Signature/prof"

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

type simCtx struct {
	ringQ  *ring.Ring
	q      uint64
	omega  []uint64
	w1     []*ring.Poly
	w2     *ring.Poly
	w3     []*ring.Poly
	A      [][]*ring.Poly
	b1     []*ring.Poly
	B0c    []*ring.Poly
	B0m    [][]*ring.Poly
	B0r    [][]*ring.Poly
	E      []int
	bar    [][]uint64
	Fpar   []*ring.Poly
	Fagg   []*ring.Poly
	M      []*ring.Poly
	Q      []*ring.Poly
	GammaP [][]uint64
	gammaP [][]uint64
	open   *lvcs.Opening
	vrf    *lvcs.VerifierState
	C      [][]uint64
	pk     *lvcs.ProverKey
}

// -------- Simulator options (for param sweeps & targeted faults) --------
type SimOpts struct {
	Ell      int  // # masked points per row (ℓ). If 0, use 1.
	Ncols    int  // |Ω|. If 0, use 8.
	Rho      int  // batching. If 0, use 1.
	FSforC   bool // derive C via FS. Default true in tests.
	TailOnly bool // force E on masked tail. Default true in tests.
	BindSqs  bool // bind Sqs in same commitment (default true).
}

func bumpAt(r *ring.Ring, p *ring.Poly, idx int, q uint64) {
	tmp := r.NewPoly()
	r.InvNTT(p, tmp)
	tmp.Coeffs[0][idx] = (tmp.Coeffs[0][idx] + 1) % q
	r.NTT(tmp, p)
}

func bumpConst(r *ring.Ring, p *ring.Poly, q uint64) { bumpAt(r, p, 0, q) }

func deepCopyOpen(o *decs.DECSOpening) *decs.DECSOpening {
	if o == nil {
		return nil
	}
	cp := *o
	if o.Indices != nil {
		cp.Indices = append([]int(nil), o.Indices...)
	}
	if o.Pvals != nil {
		cp.Pvals = make([][]uint64, len(o.Pvals))
		for i := range o.Pvals {
			cp.Pvals[i] = append([]uint64(nil), o.Pvals[i]...)
		}
	}
	if o.Paths != nil {
		cp.Paths = make([][][]byte, len(o.Paths))
		for i := range o.Paths {
			cp.Paths[i] = make([][]byte, len(o.Paths[i]))
			for j := range o.Paths[i] {
				cp.Paths[i][j] = append([]byte(nil), o.Paths[i][j]...)
			}
		}
	}
	return &cp
}

// --------------------------------------------------------------------------
// go test entry‑point
// --------------------------------------------------------------------------
func TestPACSSimulation(t *testing.T) {
	_, okLin, okEq4, okSum := buildSim(t)
	if !(okLin && okEq4 && okSum) {
		t.Fatalf("verifier rejected – some check failed")
	}
}

func TestPACSTampering(t *testing.T) {
	t.Run("LVCS/linear-map: tamper bar", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		ctx.bar[0][0] = (ctx.bar[0][0] + 1) % ctx.q
		if ctx.vrf.EvalStep2(ctx.bar, ctx.E, ctx.open.DECSOpen, ctx.C) {
			t.Fatalf("expected LVCS linear-map check to fail")
		}
	})
	t.Run("LVCS/tail-only: reject head index in E", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		Ehead := append([]int(nil), ctx.E...)
		Ehead[0] = 0
		openHead := lvcs.EvalFinish(ctx.pk, Ehead)
		if ctx.vrf.EvalStep2(ctx.bar, Ehead, openHead.DECSOpen, ctx.C) {
			t.Fatalf("expected EvalStep2 to reject head index")
		}
	})
	t.Run("LVCS/binding: mismatch between E and open.Indices", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		bad := deepCopyOpen(ctx.open.DECSOpen)
		if len(bad.Indices) > 0 {
			bad.Indices[0]++
		}
		if ctx.vrf.EvalStep2(ctx.bar, ctx.E, bad, ctx.C) {
			t.Fatalf("expected EvalStep2 to reject mismatched E vs open.Indices")
		}
	})
	t.Run("DECS/Merkle: tamper Pvals in opening", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		bad := deepCopyOpen(ctx.open.DECSOpen)
		bad.Pvals[0][0] = (bad.Pvals[0][0] + 1) % ctx.q
		if ctx.vrf.EvalStep2(ctx.bar, ctx.E, bad, ctx.C) {
			t.Fatalf("expected Merkle/masked check to fail")
		}
	})
	t.Run("DECS/Merkle: tamper Merkle path bytes", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		if len(ctx.open.DECSOpen.Paths) == 0 || len(ctx.open.DECSOpen.Paths[0]) == 0 {
			t.Skip("no path to tamper")
		}
		bad := deepCopyOpen(ctx.open.DECSOpen)
		bad.Paths[0][0][0] ^= 0x01
		if ctx.vrf.EvalStep2(ctx.bar, ctx.E, bad, ctx.C) {
			t.Fatalf("expected Merkle verification to fail")
		}
	})
	t.Run("Eq4: tamper Q", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		bumpConst(ctx.ringQ, ctx.Q[0], ctx.q)
		ok := checkEq4OnOpening(ctx.ringQ, ctx.Q, ctx.M, ctx.open, ctx.Fpar, ctx.Fagg, ctx.GammaP, ctx.gammaP, ctx.omega)
		if ok {
			t.Fatalf("expected Eq.(4) check to fail")
		}
	})
	t.Run("Eq4: tamper gammaPrime", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		ctx.gammaP[0][0] = (ctx.gammaP[0][0] + 1) % ctx.q
		ok := checkEq4OnOpening(ctx.ringQ, ctx.Q, ctx.M, ctx.open, ctx.Fpar, ctx.Fagg, ctx.GammaP, ctx.gammaP, ctx.omega)
		if ok {
			t.Fatalf("expected Eq.(4) check to fail")
		}
	})
	t.Run("SumΩ: tamper Q constant", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		bumpConst(ctx.ringQ, ctx.Q[0], ctx.q)
		if VerifyQ(ctx.ringQ, ctx.Q, ctx.omega) {
			t.Fatalf("expected ΣΩ check to fail")
		}
	})
	t.Run("LVCS/degree: tamper R_k degree", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		if ctx.ringQ.N-1 <= decs.Degree {
			t.Skip("ring dimension too small to exceed degree bound")
		}
		R := make([]*ring.Poly, len(ctx.vrf.R))
		for i, p := range ctx.vrf.R {
			R[i] = p.CopyNew()
		}
		coeff := ctx.ringQ.NewPoly()
		ctx.ringQ.InvNTT(R[0], coeff)
		idx := decs.Degree + 1
		coeff.Coeffs[0][idx] = (coeff.Coeffs[0][idx] + 1) % ctx.q
		ctx.ringQ.NTT(coeff, R[0])
		if ctx.vrf.CommitStep2(R) {
			t.Fatalf("expected degree bound to fail")
		}
	})
	t.Run("FullPACS: w3 != w1*w2", func(t *testing.T) {
		t.Parallel()
		ctx, _, _, _ := buildSim(t)
		bumpConst(ctx.ringQ, ctx.w3[0], ctx.q)
		if VerifyFullPACS(ctx.ringQ, ctx.w1, ctx.w2, ctx.w3, ctx.A, ctx.b1, ctx.B0c, ctx.B0m, ctx.B0r) {
			t.Fatalf("expected FullPACS to fail on w3≠w1·w2")
		}
	})
}

func TestPACSParamGrid(t *testing.T) {
	cases := []SimOpts{
		{Ell: 1, Ncols: 8, Rho: 1, FSforC: true, TailOnly: true, BindSqs: true},
	}
	for i, o := range cases {
		o := o
		t.Run(fmt.Sprintf("case-%d", i), func(t *testing.T) {
			t.Parallel()
			_, okLin, okEq4, okSum := buildSimWith(t, o)
			if !(okLin && okEq4 && okSum) {
				t.Fatalf("verifier rejected for opts %+v", o)
			}
		})
	}
}

func TestEq4TamperMaskOnly(t *testing.T) {
	ctx, _, _, _ := buildSim(t)
	bumpConst(ctx.ringQ, ctx.M[0], ctx.q)
	if checkEq4OnOpening(ctx.ringQ, ctx.Q, ctx.M, ctx.open, ctx.Fpar, ctx.Fagg, ctx.GammaP, ctx.gammaP, ctx.omega) {
		t.Fatalf("Eq.(4) should fail when M is tampered")
	}
}

func TestOmegaRejectsDuplicates(t *testing.T) {
	par, _ := signer.LoadParams("../Parameters/Parameters.json")
	ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})
	q := ringQ.Modulus[0]
	omega := []uint64{1 % q, 2 % q, 1 % q}
	if err := checkOmega(omega, q); err == nil {
		t.Fatalf("checkOmega must reject duplicates")
	}
}

func TestPACSDeterminism(t *testing.T) {
	_, a1, b1, c1 := buildSim(t)
	_, a2, b2, c2 := buildSim(t)
	if a1 != a2 || b1 != b2 || c1 != c2 {
		t.Fatalf("verdicts changed")
	}
}

func BenchmarkBuildSim(b *testing.B) {
	for i := 0; i < b.N; i++ {
		buildSimWith(nil, SimOpts{Ell: 1, Ncols: 8, Rho: 1, FSforC: true, TailOnly: true, BindSqs: true})
	}
}

func buildSim(t *testing.T) (*simCtx, bool, bool, bool) {
	return buildSimWith(t, SimOpts{Ell: 1, Ncols: 8, Rho: 1, FSforC: true, TailOnly: true, BindSqs: true})
}

func buildSimWith(t *testing.T, o SimOpts) (*simCtx, bool, bool, bool) {
	defer prof.Track(time.Now(), "buildSimWith")
	if o.Ell == 0 {
		o.Ell = 1
	}
	if o.Ncols == 0 {
		o.Ncols = 8
	}
	if o.Rho == 0 {
		o.Rho = 1
	}
	// ------------------------------------------------------------- parameters
	par, _ := signer.LoadParams("../Parameters/Parameters.json")
	ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})
	q := ringQ.Modulus[0]

	// ------------------------------------------------------------- witnesses
	w1, w2, w3 := BuildWitnessFromDisk() // helper in another PIOP file
	A, b1, B0c, B0m, B0r := loadPublicTables(ringQ)

	// --- NEW: remake signature rows as coefficient-packing rows over Ω -------------
	ell := o.Ell
	// build the evaluation grid Ω (size ncols used below)
	ncols := o.Ncols
	px := ringQ.NewPoly()
	px.Coeffs[0][1] = 1
	pts := ringQ.NewPoly()
	ringQ.NTT(px, pts)
	omega := pts.Coeffs[0][:ncols]
	if err := checkOmega(omega, q); err != nil {
		fmt.Println("[Ω-check] ", err)
		return nil, false, false, false
	}
	S0 := uint64(len(omega))
	S0inv := modInv(S0, q)

	// length of signature block
	mSig := len(w1) - len(B0m) - len(B0r)

	// Rebuild top mSig rows: P_t(ω_j) = a_{t,j} (coefficient packing + blinding)
	for t := 0; t < mSig; t++ {
		coeff := ringQ.NewPoly()
		ringQ.InvNTT(w1[t], coeff) // coefficient vector of the ring poly
		vals := make([]uint64, len(omega))
		for j := 0; j < len(omega); j++ {
			// **Coefficient packing**: per-column value is the coefficient a_{t,j}
			vals[j] = coeff.Coeffs[0][j] % q
		}
		w1[t] = buildValueRow(ringQ, vals, omega, ell) // deg ≤ s+ell-1 row poly
	}

	// Rebuild message and x0 rows as **column-constant** packing rows.
	for i := 0; i < len(B0m); i++ {
		tmp := ringQ.NewPoly()
		ringQ.InvNTT(w1[mSig+i], tmp)
		c := tmp.Coeffs[0][0] % q
		vals := make([]uint64, len(omega))
		for j := range vals {
			vals[j] = c
		}
		w1[mSig+i] = buildValueRow(ringQ, vals, omega, ell)
	}
	off := mSig + len(B0m)
	for i := 0; i < len(B0r); i++ {
		tmp := ringQ.NewPoly()
		ringQ.InvNTT(w1[off+i], tmp)
		c := tmp.Coeffs[0][0] % q
		vals := make([]uint64, len(omega))
		for j := range vals {
			vals[j] = c
		}
		w1[off+i] = buildValueRow(ringQ, vals, omega, ell)
	}

	// Recompute w3 = w1 * w2 using the updated packing rows.
	for i := 0; i < len(w1); i++ {
		ringQ.MulCoeffs(w1[i], w2, w3[i])
	}

	// β and radix R=2^w
	beta := par.Beta
	if beta*beta >= q {
		beta = uint64(math.Sqrt(float64(q - 1)))
	}
	wbits := 12
	spec := NewBoundSpec(q, beta, wbits, ell)

	// Allocate witness columns for decomposition
	cols := appendDecompositionColumns(ringQ, spec.LS, spec.W)

	// carry width Wc = ceil(log2(|Ω|)) + 1
	Wc := 1
	for (1 << (Wc - 1)) < len(omega) {
		Wc++
	}
	glob := appendGlobalCarrys(ringQ, spec.LS, Wc)
	slack := appendGlobalSlack(ringQ, spec.LS, spec.W)

	// Append all new columns to w1 so they are committed
	origW1Len := len(w1)
	for _, p := range cols.D {
		w1 = append(w1, p)
	}
	for _, p := range cols.T {
		w1 = append(w1, p)
	}
	for l := 0; l < spec.LS; l++ {
		for _, b := range cols.Bit[l] {
			w1 = append(w1, b)
		}
	}
	for _, p := range glob.C {
		w1 = append(w1, p)
	}
	for l := 0; l < len(glob.CBits); l++ {
		for _, b := range glob.CBits[l] {
			w1 = append(w1, b)
		}
	}
	for _, p := range slack.D {
		w1 = append(w1, p)
	}
	for l := 0; l < len(slack.DBits); l++ {
		for _, b := range slack.DBits[l] {
			w1 = append(w1, b)
		}
	}

	// Fill all new columns via interpolation over Ω
	Sqs, err := ProverFillIntegerL2(ringQ, w1, mSig, spec, cols, glob, slack, omega, ell, S0, S0inv)
	if err != nil {
		panic(err)
	}
	// Bind Sqs in the same commitment as T,D,Bit: needed for T0 - Sqs = 0.
	if o.BindSqs {
		w1 = append(w1, Sqs)
	}

	// (tight carry width): U_C = ceil(log2 |Ω|) + 1
	{
		s := len(omega)
		Wc = 1
		for (1 << uint(Wc-1)) < s {
			Wc++
		}
	}

	// ---------------------------------------------------------- LVCS.Commit
	rows := columnsToRows(ringQ, w1, w2, w3, ell, omega)
	root, pk, _ := lvcs.CommitInit(ringQ, rows, ell)

	vrf := lvcs.NewVerifier(ringQ, len(rows), decs.Eta, ncols)
	Gamma := vrf.CommitStep1(root)
	Rpolys := lvcs.CommitFinish(pk, Gamma)
	if !vrf.CommitStep2(Rpolys) {
		fmt.Println("[deg‑chk] R failed")
		return nil, false, false, false
	}

	// ---------------------------------------------------------- coefficient C
	// We prove that  v = C·rows  where C is a *1×r* public matrix.
	rRows := len(rows)
	var C [][]uint64
	if o.FSforC {
		fsC := newFSRNG("C", root[:])
		C = sampleFSMatrix(1, rRows, q, fsC)
	} else {
		C = sampleRandMatrix(1, rRows, q)
	}

	// ------------------------------------------------------- compute  bar[k][i]
	bar := lvcs.EvalInit(ringQ, pk, C)

	// ------------------------------------------------------- PACS batching
	rho := o.Rho

	FparCore := buildFpar(ringQ, w1[:origW1Len], w2, w3)
	FparDec := buildFparIntegerDecomp(ringQ, Sqs, spec, cols)
	FparCouple := BuildFparSqsCoupling(ringQ, Sqs, w1[:mSig])
	FparCarr := buildFparGlobCarryBits(ringQ, S0inv, Wc, glob)
	FparSlackB := buildFparGlobSlackBits(ringQ, S0inv, spec.W, slack)
	Fpar := append(FparCore, FparDec...)
	Fpar = append(Fpar, FparCouple...)
	Fpar = append(Fpar, FparCarr...)
	Fpar = append(Fpar, FparSlackB...)

	theta := BuildThetaPrimeSet(ringQ, A, b1, B0c, B0m, B0r, omega)
	FaggBBS := buildFaggOnOmega(ringQ, w1[:origW1Len], w2, theta, mSig)
	FaggInt := buildFaggIntegerSumDelta(ringQ, spec, S0, S0inv, cols, glob, slack, omega)
	Fagg := append(FaggBBS, FaggInt...)

	totalParallel := len(Fpar)
	totalAgg := len(Fagg)
	// Fiat–Shamir bind Γ′, γ′ to transcript material:
	// root || hash(R) || C || Ω
	Rh := sha256.Sum256(polysToBytes(Rpolys))
	fsMat := bytesU64Mat(C)
	fsOmg := bytesU64Vec(omega)
	fsGamma := newFSRNG("GammaPrime", root[:], Rh[:], fsMat, fsOmg)
	fsGammaSc := newFSRNG("gammaPrime", root[:], Rh[:], fsMat, fsOmg)
	GammaP := sampleFSMatrix(rho, totalParallel, q, fsGamma)
	gammaP := sampleFSMatrix(rho, totalAgg, q, fsGammaSc)

	fmt.Printf("→ parallel rows: %d; aggregated rows: %d; witness cols: %d\n", totalParallel, totalAgg, len(w1))

	dQ := len(w1) + ell - 1
	sumFpar := sumPolyList(ringQ, Fpar, omega)
	sumFagg := sumPolyList(ringQ, Fagg, omega)
	M := BuildMaskPolynomials(ringQ, rho, dQ, omega, GammaP, gammaP, sumFpar, sumFagg)
	Q := BuildQ(ringQ, M, Fpar, Fagg, GammaP, gammaP)

	// --------------------------------------------------------- verifier picks E
	var E []int
	if o.TailOnly {
		if ell > 1 {
			E = make([]int, ell)
			for i := 0; i < ell; i++ {
				E[i] = ncols + i
			}
		} else {
			E = []int{ncols}
		}
	} else {
		E = vrf.ChooseE(ell, ncols)
	}

	// --------------------------------------------------------- opening & check
	open := lvcs.EvalFinish(pk, E)

	ctx := &simCtx{
		ringQ:  ringQ,
		q:      q,
		omega:  omega,
		w1:     w1,
		w2:     w2,
		w3:     w3,
		A:      A,
		b1:     b1,
		B0c:    B0c,
		B0m:    B0m,
		B0r:    B0r,
		E:      E,
		bar:    bar,
		Fpar:   Fpar,
		Fagg:   Fagg,
		M:      M,
		Q:      Q,
		GammaP: GammaP,
		gammaP: gammaP,
		open:   open,
		vrf:    vrf,
		C:      C,
		pk:     pk,
	}

	okLin := vrf.EvalStep2(bar, E, open.DECSOpen, C)
	okEq4 := checkEq4OnOpening(ringQ, Q, M, open, Fpar, Fagg, GammaP, gammaP, omega)
	okSum := VerifyQ(ringQ, Q, omega)

	return ctx, okLin, okEq4, okSum
}

// ============================================================================
// Helpers
// ============================================================================

func columnsToRows(r *ring.Ring, w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly, ell int, omega []uint64) [][]uint64 {
	defer prof.Track(time.Now(), "columnsToRows")
	s := len(w1)
	ncols := len(omega)
	rows := make([][]uint64, s+2)
	q := r.Modulus[0]

	// Row 0..s-1: for each witness column k, evaluate w1[k](ω_j) for all j.
	tmp := r.NewPoly()
	for k := 0; k < s; k++ {
		rows[k] = make([]uint64, ncols)
		r.InvNTT(w1[k], tmp) // coeff domain of w1[k]
		for j := 0; j < ncols; j++ {
			rows[k][j] = EvalPoly(tmp.Coeffs[0], omega[j]%q, q)
		}
	}

	// Row s: w2(ω_j)
	r.InvNTT(w2, tmp)
	rows[s] = make([]uint64, ncols)
	for j := 0; j < ncols; j++ {
		rows[s][j] = EvalPoly(tmp.Coeffs[0], omega[j]%q, q)
	}

	// Row s+1: per‑column product w3[col](ω_col) on the diagonal; 0 elsewhere.
	rows[s+1] = make([]uint64, ncols)
	for col := 0; col < ncols && col < len(w3); col++ {
		r.InvNTT(w3[col], tmp)
		rows[s+1][col] = EvalPoly(tmp.Coeffs[0], omega[col]%q, q)
	}

	return rows
}

// Eq.(4) consistency on each opened index (unchanged)
func checkEq4OnOpening(r *ring.Ring, Q, M []*ring.Poly, op *lvcs.Opening,
	Fpar []*ring.Poly, Fagg []*ring.Poly, GammaP [][]uint64, gammaP [][]uint64, omega []uint64) bool {
	defer prof.Track(time.Now(), "checkEq4OnOpening")

	q := r.Modulus[0]
	tmp := r.NewPoly()
	for i, idx := range op.DECSOpen.Indices {
		j := idx - 1
		w := omega[j]
		r.InvNTT(Q[i], tmp)
		lhs := EvalPoly(tmp.Coeffs[0], w, q)

		rhs := evalAt(r, M[i], w)
		for t := range GammaP[i] {
			g := GammaP[i][t]
			f := evalAt(r, Fpar[t], w)
			rhs = modAdd(rhs, modMul(g, f, q), q)
		}
		for t := 0; t < len(Fagg); t++ {
			g := gammaP[i][t]
			f := evalAt(r, Fagg[t], w)
			rhs = modAdd(rhs, modMul(g, f, q), q)
		}

		if lhs != rhs {
			return false
		}
	}
	return true
}

// ------- small utilities  ---------------------------------------
func firstColumn(mat [][]uint64) [][]uint64 {
	out := make([][]uint64, len(mat))
	for i := range mat {
		if len(mat[i]) > 0 {
			out[i] = []uint64{mat[i][0]}
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
	defer prof.Track(time.Now(), "loadPublicTables")

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
