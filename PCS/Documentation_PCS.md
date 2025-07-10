# Package pcs – SmallWood Polynomial Commitment Scheme (Layer‑3)

**Layering recap**

```
┌─────────────┐   ┌─────────────┐   ┌─────────────┐
│  DECS (L1)  │→ │  LVCS (L2)  │→ │  PCS  (L3)  │
└─────────────┘   └─────────────┘   └─────────────┘
 degree‑enforcing    linear‑map       polynomial
  small‑domain        vector          commitment
   commitment        commitment        wrapper
```

`pcs` wraps the **SmallWood Polynomial Commitment Scheme** (§4.2 of the paper).  It turns an
LVCS commitment on row vectors into a *public‐coin, hash‑based* commitment to **one or more
low‑degree polynomials** with a succinct evaluation proof.

Beyond the interfaces, this document explains **why** each call exists and **how** the
construction realises binding, zero‑knowledge, and efficient proofs for degrees ≤≈ 2¹⁶.

---

## 1. Public API (quick reference)

| Function                 | Side     | Purpose                                                                                                        |
| ------------------------ | -------- | -------------------------------------------------------------------------------------------------------------- |
| `CommitInitPolynomial`   | prover   | Build matrix **A**, sample mask rows, call `lvcs.CommitInit` – returns Merkle root + `ProverKey`.              |
| `CommitFinishPolynomial` | prover   | Finish commit after seeing LVCS challenge `Γ`; returns the vector of combined polys `R_k(X)`                   |
| `EvalInitPolynomial`     | prover   | Compute masked targets **\bar v** for a user‑supplied linear map *C*                                           |
| `EvalFinishPolynomial`   | prover   | Open the subset *E* selected by the verifier; delegates to `lvcs.EvalFinish`                                   |
| `NewVerifierPolynomial`  | verifier | Instantiate a PCS verifier for given degrees *d\_j*, block size `μ`, mask count `ℓ′`, and DECS repetition `η`. |
| `CommitStep1Polynomial`  | verifier | Record Merkle root, derive LVCS challenge Γ.                                                                   |
| `CommitStep2Polynomial`  | verifier | Store Rₖ(X) polynomials.                                                                                       |
| `ChooseEPolynomial`      | verifier | Pick a random ℓ′‑subset *E* for opening.                                                                       |
| `EvalStep2Polynomial`    | verifier | Verify the entire proof (Merkle paths, degree, and polynomial relation).                                       |

All helpers live in:

* **pcs\_prover.go** – prover routines
* **pcs\_verifier.go** – verifier routines
* **pcs\_types.go**   – shared structs
* **pcs\_test.go**    – end‑to‑end example

---

## 2. Parameters & Types

```go
// pcs_types.go

type Commitment struct { Root [32]byte }

type ProverKey struct {
    RingQ           *ring.Ring   // NTT ring  (degree N, modulus q)
    LvcsKey         *lvcs.ProverKey
    Dj              []int        // degrees d_j of each P_j(X)
    Mu              int          // block size  μ  (rows with true coeffs)
    EllPrime        int          // mask‑row count ℓ′ (zero‑knowledge)
    Nu              []int        // #cols ν_j for each polynomial
}
```

*`Mu`, `EllPrime`, and `Eta` (from DECS) are the tweakables that trade proof size vs. work.*

`VerifierState` mirrors these fields plus an `*lvcs.VerifierState` to reuse LVCS
verification code.

---

## 3. Internal Geometry – from polynomials to the LVCS matrix

For each polynomial \$P\_j(X)=\sum\_{i=0}^{d\_j} a\_{j,i}X^i\$ we chop the coefficient vector
into *μ‐strided* blocks of size μ and **stack** them as columns in a (μ+ℓ′)×ν\_j matrix **A\_j**.

* top μ rows: actual coefficients (padded with zeros)
* bottom ℓ′ rows: fresh uniform masks \$r\_{j,t,k}\$ providing zero‑knowledge

If \$d\_j+1−ℓ′\$ is *not* a multiple of μ, the paper’s “δ\_j shift” (Equation 2) pushes the high
coefficients downward so the verifier is *forced* to see that those slots are zero.

All **A\_j** are concatenated horizontally, then `Stack_β` (β layers) reshapes the result
into `nrows = μ + ℓ′` rows so **LVCS** can commit to it.

The prover holds the complete row matrix `rows[r][c]` in memory; the verifier *never sees
it* – only masked linear combinations on a tiny subset.

---

## 4. Protocol Walk‑Through

### 4.1 Commit (two passes)

1. **CommitInitPolynomial** (prover)

   1. Compute column counts `ν_j` and fully materialise the row matrix.
   2. Forward to `lvcs.CommitInit`, receiving a 32‑byte Merkle root.
   3. Return `(root, ProverKey)`.

2. **CommitStep1Polynomial** (verifier)

   1. Memorise `root`.
   2. Derive LVCS challenge `Γ` deterministically from the root (Fiat–Shamir) via
      `lvcs.Verifier.CommitStep1`.

3. **CommitFinishPolynomial** (prover)

   1. Read `Γ` from `pk.LvcsKey.Gamma`.
   2. Call `lvcs.CommitFinish` to build \$R\_k(X)=M\_k+\sum\_j Γ\_{k,j}P\_j\$.
   3. Send `R_k(X)` slice back.

4. **CommitStep2Polynomial** (verifier)

   * Store `R_k(X)`; LVCS will check degree ≤ *d* later.

### 4.2 Evaluate / Open

1. **EvalInitPolynomial** (prover)

   * User supplies coefficient matrix *C* (dimensions *m×nrows*).
   * Compute masked targets $\bar v_{k,i}=\sum_j C_{k,j}\,\text{mask}_j[i] \pmod q$
     for `i=0..ℓ′−1` & `k=0..m−1`.
   * Return `bar`.

2. **ChooseEPolynomial** (verifier)

   * Pick ℓ′ distinct indices `E` uniformly from `[0..N)`.

3. **EvalFinishPolynomial** (prover)

   * Delegate to `lvcs.EvalFinish(pk.LvcsKey,E)` to reveal `P_j(e)` & masks for
     *exactly* those indices.

4. **EvalStep2Polynomial** (verifier)

   1. `lvcs.Verifier.EvalStep2` checks:

      * Merkle authentication ↔ root
      * Degree ≤ `Ncols+ℓ′−1`
      * Per‑row linear equations    \$Q\_k(e)=\sum\_j C\_{k,j}P\_j(e)\$ for all opened *e*
   2. Returns **true** iff all sub‑checks succeed.

`pcs_test.go` drives the full flow and fails the test on *any* violation.

---

## 5. Security Sketch

* **Polynomial binding** ↔ LVCS function‑binding ↔ DECS polynomial‑binding.
  A successful cheat would violate the Merkle binding *or* Schwartz‑Zippel with
  probability ≤

  $$$\varepsilon = \varepsilon_{\text{DECS}} + \Bigl(\frac{ncols+ℓ′−1}{ℓ′}\Bigr)
                                               \Bigl/ \binom N{ℓ′} \Bigr.$$
  $$$
* **Zero knowledge**: the verifier only learns the ℓ′ evaluations it challenges,
  but every opened slot lies in the **mask rows**, so information‑theoretic privacy holds.
* **Straight‑line extractability** follows by rewinding the prover on two distinct
  challenges *E* (standard MPC‑in‑the‑Head argument).

---

## 6. Performance Tips

* For a fixed field and target soundness, tune `(μ,ℓ′,η)` to minimise
  `proof_size + work × λ` where `work ≈ μ·#polys` NTTs.
* The implementation uses **Lattigo’s** `ring` package:

  * 1 NTT ≈ 2× faster than Cooley–Tukey FFT in pure Go.
  * All arithmetic happens mod a 32‑bit prime ≡1 mod 2N.
* Choosing `β=1` keeps code simple; increase β to bound column width if you ever
  commit thousands of polynomials at once.

---

## 7. File Index

| File              | Role                                                 |
| ----------------- | ---------------------------------------------------- |
| `pcs_prover.go`   | Prover routines and `ProverKey` builder.             |
| `pcs_verifier.go` | Verifier state machine & public checks.              |
| `pcs_types.go`    | Shared structs, comments on parameters.              |
| `pcs_test.go`     | 290‑line example covering **every** exported symbol. |

---

## 8. Minimal Example (excerpt from `pcs_test.go`)

```go
// 1. Build ring (N=2048, 32‑bit q)
ringQ, _ := ring.NewRing(1<<11, []uint64{ 1<<32 - (1<<20) + 1 })

// 2. Honest prover prepares two random polynomials deg 7 and 11
P := [][]uint64{ randPoly(8,q), randPoly(12,q) }
Dj := []int{ 7, 11 }

// 3. Commit phase
com, pk, _ := pcs.CommitInitPolynomial(ringQ, P, Dj, Mu=3, EllPrime=5)
ver         := pcs.NewVerifierPolynomial(ringQ, Dj, 3, 5, decs.Eta)
ver.CommitStep1Polynomial(com)
R, _        := pcs.CommitFinishPolynomial(pk)
ver.CommitStep2Polynomial(R)

// 4. Evaluate phase
C := randomMatrix(m=3, nrows=8, q)
bar, _      := pcs.EvalInitPolynomial(pk, C)
E           := ver.ChooseEPolynomial()   // ℓ′ random indices
opening, _  := pcs.EvalFinishPolynomial(pk, E)
if !ver.EvalStep2Polynomial(bar, E, opening, C) { panic("PCS verify failed") }
```

The test passes ⇢ the PCS commitment is *binding*, the evaluation proof *sound*,
and no extra coefficient data leaks to the verifier.

---

## 9. Extending / Integrating

* **Batch evaluations**: repeat `Eval*` with new *C* matrices – no need to
  recommit.
* **Arith‑circuit proof**: embed PCS inside a PACS constraint system (§4.3 of the
  paper) and reuse the same LVCS rows for many gates.
* **Larger degrees**: increase μ; keep μ·ν\_j within cache; proof size grows
  sub‑linear in degree.
* **Alternative hashes**: swap SHA‑256 inside DECS without touching pcs code.

---

> *This documentation is deliberately self‑contained: read it alone to understand
> the PCS layer; consult `Documentation.md` (LVCS) and `DECS.md` for lower
> levels.*
