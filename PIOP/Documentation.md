# Package PIOP – PACS‑aware Quadratic Gate Prototype

This document maps **every public symbol** in `quadratic_gate.go`,
`prover_helper.go`, and `PACS_Statement.go` to the *theoretical
objects* that appear in the SmallWood paper (version 2025‑1085).  Keep
this next to the code when you review or extend the prototype: you can
jump from a line of Go to the exact symbol in the paper and vice‑versa.

---

## 0. Bird’s‑eye view

```
smallwood.pdf      Go code (this repo)
──────────────     ─────────────────────────────────────
§5  PACS PIOP      BuildWitness , BuildGH , VerifyGH
§4  LVCS/PCS       pcs/  (separate package – imported here)
§3  DECS           decs/ (separate package – imported by pcs)
```

Our **PIOP layer** sits *on top of* the PCS commitment layer:

1. The prover starts with a *witness triple* `(w1,w2,w3)` that
   encodes a vSIS signature and auxiliary randomness.
2. Two families of constraints (parallel `f` and aggregated `f′`) are
   batched into polynomials **G** and **H**.
3. The verifier samples two random vectors `(γ,δ)` and accepts **iff**
   `G(X) ≡ 0` and `H(X) ≡ 0` in `F_q[X]`.

Everything below explains **how each Go helper realises one of those
three steps** and lists *unused helpers* that can be deleted later.

---

## 1. Mapping table (code ↔ paper)

| Paper symbol | Informal meaning                                     | Concrete object / function                                   |
| ------------ | ---------------------------------------------------- | ------------------------------------------------------------ |
| `w`          | witness matrix (n rows × s cols)                     | slices `w1 []*ring.Poly`, `w2 *ring.Poly`, `w3 []*ring.Poly` |
| `P_i(X)`     | blinded row polynomials (§6.1)                       | `BuildRowPolynomial`                                         |
| `M_i(X)`     | masking polynomials (§6.1, Eq.(3))                   | `BuildMaskPolynomials`                                       |
| `Θ_{j,i,k}`  | public constants inside *parallel* constraints       | *all zero* – produced by `BuildThetaZeros`                   |
| `Θ′_{j,i,k}` | public constants inside *aggregated* constraints     | interpolating polys from `BuildThetaPrimeSet`                |
| `f_j`        | parallel quadratic gate                              | hard‑coded as `w3_k - w1_k·w2` inside `BuildGH` (G‑part)     |
| `f′_j`       | aggregated linear equation (Eq.(1) in code comments) | implemented in `BuildGH` (H‑part)                            |
| `Γ, γ_i,j`   | LVCS challenge matrix & scalars                      | inside pcs/lvcs (opaque to PIOP layer)                       |
| `γ_i`        | Fiat–Shamir weights for parallel batch               | slice `gamma []uint64` drawn in `VerifyGHFromDisk`           |
| `δ_j`        | Fiat–Shamir weights for aggregated batch             | slice `delta []uint64` drawn in `VerifyGHFromDisk`           |
| `G(X)`       | batched poly for parallel constraints                | `BuildGH` returns it                                         |
| `H(X)`       | batched poly for aggregated constraints              | `BuildGH` returns it                                         |
| `Q_i(X)`     | Eq.(4) – **not implemented yet**                     | **TODO** (future work; today we batch into G/H instead)      |

---

## 2. Prover workflow (all in Go helpers)

1. **Witness construction** – `BuildWitness`

   * Verifies the “proof‑friendly” equation from the vSIS signature.
   * Produces `(w1,w2,w3)` where length `k = m+|u|+|x0|`.
2. **Blinding polynomials** (optional, Section 6 of the paper)

   * `BuildRowPolynomial` lifts each row of `w` to `P_i(X)` with ℓ
     random points – *perfect zero‑knowledge*.
   * `BuildMaskPolynomials` prepares `M_i(X)` so the Σ‐at‑Ω condition
     holds.
   * **NOTE**: the current quadratic‑gate demo *skips* this step and
     works directly with coefficient vectors; it is ready for upgrade.
3. **Constraint batching** – `BuildGH`

   * Computes `G = Σ_i γ_i·(w3_i - w1_i·w2)`  (parallel f)
   * Computes `H = Σ_j δ_j·rowErr_j`           (aggregated f′)
4. **Verification** – `VerifyGH`

   * Simply evaluates `G==0 && H==0` in NTT domain.

---

## 3. Which helpers are *not* used (can delete later)?

| File                | Function       | Rationale                              |
| ------------------- | -------------- | -------------------------------------- |
| `quadratic_gate.go` | `zeroPoly`     | trivial wrapper around `ring.NewPoly`  |
|                     | `randomVector` | only used in tests; keep if you like   |
| `PACS_Statement.go` | `BuildTheta`   | alias to `BuildThetaZeros` – redundant |

Delete them once the code stabilises; they add no semantic value.

---

## 4. Planned extensions

1. **Full PACS proof** – Implement `Q_i(X)` as in Equation (4) and add
   the Ω‑sum check (7); wire it to LVCS so only ℓ evaluations leak.
2. **Generic f / f′ registry** – Make `BuildGH` call a user‑supplied
   callback instead of hard‑coding the quadratic & linear forms.
3. **Replace toy data loaders** – remove JSON fixtures and plug real key
   generation + message hashing.
4. **`BuildQFromDisk` helper** – completed July 08.  Reconstructs all
   data from JSON, samples Γ′/γ′, calls `BuildQ`, and returns the slice
   `[]*ring.Poly` ready for commitment.

---
### What **`PACS_Simulation.go`** is and why it exists

`PACS_Simulation.go` is a *self-contained integration test* that walks through a **single, fully interactive execution** of the SmallWood PACS protocol on the toy “quadratic-gate” instance shipped with the repository.
It does **not** invent stand-in values; every public constant, witness component and ring parameter is pulled from the same JSON fixtures used by the real prover:

* `Parameters/Parameters.json` — ring degree *N* and modulus *q*
* `public_key/public_key.json` — the public matrix *A* already in NTT form
* `Parameters/Bmatrix.json` — the row-wise public constants for the aggregated constraints
* `Signature/Signature.json` — the message block *u*, masks x₀/x₁ and the lattice signature vector *s*

With those files it reproduces exactly the witness that the genuine verifier should see “under the hood”, then re-enacts the three protocol layers:

```
DECS  (degree-enforcing, small-domain commitment)
       ↓
LVCS  (linear-map commitment)
       ↓
PACS  (parallel + aggregated constraint system)
```

The code does not try to *reimplement* DECS or the Merkle tree—that is already tested elsewhere.
Instead it focuses on the parts that are **unique to PACS**:

1. building the parallel and aggregated constraint polynomials **F** and **F′**;
2. sampling the Fiat–Shamir matrices **Γ′** (polynomials) and **γ′** (scalars);
3. constructing each masking-row polynomial **Qᵢ(X)** exactly as in Equation (4);
4. re-checking

   * the point-wise consistency of (4) on a fresh evaluation point *e*, and
   * the Σ-at-Ω vanishing test of Equation (7);
5. collecting every message into an explicit, human-readable **transcript**.

If all tests pass the script prints:

```
Verifier ➜ ACCEPT – all checks passed.
```

Any failure is flagged as **REJECT** and the surrounding unit test (`TestPACSSimulation`) panics, so a continuous-integration run will turn red.

---

### High-altitude structure

```
┌────────────────────────────────────────┐
│ RunPACSSimulation()                    │
│  0  load ring parameters (N, q)        │
│  1  BuildWitnessFromDisk()             │
│  2  Ω := {1, …, s}                     │
│  3  buildFpar / buildFagg              │
│  4  sample Γ′, γ′  (verifier drives)   │
│  5  Q := BuildQ( … )                   │
│  6  verifier checks & pretty print     │
└────────────────────────────────────────┘
```

*Steps 0–1 — context & witness*
`BuildWitnessFromDisk` uses exactly the same helper pipeline as the genuine prover (Gauss sampling vectors, lifting to the NTT ring, sanity-checking the proof-friendly “row equation”) and returns the triplet **w₁, w₂, w₃** that will feed the PACS constraints.

*Step 2 — evaluation domain*
For the toy instance the number of witness columns *s* is small (≈ 32).
The simulation uses the first *s* integer roots Ω = {1,2,…,s} because they are guaranteed to be in the ring’s multiplicative subgroup and are exactly the points used when the witness rows were interpolated into polynomials.

*Step 3 — constraint polynomials*

* **Fpar**   one polynomial per column *k* realising the quadratic gate
  *Fparₖ(X) = w₃ₖ − w₁ₖ·w₂*

* **Fagg**  one polynomial per matrix row *j* realising the aggregated gate
  *F′ⱼ(X) = (b₁⊙A)s − (A·s)x₁ − B₀(1;u;x₀)*

Those two helpers are already part of the repository; the simulation calls them unmodified.

*Step 4 — Fiat–Shamir randomness*
Because we are debugging an *interactive* protocol the verifier, not the prover, draws the randomness:

* `sampleRandPolys` ⇒ **Γ′**<sub>i,j</sub>(X) for each masking row *i* and each parallel constraint *j* (degree ≤ s−1)
* `sampleRandMatrix` ⇒ **γ′**<sub>i,j</sub> field scalars for the aggregated batch

The degree of **Q** is fixed as *dQ = s + ℓ − 1* with ℓ = 1 (one random evaluation per row) so that Equation (3) of the paper is satisfied.

*Step 5 — building Q*
A single call to `BuildQ` does the heavy algebra:

```
Qᵢ(X) = Mᵢ(X)                              ← mask
       + Σ_j Γ′ᵢ,j(X) · Fpar_j(X)           ← parallel mix
       + Σ_j γ′ᵢ,j   · Fagg_j(X)            ← aggregated mix
```

*Step 6 — verifier checks*

1. **Point-wise consistency** of Equation (4) at `testPoint`.
   The code in `verifyRelationsOnE` converts everything back to coefficient form (one inverse-NTT per poly), evaluates both sides and bails if any row mismatches.
2. **Σ-at-Ω test** of Equation (7) via `VerifyQ`.
   It sums each `Qᵢ(ω)` for all ω ∈ Ω in the coefficient domain; a non-zero sum means the prover cheated.
3. **Merkle integrity** is marked as *true* here because the script is focused on algebra; the Merkle paths are unit-tested in `decs_merkle_test.go`.

A short `transcript` struct stores those booleans plus the important public values (Merkle root, first coefficients of Γ′, hash of the combined Rₖ polynomials, chosen opening set E).
`prettyPrintTranscript` dumps them in a neat order, so a developer can paste the output into an issue or a log.

---

### Key helpers the simulation re-uses

| helper                    | responsibility                              |
| ------------------------- | ------------------------------------------- |
| `BuildWitnessFromDisk`    | rebuild witness exactly like the prover     |
| `buildFpar` / `buildFagg` | implement the two constraint families       |
| `BuildMaskPolynomials`    | generate the *Mᵢ* rows satisfying ΣΩ Mᵢ = 0 |
| `BuildQ`                  | implement Equation (4) verbatim             |
| `VerifyQ`                 | implement Equation (7) verbatim             |

None of those helpers are changed by the simulation; you therefore exercise the very same paths the library will take in production.


