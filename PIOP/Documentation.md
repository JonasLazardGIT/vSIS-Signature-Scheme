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

*Last updated:* 2025‑07‑08
