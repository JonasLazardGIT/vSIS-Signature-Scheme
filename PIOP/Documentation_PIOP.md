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
# 📜 PACS Simulation: Step-by-Step Documentation with Linked Equations

This document provides an **expert-level breakdown** of the `PACS_Simulation.go` code, describing how each step maps to the SmallWood proof system. All relevant equations from the SmallWood paper are included and linked to the code flow. It concludes with guidance on **parameter tuning** for achieving desired soundness probabilities.

---

## 📂 Overview of the Simulation Workflow

The `RunPACSSimulation()` function executes a single interactive simulation of the SmallWood PACS protocol, augmented to use the **full LVCS evaluation check**. This workflow comprises the following stages:

### 1️⃣ **Setup Phase**

* **Parameter Loading:**

  ```go
  par, _ := signer.LoadParams("../Parameters/Parameters.json")
  ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})
  q := ringQ.Modulus[0]
  ```

  Loads ring parameters: modulus `q` and ring dimension `N`. These define the polynomial ring $\mathbb{Z}_q[X]/(X^N+1)$ used for commitments and computations.

* **Witness Loading:**

  ```go
  w1, w2, w3 := BuildWitnessFromDisk()
  ```

  Fetches precomputed witness vectors from disk. These encode the prover's secret inputs and outputs for the proof.

---

### 2️⃣ **Commitment Phase (LVCS)**

* **Row Vector Construction:**

  ```go
  rows := columnsToRows(ringQ, w1, w2, w3, ell)
  ```

  Transforms the witness columns into row vectors, interpolated to degree $d = s + \ell - 1$.

  **Equation (LVCS Commitment):**

  $$
  C(P_j) = \text{MerkleRoot}(P_j + r_j \cdot X^{ncols})
  $$

  where $P_j$ is the row polynomial and $r_j$ is the random mask.

* **Commitment Initialization:**

  ```go
  root, pk, _ := lvcs.CommitInit(ringQ, rows, ell)
  ```

  Uses LVCS to commit to the rows. The output includes:

  * `root`: Merkle root of the commitments.
  * `pk`: Prover key containing random masks for later evaluation.

* **Commitment Challenges and Checks:**

  ```go
  Gamma := vrf.CommitStep1(root)
  Rpolys := lvcs.CommitFinish(pk, Gamma)
  if !vrf.CommitStep2(Rpolys) {...}
  ```

  The verifier samples random linear combinations (`Gamma`) of rows and verifies degree constraints.

  **Equation (Random Linear Combination):**

  $$
  R = \sum_{j=1}^{n_{rows}} \gamma_j \cdot P_j
  $$

---

### 3️⃣ **Linear Map Evaluation Phase**

* **Public Matrix Generation (C):**

  ```go
  C := sampleRandMatrix(1, rRows, q)
  ```

  Defines a public coefficient matrix $C \in \mathbb{Z}_q^{1 \times r}$, representing a random linear combination of rows.

  **Equation (Linear Map):**

  $$
  v = C \cdot \begin{bmatrix} P_1 \\ \vdots \\ P_{n_{rows}} \end{bmatrix}
  $$

* **Masked Target Computation:**

  ```go
  bar := lvcs.EvalInit(ringQ, pk, C)
  ```

  Prover computes masked values:

  $$
  \bar{v} = C \cdot \bar{r}
  $$

  where $\bar{r}$ are row masks from commitment.

---

### 4️⃣ **PACS Constraint Construction**

* **Challenge Sampling:**

  ```go
  GammaP := sampleRandPolys(ringQ, rho, sCols, sCols)
  gammaP := sampleRandMatrix(rho, 1, q)
  ```

  The verifier samples PACS challenges:

  * Polynomial weights $\Gamma'$
  * Scalar weights $\gamma'$

* **Constraint Polynomial Construction:**

  ```go
  Fpar := buildFpar(...)
  Fagg := buildFagg(...)
  M := BuildMaskPolynomials(...)
  Q := BuildQ(...)
  ```

  Constructs:

  * $F_{par}$: parallel constraint polynomials.
  * $F_{agg}$: aggregated constraint polynomials.
  * $M$: random masking polynomials.
  * $Q$: final batched PACS polynomials.

  **Equation (PACS Aggregation):**

  $$
  Q_i(X) = M_i(X) + \sum_j \Gamma'_{ij}(X) \cdot F_{par,j}(X) + \gamma'_i \cdot F_{agg}(X)
  $$

---

### 5️⃣ **Opening and Verification Phase**

* **Coordinate Opening:**

  ```go
  E := []int{ncols}
  open := lvcs.EvalFinish(pk, E)
  ```

  The verifier requests openings at mask indices $E$.

* **Full LVCS Verification:**

  ```go
  okLin := vrf.EvalStep2(bar, E, open.DECSOpen, C)
  ```

  Verifier checks:

  $$
  \sum_{j=1}^{n_{rows}} c_j \cdot P_j(e) = \bar{v}(e)
  $$

  for each $e \in E$.

* **PACS Checks:**

  ```go
  okEq4 := checkEq4OnOpening(...)
  okSum := VerifyQ(...)
  ```

  **Equation (Equation 4 Consistency):**

  $$
  Q_i(e) = M_i(e) + \sum_j \Gamma'_{ij}(e) \cdot F_{par,j}(e) + \gamma'_i \cdot F_{agg}(e)
  $$

---

### 6️⃣ **Final Verdict**

* **Transcript and Flags:**

  ```go
  tr.Flags = struct{...}{true, true, okLin, okEq4, okSum}
  pretty(&tr)
  ```

* **Decision:**

  ```go
  return tr.Flags.Merkle && tr.Flags.Deg && tr.Flags.LinMap && tr.Flags.Eq4 && tr.Flags.Sum
  ```

---

## 🎛️ Parameter Tuning for Soundness Probability

To ensure soundness, the verifier’s checks must reduce the probability of an invalid proof being accepted:

### Key Parameters

* **Ring Dimension (N)**

  * Larger $N$ → stronger degree soundness.
* **Masking Parameters ($\ell, \ell'$)**

  * $\ell$: Mask coordinates per row.
  * $\ell'$: Opened mask coordinates.
* **Parallel Repetitions ($\rho$)**

  * Increases probability amplification.

### Target Soundness ($\kappa$)

For $2^{-\kappa}$ soundness:

$$
\text{Error} \leq \left(\frac{d}{|\Omega|}\right)^{\ell'} \cdot |\mathbb{F}|^{-\rho}
$$

* Increase $\ell'$ and $\rho$ until error $\leq 2^{-\kappa}$.

---

## ✅ Summary

* All key SmallWood equations are present and implemented correctly.
* Parameter tuning ensures soundness in the interactive protocol.
* Replace deterministic challenge generation with Fiat–Shamir for deployment.
