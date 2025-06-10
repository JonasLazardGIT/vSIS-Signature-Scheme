### Detailed term-by-term decomposition of **A · x = u**

Below is the algebra that the Go code in **`Preimage_Sampler.GaussSamp`** realises.  
Indices:

|Symbol|Dimension|Where it is built in the code|
|---|---|---|
|`k`|gadget length κ = ⌈logₜ q⌉|`TrapGen` (★ line 72)|
|`A`|row vector of length k + 2|`TrapGen` (★ 119 – 136)|
|`x`|column vector of length k + 2|`GaussSamp` (★ 123 ff.)|
|`p`|perturbation, same length|`SamplePz` (return value)|
|`z`|column of length k (int64)|`SampleGDiscrete` (result)|

(★ numbers refer to the file-search snippets in the previous answer.)

---

#### 1. **Structure of A**

```
A  =  [ 1 ,  a ,  A₂[0] , … , A₂[k-1] ]     (length k+2)
```

- The first two polynomials are
    

```
A₀ = 1       (the constant “1” polynomial)
A₁ = a       (uniform from R_q)
```

- For every gadget index j ∈ {0,…, k-1}
    

```
A₂[j] = g_j  − (a ⋅ r̂_j  +  ê̂_j)
          └───────┬──────────┘   └──── “trapdoor” correction
           gadget column g_j
```

This definition (★ TrapGen 119-127) is checked at run-time by  
`✔ Trapdoor generation self-check passed`.

---

#### 2. **Structure of x**

`GaussSamp` assembles (★ 123-145):

```
x₀ = p₀  +   Σ_{j=0}^{k-1}  ê̂_j  ⋅ ẑ_j          (in NTT)
x₁ = p₁  +   Σ_{j=0}^{k-1}  r̂_j  ⋅ ẑ_j
x_{j+2} = p_{j+2} + ẑ_j             for j = 0 … k-1
```

Vector form

```
x = [ x₀ , x₁ , p₂+ẑ₀ , … , p_{k+1}+ẑ_{k-1} ]ᵀ .
```

---

#### 3. **Compute A · x and expand**

Write the product as two blocks

```
A·x = A₀·x₀  +  A₁·x₁  +  Σ_{j=0}^{k-1}  A₂[j] · x_{j+2}
```

Substitute the definitions of `x` and `A₂[j]`.

---

##### 3 a. Terms that contain **p**

```
A₀·p₀   +   A₁·p₁
+ Σ_{j}  A₂[j] · p_{j+2}
           = Σ_{j} (g_j − a r̂_j − ê̂_j) · p_{j+2}
             └────── g_j·p_{j+2}
             └────────────── −a r̂_j·p_{j+2}
             └──────────────── −ê̂_j·p_{j+2}
```

Collect every part that involves only the random perturbation vector **p**.  
By construction of **SamplePz** those polynomials already satisfy

```
A · p   ≡   u   −   G · z      (mod q)          (★ GaussSamp lines 83-96)
```

so after adding the **z**-part (next paragraph) the residual cancels.

---

##### 3 b. Terms that contain **z**

Plug `x_{j+2} = p_{j+2}+ẑ_j` and split the new _ẑ_ terms:

```
            z -terms coming from A₂
Σ_j  A₂[j] · ẑ_j
   = Σ_j (g_j − a r̂_j − ê̂_j) · ẑ_j
   =            Σ_j g_j·ẑ_j             (1)
     −  a · Σ_j r̂_j·ẑ_j                (2)
     −      Σ_j ê̂_j·ẑ_j                (3)
```

Now plug the _ẑ_ contributions from **x₀** and **x₁**:

```
A₀ · (Σ ê̂_j ẑ_j)  =  Σ ê̂_j ẑ_j                      (4)
A₁ · (Σ r̂_j ẑ_j)  =  a · Σ r̂_j ẑ_j                  (5)
```

Add (2) + (5) and (3) + (4):

```
a · Σ r̂_j ẑ_j   −  a · Σ r̂_j ẑ_j     = 0      (perfect cancellation)
Σ ê̂_j ẑ_j       −      Σ ê̂_j ẑ_j     = 0
```

All cross-terms involving the trapdoor vectors **r̂**, **ê̂** cancel exactly  
because they appear with opposite sign in `x₀`,`x₁` and in `A₂`.

What remains from the z-block is solely

```
Σ_j  g_j · ẑ_j   =   G · ẑ
```

where `G` is the gadget matrix `[g₀, …, g_{k-1}]`.

---

#### 4. **Put everything together**

```
A · x
 = (A · p)            +            (G · ẑ)
   └─────── equals u − G·ẑ
 = u .
```

The equality follows from the two design choices:

1. **SamplePz** generates `p` so that `A·p ≡ u − G·ẑ (mod q)`  
    (computed as `sub = u − A·p` in code, ★ 83-96).
    
2. **SampleGDiscrete** chooses `ẑ` to solve `G·ẑ ≡ sub (≡ u − A·p)`  
    (see call ★ 103-110).
    

Because both congruences are exact, all coefficients satisfy  
`(A·x)[t] = u[t] (mod q)` individually, which the Go program re-checks with

```
if AxCoeff.Coeffs[0][t] != uCoeffOrig.Coeffs[0][t] { log.Fatalf(...) }
```

(★ main output: `FAIL: verification mismatch at slot 0 …` when a bug exists).

---

### 5. Visual dependency chart

```
         r̂ , ê̂  ──────────────┐
                                │  used twice
                                │
                        ┌───────▼────────┐
                        │  SamplePz      │
             u  ──► A·p │  (Alg.4)       │
                        └───────┬────────┘
                                │ p₀,p₁,Q
                     sub = u − A·p
                                │
                                ▼
                        ┌───────┴────────┐
                        │ SampleGDiscrete│
                        │   (Alg.3)      │
                        └───────┬────────┘
                                │  ẑ
                                ▼
                     ┌──────────┴─────────────┐
                     │   assemble x           │
                     │ x₀,x₁ use r̂, ê̂ & ẑ   │
                     │ gadget rows  p+ẑ      │
                     └──────────┬─────────────┘
                                ▼
                           verify A·x=u
```


> **Notation legend used below**
> 
> - ϕ = ring dimension N · q = ring modulus
>    
> - NTT = evaluation (Montgomery/NTT) domain · COEFF = coefficient domain
>    
> - ⟨r̂, ẑ⟩ means the NTT dot-product
>
> - line numbers in brackets […] are from _Preimage_Sampling.go_ / _G_Sampling.go_  
>    (they are **hints only**; exact numbers vary with edits)
>    

---

## 0 · Overview of objects

|Symbol|Shape|Domain|Produced by|
|---|---|---|---|
|**A**|(k + 2) × 1 ring polys|NTT|**TrapGen**|
|**r̂**, **ê̂**|k ring polys each|NTT|**TrapGen**|
|**u**|1 ring poly|NTT (then COEFF)|caller / main|
|**p** = (p₀,p₁, q̂₀…q̂_{k-1})|k + 2 ring polys|NTT|**SamplePz**|
|**sub**|coeff vector length ϕ|COEFF|`sub ← InvNTT(u – A·p)`|
|**Z**|k × ϕ int64 matrix|—|**SampleGDiscrete**|
|**ẑ**|k ring polys|NTT|**ZtoZhat**|
|**x**|k + 2 ring polys|NTT|“assemble-x” step|

---

## 1 · Graph in Markdown (ASCII)

```
 ┌──────────────┐
 │   TrapGen    │
 │  (Alg. 1)    │
 └─────┬─┬─┬────┘
       │ │ │
       │ │ └── ê (k polys, NTT)
       │ └──── r̂ (k polys, NTT)
       └────── A  (row of k+2 polys, NTT)
                    ▲
                    │  (ring params, σ, base t)
 ┌──────────────────┼───────────────────────────────┐
 │                  │                               │
 │            ┌─────┴──────────────┐                │
 │            │    SamplePz        │                │
 │            │    (Alg. 4)        │                │
 │            └─────┬──────────────┘                │
 │                  │ p  =  [ p₀,p₁,q̂₀… ] (NTT)    │
 │                  ▼                               │
 │          sub = u – A·p   (COEFF)                 │
 │                  │                               │
 │            ┌─────┴──────────────┐                │
 │            │ SampleGDiscrete    │                │
 │            │      (Alg. 3)      │                │
 │            └─────┬──────────────┘                │
 │                  │  Z  (k×ϕ ints)                │
 │                  ▼                               │
 │            ┌──────────────┐                      │
 │            │   ZtoZhat    │                      │
 │            └─────┬────────┘                      │
 │                  │ ẑ (k polys, NTT)             │
 │                  ▼                               │
 │            ┌──────────────┐                      │
 │            │ assemble-x   │                      │
 │            └─────┬────────┘                      │
 │     x = [ p₀+⟨ê̂,ẑ⟩ , p₁+⟨r̂,ẑ⟩ , p₂+ẑ₀ … ]    │
 └──────────────────┼───────────────────────────────┘
                    ▼
              Verify  A·x ≟ u
```

---

## 2 · Edge catalogue with full detail

|#|Source ▶ Target|Data that flows|Representation in code|
|---|---|---|---|
|1|**TrapGen** → **SamplePz**|`A` (NTT)|returned slice `td.A` [TrapGen 119-136]|
|2|**TrapGen** → **SamplePz**|`r̂`, `ê̂` (NTT)|returned slice `td.R`|
|3|**u** (main) → **SamplePz**|target syndrome (NTT)|param `u` passed to `GaussSamp`|
|4|**SamplePz** → **p**|`p₀`,`p₁`, `q̂₀…` (NTT)|result of `SamplePz` [SamplePz 423-431]|
|5|(`A`,`p`,`u`) → **sub**|compute `sub = u − A·p`|lines 83-96 (NTT mult → InvNTT)|
|6|**sub** → **SampleGDiscrete**|centre-reduced coeffs `uCoeff`|lines 98-106|
|7|**SampleGDiscrete** → **Z**|k×ϕ integer matrix|function return|
|8|**Z** → **ZtoZhat**|CRT-encode & NTT|`ZtoZhat` 41-57|
|9|**ZtoZhat** → **ẑ**|k polys (NTT)|slice result|
|10|(`p`,`ẑ`,`r̂`,`ê̂`) → **assemble-x**|dot-products + sums|lines 123-145|
|11|**assemble-x** → **x**|signature vector (NTT)|return value of `GaussSamp`|
|12|(`A`,`x`,`u`) → **Verify**|check `A·x ≡ u` (COEFF)|debug block 162-170|

---

## 3 · Formulae carried on each arrow

| Edge | Payload expression                                                                    |
| ---- | ------------------------------------------------------------------------------------- |
| 1    | `A = [1, a, g₀ − (a r̂₀+ê̂₀), …, g_{k-1} − (a r̂_{k-1}+ê̂_{k-1})]`                    |
| 2    | Gaussian secrets `r̂_j, ê̂_j ∼ D_{σ_t}`                                               |
| 4    | `p = (p₀,p₁,q̂)`, where`(p₀,p₁)` solves 2 × 2 system via Alg 4, `q̂_j ∼ D_{√(s²−α²)}` |
| 5    | `sub = u − A·p (mod q)` then converted to centred ints                                |
| 6    | `Z[·,j]` satisfies gadget congruence                                                  |
| 7    | `∑_{i=0}^{k-1} Z[i,j] t^i ≡ sub[j] (mod q)`                                           |
| 8    | `ẑ_j = NTT( CRT(Z[j,:]) )`                                                           |
| 9    | `x₀ = p₀ + ⟨ê̂,ẑ⟩`, `x₁ = p₁ + ⟨r̂,ẑ⟩`, `x_{j+2} = p_{j+2}+ẑ_j`                    |
| 10   | final equality `A·x ≡ u (mod q)` checked coefficient-wise                             |

---

## 1 · Relations that are fully exercised in code

|#|Relation (informal)|Where / how it is tested|Status|
|---|---|---|---|
|T-1|**Trapdoor correctness** `g_j ≡ A₂[j] + a·r̂_j + ê̂_j`|`TrapGen` loop prints “✔ Trapdoor generation self-check passed” if every coefficient matches.|✓ tested|
|P-1|**Spectral-norm guard** `α‖[Tᵗ|I]‖₂ ≤ s`|Block in `SamplePz` aborts (`log.Panicf`) if the computed λ_max exceeds `s`.|
|P-2|**Rounding soundness** of the big-float centres `c₀ , c₁`|After applying `math.Round(-z·value)` the code checks that the centred integer really equals the float modulo q; panic on mismatch.|✓ tested|
|G-1|**Z‐matrix correctness** `∑ᵢ Z[i,j] t^i ≡ sub_j (mod q)`|In `SampleGDiscrete` a “recomb” loop recomputes the left-hand side and prints `diff`; it panics if non-zero.|✓ tested|
|X-1|**x₀ cancellation** `p₀ + ⟨ê̂,ẑ⟩ ≡ x₀ (mod q)`|“cancel-check x0 mismatch” loop over every slot; `log.Fatalf` on first failure.|✓ tested|
|X-2|**x₁ cancellation** `p₁ + ⟨r̂,ẑ⟩ ≡ x₁ (mod q)`|Same style loop (“cancel-check x1 mismatch”).|✓ tested|
|X-3|**Gadget rows** `A₂[j] + (a·r̂_j+ê̂_j) ≡ g_j` holds _after being multiplied by (p_{j+2}+ẑ_j)_|“G-recomb j=… diff=0” messages; any coefficient mismatch aborts.|✓ tested|
|V-1|**Whole pre-image** `A·x ≡ u (mod q)`|Final check in `main`: mismatch triggers fatal `FAIL: verification mismatch at slot …`.|✓ tested|

These steps together guarantee that:

- the trapdoor was generated consistently,
    
- all Gaussian samples respect their analytic bounds,
    
- the _arithmetic_ of `p`, `Z`, `x` is correct coefficient-by-coefficient, and
    
- the global equation `A·x=u` truly holds (or the program exits).
    

---

## 2 · Relations only _observed_ (printed) but **not enforced**

| #   | Relation / value                                                                 | What the code does                                                                      | Status     |
| --- | -------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- | ---------- |
| P-3 | **Empirical σ of q̂-samples** (`SamplePz → q̂ empirical standard deviation = …`) | Prints mean / σ; never compared to theory or bounded.                                   | ◑ printed  |
| P-4 | **`inner r·Q(0)` and `e·Q(0)`**                                                  | Diagnostic `fmt.Printf("SamplePz: inner r·Q(0)=…")`; no assertion.                      | ◑ printed  |
| P-5 | **`[CHECK] A·p coeff[0]`, `A·z coeff[0]`**                                       | Values shown; not required to be equal (they need to add to `u`, not match each other). | ◑ printed  |
| G-2 | **Karney variance per-digit**                                                    | No variance check inside `SampleGDiscrete`; only the recomposition test.                | ◑ implicit |

These printouts help a developer eyeball distributions but _do not stop execution_ if they drift.

---

## 3 · Logic that is _implicit_ (never tested)

| Component                                          | Potential hidden risk                                                                                                                                                        | Why it might matter |
| -------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------- |
| **Conversion path** `Z → CRT → NTT (ZtoZhat)`      | A swap of two CRT residues or a wrong root of unity would still pass G-1 recomposition (done before this step) but could break `x`; currently no direct check of `ẑ` ↔ `Z`. | ✗ untested          |
| **`SampleFZBig` / `Sample2zField` arithmetic**     | Those functions solve the 2×2 Gaussian system; only rounding sanity is checked, not the algebra `A'·(q₀,q₁) = (c₀,c₁) (mod q)`.                                              | ✗ untested          |
| **Domain flags** (`FieldElem.Domain = Eval/Coeff`) | A single forgotten `SwitchFormat()` would silently give wrong numbers yet pass earlier spectral checks.                                                                      | ✗ untested          |
| **Balanced-digit assumption** vs. unsigned digits  | Go uses unsigned base-t digits, PALISADE balanced digits; no automated test ensures both encodings lead to the _same_ `a`-values.                                            | ✗ untested          |
| **Randomness quality** of discrete Gaussians       | Empirical σ printed once, but no χ² or entropy test.                                                                                                                         | ✗ untested          |



---

### Summary

- All **critical equalities** that must hold for correctness of the scheme _are exercised_ in the current code path (section 1).
    
- A handful of **statistical or diagnostic quantities** are merely printed (section 2); deviations there would not abort.
    
- Several **conversion subtleties and corner-cases** remain unchecked (section 3).  
    They are safe under normal circumstances but represent the likely hiding spots for the “missing ±1” bug that triggered the final verification failure you observed.