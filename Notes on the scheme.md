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

Below every symbol that appears in the derivation is pinned down by an explicit **equation** or **constraint** showing how it is coupled to the target *u*, the public matrix *A*, the gadget *G*, or the trapdoor $(\hat r,\hat e)$.

---

### 1 Global setting

*Ring*. $R_q = \mathbb{Z}_q[X]/(X^N+1)$ with modulus $q$ and dimension $N$.
*Gadget*. Pick a radix $t\;(=2\text{ or }3\text{ most often})$ and set

$$
G \;=\;[\,g_0,\dots,g_{k-1}\,],\quad
g_j \;=\;t^{\,j}\;\in R_q,
\qquad k=\lceil\log_t q\rceil .
$$

---

### 2 Trapdoor and public key

The trapdoor is the ordered pair

$$
T \;=\;(\hat r,\hat e), \qquad
\hat r,\hat e \in R_q^{k}.
$$

Choose a fresh $$a\overset{\$}{\leftarrow} R_q$$ and publish the **1 × (k+2)** vector

$$
A \;=\;[1,\;a,\;A_2[0],\dots,A_2[k\!-\!1]],
$$

$$
\boxed{A_2[j] = g_j - (a\hat r_j+\hat e_j)}\qquad(j=0,\dots,k-1).
$$

Because $A_2$ is built from $G$ and $T$, anyone who knows $T$ can later cancel those two hidden terms.

---

### 3 Target instance

The signer/verifier supplies the target polynomial

$$
u \in R_q.
$$

---

### 4 Joint perturbation–gadget solution

The sampler first draws a discrete-Gaussian vector

$$
\hat z = (\hat z_0,\dots,\hat z_{k-1})\in R_q^{k},
$$

then finds a small-norm **perturbation**

$$
p=(p_0,p_1,p_2,\dots,p_{k+1})\in R_q^{k+2}
$$

by solving the linear congruence

$$
\boxed{A\cdot p \;=\; u - G\cdot\hat z \pmod{q}}\qquad(\star)
$$

with the additional requirement that the coefficients of $p$ lie in a narrow Gaussian; the concrete SearchPz routine fulfils both goals simultaneously.

Equation $(\star)$ ties $p$ directly to the public key $A$, the anonymous target $u$ and the (still private) gadget solution $\hat z$.

---

### 5 Final preimage vector

Using the trapdoor, assemble

$$
\begin{aligned}
x_0 &= p_0 + \sum_{j=0}^{k-1} \hat e_j\hat z_j,\\[2mm]
x_1 &= p_1 + \sum_{j=0}^{k-1} \hat r_j\hat z_j,\\[2mm]
x_{j+2} &= p_{j+2} + \hat z_j \quad(j=0,\dots,k-1).
\end{aligned}
$$

Compactly,

$$
x \;=\;\bigl[x_0,\,x_1,\,x_2,\dots,x_{k+1}\bigr]^{\!\top}\in R_q^{k+2}.
$$

Every coordinate is now a **deterministic function** of
$(u,\,A,\,G,\,\hat r,\,\hat e,\,\hat z,\,p)$.

---

### 6 Verification identity

Plugging the definitions of $A$ and $x$ gives

$$
A\cdot x
   \;=\;
     \underbrace{A\cdot p}_{u-G\hat z}
     \;+\;
     \underbrace{\bigl(G\hat z - a\!\sum_j\hat r_j\hat z_j -\sum_j\hat e_j\hat z_j\bigr)}_{\text{added in }A_2\text{ block}}
     \;+\;
     \underbrace{\bigl(a\!\sum_j\hat r_j\hat z_j +\sum_j\hat e_j\hat z_j\bigr)}_{\text{added in }x_0,x_1}.
$$

The two highlighted blocks cancel term-wise, leaving

$$
\boxed{A\cdot x \;=\; u \pmod{q}}.
$$

---

### 7 Parameter map (quick lookup)

| Symbol      | Equation / Constraint                             | Depends on       |
| ----------- | ------------------------------------------------- | ---------------- |
| $k$         | $\lceil\log_t q\rceil$                            | $q,t$            |
| $g_j$       | $t^{\,j}$                                         | $t$              |
| $G$         | $[g_0,\dots,g_{k-1}]$                             | $t,q$            |
| $A_0$       | $1$                                               | —                |
| $A_1$       | $$a\overset{\$}{\leftarrow}R_q$$                  | randomness       |
| $A_2[j]$    | $g_j-(a\hat r_j+\hat e_j)$                        | $G,\,a,\,T$      |
| $p$         | unique small-norm solution of $A p = u - G\hat z$ | $u,A,G,\hat z$   |
| $\hat z$    | discrete-Gaussian solving $G\hat z \equiv u-Ap$   | $u,A,G,p$        |
| $x_0$       | $p_0+\sum \hat e_j\hat z_j$                       | $p,\,T,\hat z$   |
| $x_1$       | $p_1+\sum \hat r_j\hat z_j$                       | $p,\,T,\hat z$   |
| $x_{j+2}$   | $p_{j+2}+\hat z_j$                                | $p,\hat z$       |
| $x$         | concatenation of the above                        | all of the above |
| Final check | $A x \equiv u$                                    | public data only |

Every parameter is therefore expressed either **directly** in terms of $(u,A)$ or **indirectly** through the trapdoor pair and the gadget relation, satisfying the requested linkage.
