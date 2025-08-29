
## 1. Shared Types and Constants (`decs_types.go`)

```go
// DECSOpening holds the prover’s opening in DECS.Eval.
type DECSOpening struct {
  Indices []int      // E ⊆ {0,…,N−1}: which coordinates we open
  Pvals   [][]uint64 // P₀(e),…,P_{r−1}(e) for each e∈E
  Mvals   [][]uint64 // M₀(e),…,M_{η−1}(e) for each e∈E
  Paths   [][][]byte // Merkle authentication paths for each e
  Nonces  [][]byte   // ρ_e used when hashing leaf e
}

// protocol-wide parameters:
const (
  Eta        = 2    // η, number of mask polynomials
  Degree     = 4095 // d, maximum degree for P and M
  NonceBytes = 16   // λ, bytes in each per-leaf nonce ρ
)
```

* **`DECSOpening`**: captures everything the prover reveals in DECS.Eval:

  * the chosen evaluation indices `Indices`,
  * the raw polynomial values `Pvals` and `Mvals`,
  * the Merkle “paths” proving those values were committed under the root,
  * and the random nonces used per-leaf to prevent second-preimage attacks on the leaf hash.
* **`Eta`, `Degree`, `NonceBytes`** exactly mirror the DECS paper’s parameters:
  η masks of degree ≤ d, and λ-byte randomness.

---

## 2. Prover (`decs_prover.go`)

```go
type Prover struct {
  ringQ  *ring.Ring    // NTT context over Z/qZ
  P      []*ring.Poly  // r input polynomials P₀,…,P_{r−1} in coeff form
  M      []*ring.Poly  // η random mask polynomials (coeff form)
  nonces [][]byte      // per-leaf nonces ρ₀…ρ_{N−1}
  mt     *MerkleTree   // authenticated over leaves u₀…u_{N−1}
  Pvals  []*ring.Poly  // NTT(P)
  Mvals  []*ring.Poly  // NTT(M)
  root   [32]byte      // Merkle root
  R      []*ring.Poly  // η final polynomials R₀…R_{η−1} (coeff form)
}
```

* **`ringQ`** encapsulates modulus q, ring dimension N, primitive root for the NTT, etc.
* **`P`** is the vector of input polynomials ∈ (F\[X]≤d)ʳ.
* **`M`** are the prover’s random degree-d masks.
* **`Pvals`/`Mvals`** store the evaluation-form (NTT) of P/M so we can cheaply read off P(e) and M(e).

---

### `NewProver(ringQ, P)`

Allocates a `Prover` struct binding us to a particular ring and witness polynomials P₀…P\_{r−1}.

---

### `CommitInit() (root, error)`

Implements **DECS.Commit Step 1**:

1. **Sample masks**
   We draw η polynomials M₀…M\_{η−1} uniformly at random from degree-≤ d, then zero out all coefficients above `Degree`.

2. **NTT-transform**
   Transform each P\_j and M\_k from coefficient to evaluation form (`Pvals`, `Mvals`). After this,

   $$
     Pvals[j].Coeffs[0][i] = P_j(e_i),\quad Mvals[k].Coeffs[0][i] = M_k(e_i)
   $$

   where `e_i` is the i-th NTT root of unity.

3. **Build leaves**
   For each $i\in[0,N)$ we form a buffer

   $$
     \text{buf} = \bigl(P_0(e_i)\|…\|P_{r-1}(e_i)\|\;M_0(e_i)\|…\|M_{η-1}(e_i)\|\;i\|\;\rho_i\bigr)
   $$

   encoding each $P_j(e_i)$ and $M_k(e_i)$ as 64-bit little-endian integers, the index $i$ as 32-bit little-endian, and appending a fresh $\rho_i\in\{0,1\}^{\lambda}$. We hash this buffer via SHA-256 to get a 32-byte leaf.

4. **Merkle tree over leaves**
   We pad the leaf array to a power of two (filling missing leaves with `H(nil)`), then build all upper layers by hashing each pair of siblings. The root is saved in `pr.root`.

Finally, `CommitInit` returns that 32-byte root, which the prover sends to the verifier.

---

### `CommitStep2(Gamma) []*ring.Poly`

Implements **DECS.Commit Steps 2–3**:

* **Input**: the challenge matrix $\Gamma\in\F_q^{\eta\times r}$ supplied by the verifier.
* **Goal**: compute, for each $k=0…\eta-1$,

  $$
    R_k(X) \;=\; M_k(X)\;+\;\sum_{j=0}^{r-1}\Gamma_{k,j}\,P_j(X)\,.
  $$
* **Implementation**:

  1. Inverse-NTT `Mvals[k]` back to coefficient form → `tmp`.
  2. Copy into `R[k]`.
  3. For each j: inverse-NTT `Pvals[j]`→`tmp`, multiply `tmp` by scalar $\Gamma[k][j]$ (via `MulScalar`), add into `R[k]`.
  4. NTT `R[k]` → evaluation form for future checks.

Returns the slice `R = [R₀,…,R_{η−1}]`, which the prover sends to the verifier.

---

### `EvalOpen(E []int) *DECSOpening`

Implements **DECS.Eval Step 1**. Given a subset $E = \{e_{t}\}_{t=0}^{\ell-1}$:

* For each $t$,

  * `open.Pvals[t][j] = P_j(e_t)` (read off from `pr.Pvals[j].Coeffs[0][e_t]`),
  * `open.Mvals[t][k] = M_k(e_t)`,
  * `open.Paths[t] = pr.mt.Path(e_t)` is the Merkle authentication path,
  * `open.Nonces[t] = pr.nonces[e_t]`.

The returned `DECSOpening` bundles exactly those pieces the verifier needs to check correctness.

---

### `deriveGamma(root, η, r)`

Deterministically expands the 32-byte root into an $\eta\times r$ matrix over $\F_q$ by:

* Appending a 4-byte counter to the 32-byte root,
* Hashing with SHA-256,
* Taking the first 4 bytes as a little-endian `uint32`,
* Repeating until $\eta\cdot r$ values are produced.

This replaces an interactive coin-flip step in a single round (via Fiat–Shamir in practice).

---

## 3. Verifier (`decs_verifier.go`)

```go
type Verifier struct {
  ringQ *ring.Ring
  r, eta int
}
```

Holds the ring context and the dimensions $r$ (# of P-polys) and $\eta$.

---

### `NewVerifier(ringQ, r, eta)`

Simple constructor.

---

### `DeriveGamma(root)`

Identical to the prover’s `deriveGamma`, implements **DECS.Commit Step 2**: sampling the challenge matrix.

---

### `VerifyCommit(root, R, Gamma)`

A no-op here (in a full implementation one might check bounds on R). We record `(root, R, Γ)` for the final check.

---

### `VerifyEval(root, Γ, R, open) bool`

Implements **DECS.Eval Step 2**:

1. **Pre-NTT** each `R[k]` back to evaluation form → `Re[k]`.
2. For each opened index `e = open.Indices[t]`:

   * **Reconstruct the leaf buffer** exactly as the prover did (64-bit little-endian values, 32-bit index, nonce) and run `VerifyPath` with the provided `open.Paths[t]`. This ensures that the revealed $(P_j(e),M_k(e))$ truly came from the committed Merkle root.
   * **Check the masked relation**: compute

     $$
       \mathrm{lhs} = R_k(e)\quad(\text{from }Re[k].Coeffs[0][e]), 
       \quad
       \mathrm{rhs} \;=\; M_k(e)\;+\;\sum_{j}\Gamma_{k,j}\,P_j(e)\quad(\bmod q).
     $$

     Reject on any mismatch.

If all leaves authenticate and all relations hold, return true.

---

## 4. Merkle Tree (`merkle.go`)

A minimal, power-of-two Merkle‐tree over 32-byte leaves:

* **`BuildMerkleTree(leaves [][]byte)`**

  1. Pad leaf list to the next power of two by hashing `nil` for missing spots.
  2. `layers[0] = sha256(leaves[i])`.
  3. For each upper layer, `node = H(left || right)`.

* **`Root()`** returns `layers[last][0]`.

* **`Path(idx)`** gathers sibling hashes at each level.

* **`VerifyPath(leaf, path, root, idx)`** re-computes the path from `leaf` to `root`, taking left/right order from the low bit of `idx`.

---

## 5. End-to-End Test (`decs_test.go`)

```go
func TestDECSProtocol(t *testing.T) {
  // 1) ring.NewRing(N, moduli)
  // 2) Sample r random P_j ∈ F[X]≤d via Lattigo’s UniformSampler, zero high coeffs.
  // 3) Prover.CommitInit() → root
  // 4) Verifier.DeriveGamma(root) → Γ
  // 5) Prover.CommitStep2(Γ) → R
  // 6) Verifier.VerifyCommit(root,R,Γ)
  // 7) Pick ℓ random indices E
  // 8) Prover.EvalOpen(E)
  // 9) Verifier.VerifyEval(root,Γ,R,opening)
}
```

At each step we exercise exactly the three transcripts of DECS:

1. **CommitInit** (Prover→root)
2. **Challenge** (Verifier→Γ) + **CommitStep2** (Prover→R)
3. **EvalOpen** (Prover→open) + **VerifyEval** (Verifier→bool)

and assert the final `bool` is true.

---

### Why it all fits the DECS paper

* **Low-degree masking (M)**: enforces zero-knowledge.
* **Merkle commitment** to $(P(e),M(e))$ for all $e$.
* **Random challenge** $\Gamma$ tests linear combinations of $P$’s evaluations.
* **R-polynomials** bind the prover to a single underlying polynomial vector.
* **Open a small subset $E$** so proof size scales with $|E|$, not full $N$.

All underlying NTT, hash, and Merkle operations exactly mirror the “hash + Merkle tree + small‐domain PCS” construction in the paper.
