# Package lvcs

Linear‐map Vector Commitment Scheme (LVCS), built atop the DECS small-domain PCS.
Implements §4.1 of “SmallWood: Hash-Based Polynomial Commitments and Zero-Knowledge Arguments…”:

1. **Commit**: lift row-vectors to polynomials with random masks; commit via DECS.
2. **Eval**: open masked linear combinations on a small random subset; verify Merkle, degree, and linear checks, leaking only masked positions.

## Files

* **lvcs\_prover.go**
  Prover‐side routines: `CommitInit`, `CommitFinish`, `EvalInit`, `EvalFinish`.

* **lvcs\_verifier.go**
  Verifier‐side routines: `NewVerifier`, `CommitStep1`, `CommitStep2`, `ChooseE`, `EvalStep2`.

* **Interpolate.go**
  Helper: `interpolateRow`—Lagrange‐style interpolation over the first $ncols+\ell$ NTT points.

* **lvcs\_test.go**
  End-to-end test driving one full commit + eval.

---

## lvcs\_prover.go

```go
// ProverKey holds all state the prover needs between Commit and Eval.
type ProverKey struct {
    RingQ      *ring.Ring   // NTT ring context, for modulus and roots
    DecsProver *decs.Prover // underlying DECS prover instance

    RowData  [][]uint64 // the original row vectors r_j ∈ F_q^{ncols}
    MaskData [][]uint64 // the ℓ random mask vectors \bar r_j
}
```

### CommitInit

```go
func CommitInit(
    ringQ *ring.Ring,
    rows [][]uint64,  // r_j length = ncols = N - ℓ
    ell int,          // ℓ = mask size
) (root [32]byte, prover *ProverKey, err error)
```

1. **Sample Masks**: for each row $j$, pick $\ell$ uniform scalars in $\F_q$.
2. **Interpolate** each $(r_j, \bar r_j)$ pair to a degree-<$(ncols+ℓ)$ polynomial $P_j(X)$ via `interpolateRow`.
3. **DECS Commit**: `decs.NewProver(ringQ, polys).CommitInit()` → 32-byte Merkle root.
4. Returns `root` and a `ProverKey` containing DECS state + row/mask data.

### CommitFinish

```go
func CommitFinish(prover *ProverKey, Gamma [][]uint64) []*ring.Poly
```

* Inputs: the random DECS challenge $\Gamma\in\F_q^{η\times nrows}$.
* Calls `prover.DecsProver.CommitStep2(Gamma)` → returns the vector of $\eta$ combined polynomials $R_k(X)$.

### EvalInit

```go
func EvalInit(
    ringQ *ring.Ring,
    prover *ProverKey,
    C [][]uint64,   // coefficient matrix size m×nrows
) [][]uint64
```

* Computes masked targets
  $\bar v_k[i]=\sum_{j=1}^{nrows}C_{k,j}\,\bar r_{j}[i]\mod q$.
* Returns `bar` of shape m×ℓ.

### EvalFinish

```go
func EvalFinish(prover *ProverKey, E []int) *decs.DECSOpening
```

* Inputs: challenge set `E` (ℓ distinct indices in \[0..N)).
* Delegates to `prover.DecsProver.EvalOpen(E)` to produce the DECS opening (Merkle paths, evaluations, nonces).

---

## Interpolate.go

```go
func interpolateRow(
    ringQ *ring.Ring,
    row []uint64,    // length ncols
    mask []uint64,   // length ℓ
    ncols, ell int,
) (*ring.Poly, error)
```

Performs Lagrange interpolation of degree <$ncols+ℓ$ on the points

$$
  \{(\omega^i,\,\text{row}[i])\}_{i=0..ncols-1}\;\cup\;\{(\omega^{ncols+i},\,\text{mask}[i])\}_{i=0..ℓ-1},
$$

where $\omega$ is the primitive NTT root of unity.

1. Builds the vanishing polynomial $T(X)=\prod_{i=0..m-1}(X-\omega^i)$, $m=ncols+ℓ$.
2. For each $i$, synthetic-divide $T$ by $(X-\omega^i)$ to get $Q_i(X)$.
3. Compute barycentric denominator $\prod_{j≠i}(\omega^i-\omega^j)$ and invert mod q.
4. Accumulate $P(X)=\sum_i y_i·Q_i(X)/\denom_i$.
5. Zero-pad to full ring length.

Returns a `*ring.Poly` containing the coefficient vector.

---

## lvcs\_verifier.go

```go
type VerifierState struct {
    RingQ  *ring.Ring  // NTT context
    r, eta int         // number of rows, DECS η

    Root  [32]byte     // Merkle root from CommitInit
    Gamma [][]uint64   // DECS challenge matrix η×r
    R     []*ring.Poly // polynomials R_k(X) from CommitFinish
}
```

### NewVerifier

```go
func NewVerifier(ringQ *ring.Ring, r, eta int) *VerifierState
```

Creates a fresh verifier with no commitment recorded yet.

### CommitStep1

```go
func (v *VerifierState) CommitStep1(root [32]byte)
```

Records `root`, then runs `decs.NewVerifier(ringQ,r,eta).DeriveGamma(root)` to sample the same $\Gamma$.

### CommitStep2

```go
func (v *VerifierState) CommitStep2(R []*ring.Poly)
```

Stores the prover’s $R_k(X)$ polynomials.

### ChooseE

```go
func (v *VerifierState) ChooseE(ell int) []int
```

Returns ℓ uniform random indices in $[0..N-1]$ using Go’s `math/rand`.

### EvalStep2

```go
func (v *VerifierState) EvalStep2(
    bar [][]uint64,         // prover’s masked targets \bar v_k
    E []int,                // challenge subset ℓ
    open *decs.DECSOpening, // DECS opening (paths + values)
    C [][]uint64,           // coefficient matrix
) bool
```

1. **DECS Verify**: calls `decs.NewVerifier(...).VerifyEval(root,Gamma,R,open)`.
2. **Linear Check**: let

   * `ncols = N - ℓ`.
   * For each opened index `idx` in `open.Indices`:

     * if `idx < ncols`, skip (these leak real data).
     * else `t = idx - ncols` is the masked‐slot; check
       $\sum_j C_{k,j}\,P_j(\omega^{idx})\;\stackrel?=\;bar[k][t]$
       for all k.
3. Return `true` iff all checks pass.

---

## lvcs\_test.go

```go
func TestLVCSCommitAndEval(t *testing.T)
```

An end-to-end sanity check:

1. **Setup**

   * Choose $N=2^{11}$, prime mod ≡1 mod 2N, build `ringQ`.
   * Pick `nrows=4`, `ell=8`.
   * Sample `rows[j]` ∈ F\_q^{N-ell} uniformly.

2. **Commit Phase**

   * Prover: `root, proverKey := CommitInit(ringQ, rows, ell)`.
   * Verifier: `ver := NewVerifier(ringQ,nrows,decs.Eta); ver.CommitStep1(root)`.
   * Prover: `R := CommitFinish(proverKey, ver.Gamma)`.
   * Verifier: `ver.CommitStep2(R)`.

3. **Eval Phase**

   * Pick `m=3` and sample random coefficient matrix `C[m][nrows]`.
   * Prover: `bar := EvalInit(ringQ,proverKey,C)`.
   * Verifier: `E := ver.ChooseE(ell)`.
   * Prover: `opening := EvalFinish(proverKey, E)`.
   * Verifier: `ok := ver.EvalStep2(bar, E, opening, C)`.
   * `t.Fatalf` if `!ok`; otherwise report success.

---

## Mathematical Rationale

* **Binding**: the DECS Merkle + degree check binds the prover to *one* set of polynomials $\{P_j\}$.
* **Soundness of open**: if the prover cheats on $\bar v_k$, then for a random subset $E$, the Schwartz–Zippel bound forces detection.
* **Zero-knowledge**: only ℓ positions are opened, but those ℓ were masked by random $\bar r_j$; the real row data $r_j$ remain hidden.

Together, this realizes a zero-knowledge, binding LVCS exactly as required by §4.1.

