package pcs

import (
	"errors"
	"math/rand"
	"time"

	lvcs "vSIS-Signature/LVCS"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// CommitInitPolynomial commits the vector of polynomials P_j(X) (with coefficients P[j][i])
// under degrees Dj[j], using block‐size Mu and mask‐rows EllPrime.
// Returns a 32‐byte root plus a ProverKey for subsequent steps.
func CommitInitPolynomial(
	ringQ *ring.Ring,
	P [][]uint64, // P[j] holds the Dj[j]+1 coefficients of the j-th polynomial
	Dj []int,
	Mu, EllPrime int,
) (Commitment, *ProverKey, error) {

	npcs := len(P)
	if len(Dj) != npcs {
		return Commitment{}, nil, errors.New("pcs: length mismatch P vs Dj")
	}

	// 1) Compute ν_j for each polynomial
	Nu := make([]int, npcs)
	sumNu := 0
	for j := 0; j < npcs; j++ {
		d := Dj[j] + 1 - EllPrime
		if d <= 0 {
			Nu[j] = 1
		} else {
			Nu[j] = (d + Mu - 1) / Mu
		}
		sumNu += Nu[j]
	}

	// 2) Build the (Mu+EllPrime)×sumNu row‐matrix by concatenating each A_j
	nrows := Mu + EllPrime
	rows := make([][]uint64, nrows)
	for r := 0; r < nrows; r++ {
		rows[r] = make([]uint64, sumNu)
	}

	colOff := 0
	rand.Seed(time.Now().UnixNano())
	q0 := ringQ.Modulus[0] // assume a single modulus

	for j := 0; j < npcs; j++ {
		nu := Nu[j]

		// Fill the first Mu rows with stride‐Mu blocks of P[j]
		for r := 0; r < Mu; r++ {
			for k := 0; k < nu; k++ {
				idx := r + k*Mu
				if idx < len(P[j]) {
					rows[r][colOff+k] = P[j][idx]
				} else {
					rows[r][colOff+k] = 0
				}
			}
		}

		// Fill the EllPrime mask rows with fresh uniform scalars
		for r := 0; r < EllPrime; r++ {
			for k := 0; k < nu; k++ {
				rows[Mu+r][colOff+k] = rand.Uint64() % q0
			}
		}

		colOff += nu
	}

	// 3) Delegate to LVCS.Commit
	root, lvcsKey, err := lvcs.CommitInit(ringQ, rows, EllPrime)
	if err != nil {
		return Commitment{}, nil, err
	}

	// 4) Package prover state
	pk := &ProverKey{
		RingQ:    ringQ,
		LvcsKey:  lvcsKey,
		Dj:       Dj,
		Mu:       Mu,
		EllPrime: EllPrime,
		Nu:       Nu,
	}

	return Commitment{Root: root}, pk, nil
}

// CommitFinishPolynomial runs the second pass of the commit protocol,
// taking the challenge Γ from the verifier and returning the R_k(X) polynomials.
func CommitFinishPolynomial(pk *ProverKey) ([]*ring.Poly, error) {
	Gamma := pk.LvcsKey.Gamma
	Rpolys := lvcs.CommitFinish(pk.LvcsKey, Gamma)
	return Rpolys, nil
}

// EvalInitPolynomial computes the masked linear‐map targets
// bar[k][i] = Σ_j C[k][j]·mask[j][i], for i=0..ℓ′-1.
func EvalInitPolynomial(
	pk *ProverKey,
	C [][]uint64, // m × (Mu+EllPrime)*β rows
) ([][]uint64, error) {
	return lvcs.EvalInit(pk.RingQ, pk.LvcsKey, C), nil
}

// EvalFinishPolynomial runs the LVCS opening on the selected subset E,
// returning the DECS opening (paths, values, nonces).
func EvalFinishPolynomial(
	pk *ProverKey,
	E []int,
) (*lvcs.Opening, error) {
	open := lvcs.EvalFinish(pk.LvcsKey, E)
	return open, nil
}
