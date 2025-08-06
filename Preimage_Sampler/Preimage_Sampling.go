// main.go
package Preimage_Sampler

import (
	"fmt"
	"log"
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// SpectralBound returns the analytic upper-bound
//
//	s_max = 1.8*(t+1)*sigma^2*(sqrt(n*k)+sqrt(2n)+4.7)
func SpectralBound(n, k int, base uint64) float64 {
	const (
		dgError       = 8.27181e-25
		nMax          = 2048
		spectralConst = 1.8
	)
	sigma := math.Sqrt(math.Log(2*float64(nMax)/dgError) / math.Pi)
	sig2 := sigma * sigma
	term := math.Sqrt(float64(n*k)) + math.Sqrt(2*float64(n)) + 4.7
	return spectralConst * float64(base+1) * sig2 * term
}

// calculateParams computes σₜ=(t+1)σ and returns the spectral bound s
// used for perturbation sampling.
func CalculateParams(base uint64, n, k int) (sigmaT, s float64) {
	sigma := 3.19                    // smoothing parameter from Sec V-A1
	sigmaT = float64(base+1) * sigma // σₜ = (t+1)·σ
	s = SpectralBound(n, k, base)    // 1.3 * sigma * sigmaT * sqrt(n*k + 2*n + 4.7)
	return
}

// ZtoZhat converts the integer gadget matrix Z (κ × N) into κ polynomials
// in R_q.  Each row is
//  1. copied to coefficient form,
//  2. sent through the forward NTT,                 (still standard residues)
//  3. multiplied slot-wise by ψ^{-1},               (still standard residues)
//     where ψ = NTT(1) and ψ^{-1} is its entry-wise inverse.
//
// The result slice contains κ polynomials in the *standard* NTT frame
// (no Montgomery factor), ready for later use with ringQ.MulCoeffs.
func ZtoZhat(Z [][]int64, ringQ *ring.Ring) []*ring.Poly {
	k, N := len(Z), ringQ.N
	if k == 0 {
		log.Fatal("empty gadget vector")
	}
	q := ringQ.Modulus[0]
	bredCtx := ringQ.BredParams[0] // Barrett context for q

	// ---------------------------------------------------------------------
	// 1) Pre-compute ψ⁻¹  in standard residues
	// ---------------------------------------------------------------------
	psi := ringQ.NewPoly()
	for i := 0; i < N; i++ { // constant poly "1"
		psi.Coeffs[0][i] = 1
	}
	ringQ.NTT(psi, psi) // ψ(i)  (standard residues)

	psiInv := make([]uint64, N)
	for i := 0; i < N; i++ {
		coeff := psi.Coeffs[0][i]              // already standard
		psiInv[i] = ring.ModExp(coeff, q-2, q) // coeff^{-1}  mod q
	}

	// ---------------------------------------------------------------------
	// 2) Process each row of Z
	// ---------------------------------------------------------------------
	out := make([]*ring.Poly, k)

	for row := 0; row < k; row++ {
		if len(Z[row]) != N {
			log.Fatalf("row %d length %d ≠ %d", row, len(Z[row]), N)
		}

		p := ringQ.NewPoly() // coefficient domain
		for t := 0; t < N; t++ {
			v := Z[row][t] % int64(q) // map into [0, q)
			if v < 0 {
				v += int64(q)
			}
			p.Coeffs[0][t] = uint64(v)
		}

		ringQ.NTT(p, p) // to evaluation domain (standard)

		// slot-wise multiply by ψ⁻¹   --- still standard
		for t := 0; t < N; t++ {
			p.Coeffs[0][t] = ring.BRedConstant( // (p·ψInv) mod q
				p.Coeffs[0][t], psiInv[t], q, bredCtx)
		}
		out[row] = p // standard NTT poly
	}
	return out
}

// GaussSamp implements Alg 2 (discrete branch) from the paper.
//
//	ringQ:    the R_q ring
//	A:        public key row, length k+2
//	rHat,eHat: trapdoor polys (in EVAL) each length k
//	u:        target syndrome in EVAL
//	sigma:    σₜ from calculateParams
//	s:        “large” σ for SamplePz from calculateParams
//	base:     gadget base t
//	k:        gadget length κ
func GaussSamp(
	ringQ *ring.Ring,
	A []*ring.Poly,
	rHat, eHat []*ring.Poly,
	u *ring.Poly,
	sigma, s float64,
	base uint64,
	k int,
) []*ring.Poly {
	N := ringQ.N
	// 1) perturbation: p ∈ Z^{(k+2)×N}
	p := SamplePz(ringQ, s, (2+1)*sigma, [2][]*ring.Poly{rHat, eHat}, k+2, 256)

	// 2) compute sub = u - A·p in EVAL, then back to COEFF
	pertEval := ringQ.NewPoly()
	tmp := ringQ.NewPoly()
	for i := range p {
		ringQ.MulCoeffs(A[i], p[i], tmp)
		ringQ.Add(pertEval, tmp, pertEval)
	}
	subEval := ringQ.NewPoly()
	ringQ.Sub(u, pertEval, subEval)

	sub := ringQ.NewPoly()
	ringQ.InvNTT(subEval, sub) // sub.Coeffs now holds u - A·p in coeffs

	// 3) discrete G-sampling (Alg 3) to get Z ∈ Z^{κ×N}
	//    we only need the level-0 coefficients
	uCoeffs := make([]uint64, N)
	for j := 0; j < N; j++ {
		signed := UnsignedToSigned(sub.Coeffs[0][j], ringQ.Modulus[0])
		uCoeffs[j] = SignedToUnsigned(signed, ringQ.Modulus[0])

	}

	// fmt.Printf("uCoeffs = %v\n", uCoeffs)
	Zmat := SampleGDiscrete(ringQ, (float64(base+1) * sigma), base, uCoeffs, k)

	// Zmat is a matrix of integers Z ∈ ℤ^{κ×N} where each row is a poly in R_q

	// 4) CRT+NTT each row of Z to get zHat ∈ R_q^κ
	zHat := ZtoZhat(Zmat, ringQ)

	fmt.Println("GaussSamp: using transposed trapdoor for inner products")

	// 5) assemble x = [ p₀ + ê⋅zHat,  p₁ + r̂⋅zHat,  p₂+ẑ₀, …, p_{k+1}+ẑ_{k-1} ]
	x := make([]*ring.Poly, k+2)

	// row 0: p[0] + eHat . zHat with eHat . zHat = sum(MulCoeffsAndAdd(eHat[j], zHat[j]))
	sum0 := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		tmpez := ringQ.NewPoly()
		ringQ.MulCoeffs(eHat[j], zHat[j], tmpez)
		ringQ.Add(sum0, tmpez, sum0)
	}
	x[0] = ringQ.NewPoly()
	// In the PALISADE C++ implementation these terms are added
	// when forming the first two rows of x.
	ringQ.Add(p[0], sum0, x[0])

	// row 1: p[1] + zHat with rHat . zHat = sum(MulCoeffsAndAdd(rHat[j], zHat[j]))
	sum1 := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		tmprz := ringQ.NewPoly()
		ringQ.MulCoeffs(rHat[j], zHat[j], tmprz)
		ringQ.Add(sum1, tmprz, sum1)
	}
	x[1] = ringQ.NewPoly()
	ringQ.Add(p[1], sum1, x[1])

	// rows 2…k+1: just p[i] + zHat[i-2]
	for i := 2; i < k+2; i++ {
		x[i] = ringQ.NewPoly()
		ringQ.Add(p[i], zHat[i-2], x[i])
	}
	// -----------------------------------------------------------------------------
	// Collapse tracer – verifies each algebraic block inside A·x
	// -----------------------------------------------------------------------------
	{
		// shorthand helpers -------------------------------------------------------
		centre := func(v uint64) int64 {
			q := int64(ringQ.Modulus[0])
			x := int64(v)
			if x > q/2 {
				x -= q
			}
			return x
		}
		mustBeZero := func(name string, poly *ring.Poly) {
			coeff := ringQ.NewPoly()
			ringQ.InvNTT(poly, coeff)
			for t := 0; t < ringQ.N; t++ {
				if centre(coeff.Coeffs[0][t]) != 0 {
					log.Fatalf("[TRACE] %s mismatch at slot %d: %d",
						name, t, centre(coeff.Coeffs[0][t]))
				}
			}
			fmt.Printf("[TRACE] %s  ✔\n", name)
		}

		// -------------------------------------------------------------------------
		// Part 0 : pre-compute reusable pieces
		// -------------------------------------------------------------------------
		ApEval := ringQ.NewPoly() //  A·p  (already needed later)
		tmpEval := ringQ.NewPoly()
		for i := range p {
			ringQ.MulCoeffs(A[i], p[i], tmpEval)
			ringQ.Add(ApEval, tmpEval, ApEval)
		}

		sumEZEval := ringQ.NewPoly()
		ring.Copy(sum0, sumEZEval) // ⟨ê̂,ẑ⟩
		sumRZEval := ringQ.NewPoly()
		ring.Copy(sum1, sumRZEval) // ⟨r̂,ẑ⟩

		// Σ_j A₂[j]·ẑ_j
		A2zEval := ringQ.NewPoly()
		for j := 0; j < k; j++ {
			ringQ.MulCoeffs(A[j+2], zHat[j], tmpEval)
			ringQ.Add(A2zEval, tmpEval, A2zEval)
		}

		// -------------------------------------------------------------------------
		// 1)  A[0]·x[0]  ==  p0 + Σ ê_j·ẑ_j
		// -------------------------------------------------------------------------
		left0 := ringQ.NewPoly()
		ringQ.MulCoeffs(A[0], x[0], left0) // A0 = 1, but keep generic
		right0 := ringQ.NewPoly()
		ringQ.Add(p[0], sumEZEval, right0) // p0 + ⟨ê̂,ẑ⟩

		diff0 := ringQ.NewPoly()
		ringQ.Sub(left0, right0, diff0)
		mustBeZero("A0·x0 == p0+e·z", diff0)

		// -------------------------------------------------------------------------
		// 2)  A[1]·x[1]  ==  a·p1 + a·Σ r̂_j·ẑ_j
		// -------------------------------------------------------------------------
		left1 := ringQ.NewPoly()
		ringQ.MulCoeffs(A[1], x[1], left1)

		right1 := ringQ.NewPoly()
		ringQ.MulCoeffs(A[1], p[1], right1) // a·p1
		tmpEval = ringQ.NewPoly()
		ringQ.MulCoeffs(A[1], sumRZEval, tmpEval) // a·⟨r̂,ẑ⟩
		ringQ.Add(right1, tmpEval, right1)

		diff1 := ringQ.NewPoly()
		ringQ.Sub(left1, right1, diff1)
		mustBeZero("A1·x1 == a·p1+a·r·z", diff1)

		// // -------------------------------------------------------------------------
		// // 3)  Σ_j(g_j−a r̂_j−ê_j)·ẑ_j =  A·p − u − a·r·z − e·z
		// // -------------------------------------------------------------------------
		// rhs := ringQ.NewPoly()
		// ringQ.Sub(ApEval, u, rhs)      //  A·p − u
		// ringQ.Sub(rhs, sumEZEval, rhs) // - e·z
		// tmpEval = ringQ.NewPoly()
		// ringQ.MulCoeffs(A[1], sumRZEval, tmpEval) // a·r·z
		// ringQ.Sub(rhs, tmpEval, rhs)              // - a·r·z

		// diff2 := ringQ.NewPoly()
		// ringQ.Sub(A2zEval, rhs, diff2)
		// mustBeZero("Σ(g−ar̂−ê)·ẑ == A·p-u-a·r·z-e·z", diff2)
	}
	// -----------------------------------------------------------------------------
	// End of collapse tracer
	// -----------------------------------------------------------------------------

	// //! ---------- DEBUG B : verify Ap + Gz = u ----------
	// if true {
	// 	// helper to centre a coeff into (−q/2, q/2]
	// 	centre := func(v uint64) int64 {
	// 		q := int64(ringQ.Modulus[0])
	// 		x := int64(v)
	// 		if x > q/2 {
	// 			x -= q
	// 		}
	// 		return x
	// 	}

	// 	tmp := ringQ.NewPoly()

	// 	// 1) Compute A·x in NTT, then back to coefficients
	// 	AxEval := ringQ.NewPoly()
	// 	for i := range A {
	// 		ringQ.MulCoeffs(A[i], x[i], tmp)
	// 		ringQ.Add(AxEval, tmp, AxEval)
	// 	}

	// 	// 1b) also compute A·p and A·z separately for diagnostics
	// 	ApEval := ringQ.NewPoly()
	// 	for i := range p {
	// 		ringQ.MulCoeffs(A[i], p[i], tmp)
	// 		ringQ.Add(ApEval, tmp, ApEval)
	// 	}
	// 	ApCoeff := ringQ.NewPoly()
	// 	ringQ.InvNTT(ApEval, ApCoeff)

	// 	AzEval := ringQ.NewPoly()
	// 	// contributions from gadget columns
	// 	for j := 0; j < k; j++ {
	// 		ringQ.MulCoeffs(A[j+2], zHat[j], tmp)
	// 		ringQ.Add(AzEval, tmp, AzEval)
	// 	}
	// 	// contributions from x0,x1 parts
	// 	eDotZEvalAcc := ringQ.NewPoly()
	// 	rDotZEvalAcc := ringQ.NewPoly()
	// 	for j := 0; j < k; j++ {
	// 		ringQ.MulCoeffs(eHat[j], zHat[j], tmp)
	// 		ringQ.Add(eDotZEvalAcc, tmp, eDotZEvalAcc)
	// 		ringQ.MulCoeffs(rHat[j], zHat[j], tmp)
	// 		ringQ.Add(rDotZEvalAcc, tmp, rDotZEvalAcc)
	// 	}
	// 	ringQ.MulCoeffs(A[0], eDotZEvalAcc, tmp)
	// 	ringQ.Add(AzEval, tmp, AzEval)
	// 	ringQ.MulCoeffs(A[1], rDotZEvalAcc, tmp)
	// 	ringQ.Add(AzEval, tmp, AzEval)

	// 	AzCoeff := ringQ.NewPoly()
	// 	ringQ.InvNTT(AzEval, AzCoeff)

	// 	// 1c) verify that A·p + G·z equals the target u
	// 	Gz := CreateGadgetMatrix(ringQ, base, 1, k)
	// 	GzEval := ringQ.NewPoly()
	// 	for j := 0; j < k; j++ {
	// 		ringQ.NTT(Gz[j], Gz[j])
	// 		ringQ.MulCoeffs(Gz[j], zHat[j], tmp)
	// 		ringQ.Add(GzEval, tmp, GzEval)
	// 	}

	// 	GzCoeff_FromVector := ringQ.NewPoly()
	// 	ringQ.InvNTT(GzEval, GzCoeff_FromVector)

	// 	uCoeff := ringQ.NewPoly()
	// 	ringQ.InvNTT(u, uCoeff)

	// 	sumCheck := ringQ.NewPoly()
	// 	ringQ.Add(ApCoeff, GzCoeff_FromVector, sumCheck)
	// 	ringQ.Sub(sumCheck, uCoeff, sumCheck)

	// 	diff0 := centre(sumCheck.Coeffs[0][0])
	// 	if diff0 != 0 {
	// 		log.Printf("[CHECK] A·p coeff[0]=%d", centre(ApCoeff.Coeffs[0][0]))
	// 		log.Printf("[CHECK] G·z coeff[0]=%d", centre(GzCoeff_FromVector.Coeffs[0][0]))
	// 		log.Printf("[CHECK] u coeff[0]=%d", centre(uCoeff.Coeffs[0][0]))
	// 		log.Fatalf("[CHECK] A·p+G·z mismatch at slot 0: %d", diff0)
	// 	}
	// 	for t := 1; t < ringQ.N; t++ {
	// 		if centre(sumCheck.Coeffs[0][t]) != 0 {
	// 			log.Fatalf("[CHECK] A·p+G·z mismatch at slot %d: %d", t, centre(sumCheck.Coeffs[0][t]))
	// 		}
	// 	}
	// 	fmt.Println("[CHECK] A·p + G·z ≡ u  ✔")

	// 	// 1d) verify that the object per object assembled x indeed satisfies A·x = u
	// 	AxEvalCheck := ringQ.NewPoly()
	// 	ring.Copy(ApEval, AxEvalCheck)

	// 	blockA2 := ringQ.NewPoly()
	// 	ring.Copy(GzEval, blockA2)
	// 	ringQ.MulCoeffs(A[1], rDotZEvalAcc, tmp)
	// 	ringQ.Sub(blockA2, tmp, blockA2)
	// 	ringQ.Sub(blockA2, eDotZEvalAcc, blockA2)

	// 	blockX01 := ringQ.NewPoly()
	// 	ringQ.MulCoeffs(A[1], rDotZEvalAcc, blockX01)
	// 	ringQ.Add(blockX01, eDotZEvalAcc, blockX01)

	// 	ringQ.Add(AxEvalCheck, blockA2, AxEvalCheck)
	// 	ringQ.Add(AxEvalCheck, blockX01, AxEvalCheck)

	// 	AxCoeff := ringQ.NewPoly()
	// 	ringQ.InvNTT(AxEvalCheck, AxCoeff)

	// 	diffAx := ringQ.NewPoly()
	// 	ringQ.Sub(AxCoeff, uCoeff, diffAx)

	// 	for t := 1; t < ringQ.N; t++ {
	// 		if centre(diffAx.Coeffs[0][t]) != 0 {
	// 			log.Fatalf("[CHECK] A·x mismatch at slot %d: %d", t, centre(diffAx.Coeffs[0][t]))
	// 		}
	// 	}
	// 	fmt.Println("[CHECK] A·x ≡ u  ✔")

	// 	// 6a) accumEval = A ⋅ x in evaluation domain
	// 	accumEval := ringQ.NewPoly()
	// 	tmpEval := ringQ.NewPoly()
	// 	for i := 0; i < len(A); i++ {
	// 		ringQ.MulCoeffs(A[i], x[i], tmpEval)
	// 		ringQ.Add(accumEval, tmpEval, accumEval)
	// 	}

	// 	// 6b) Verify that A·x ≡ u mod q
	// 	for i := 0; i < ringQ.N; i++ {
	// 		if accumEval.Coeffs[0][i] != u.Coeffs[0][i] {
	// 			log.Fatalf("[GaussSamp] A·x mismatch at slot %d: got %d want %d",
	// 				i, accumEval.Coeffs[0][i], u.Coeffs[0][i])
	// 		}
	// 	}
	// }
	// //! ---------- END DEBUG B ----------
	return x
}

func Main() {
	// ----------------------------------------------------------------------------
	// 1) Parameters
	// ----------------------------------------------------------------------------
	n := 512                //64                // ring dimension (must match paper)
	qMod := uint64(8399873) //uint64(114689) // ciphertext modulus
	base := uint64(2)       // gadget base t
	k := int(math.Ceil(math.Log(float64(qMod)) / math.Log(float64(base))))

	// σₜ (sigmaT) and spectral bound s from our parameter calculation
	sigmaT, s := CalculateParams(base, n, k)
	// small Gaussian width σ = σₜ / (t+1)
	sigma := sigmaT / float64(base+1)

	// ----------------------------------------------------------------------------
	// 2) Setup ring
	// ----------------------------------------------------------------------------
	ringQ, err := ring.NewRing(n, []uint64{qMod})
	if err != nil {
		panic(err)
	}

	// ----------------------------------------------------------------------------
	// 3) Trapdoor generation
	// ----------------------------------------------------------------------------
	trap := TrapGen(ringQ, base, sigmaT)

	// ----------------------------------------------------------------------------
	// 4) Sample a random syndrome u ∈ R_q in EVAL
	// ----------------------------------------------------------------------------
	u := ringQ.NewPoly()
	for i := 0; i < ringQ.N; i++ {
		u.Coeffs[0][i] = uint64(i % int(qMod))
	}
	uEval := ringQ.NewPoly() // u in evaluation domain
	ringQ.NTT(u, uEval)

	// ----------------------------------------------------------------------------
	// 5) Gaussian preimage sampling: find x s.t. A·x ≈ u
	// ----------------------------------------------------------------------------
	x := GaussSamp(ringQ, trap.A, trap.R[0], trap.R[1], uEval, sigma, s, base, k)

	// ----------------------------------------------------------------------------
	// 6) Verification: compute A·x and check exact equality mod q
	// ----------------------------------------------------------------------------
	// 6a) accumEval = A ⋅ x in evaluation domain
	accumEval := ringQ.NewPoly()
	tmpEval := ringQ.NewPoly()
	for i := 0; i < len(trap.A); i++ {
		ringQ.MulCoeffs(trap.A[i], x[i], tmpEval)
		ringQ.Add(accumEval, tmpEval, accumEval)
	}

	// 6b) Verify that A·x ≡ u mod q
	for i := 0; i < ringQ.N; i++ {
		if accumEval.Coeffs[0][i] != uEval.Coeffs[0][i] {
			log.Fatalf("A·x mismatch at slot %d: got %d want %d",
				i, accumEval.Coeffs[0][i], uEval.Coeffs[0][i])
		}
	}

	fmt.Println("SUCCESS: exact A·x ≡ u mod q")
}
