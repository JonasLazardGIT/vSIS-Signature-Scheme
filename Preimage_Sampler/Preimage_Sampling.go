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

func ZtoZhat(Z [][]int64, ringQ *ring.Ring) []*ring.Poly {
	k, N := len(Z), ringQ.N
	if k == 0 {
		log.Fatal("empty gadget vector")
	}

	q := ringQ.Modulus[0]

	// psiInv holds the coefficient-wise inverse of NTT(1).
	psi := ringQ.NewPoly()
	for i := 0; i < N; i++ {
		psi.Coeffs[0][i] = 1
	}
	ringQ.NTT(psi, psi)
	psiInv := make([]uint64, N)
	for i := 0; i < N; i++ {
		standard := ring.InvMFormConstant(psi.Coeffs[0][i], q, ringQ.MredParams[0])
		invStd := ring.ModExp(standard, q-2, q)
		psiInv[i] = ring.MForm(invStd, q, ringQ.BredParams[0])
	}

	out := make([]*ring.Poly, k)

	for row := 0; row < k; row++ {
		if len(Z[row]) != N {
			log.Fatalf("row %d length %d ≠ %d", row, len(Z[row]), N)
		}

		p := ringQ.NewPoly() // coefficient domain
		for t := 0; t < N; t++ {
			v := Z[row][t] % int64(q) // canonical 0…q-1
			if v < 0 {
				v += int64(q)
			}
			p.Coeffs[0][t] = uint64(v)
		}
		ringQ.NTT(p, p) // to evaluation domain

		// multiply by psi^{-1} coefficient-wise so that G·zHat matches
		for t := 0; t < N; t++ {
			p.Coeffs[0][t] = ring.MRedConstant(p.Coeffs[0][t], psiInv[t], q, ringQ.MredParams[0])
		}

		out[row] = p
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
		ringQ.MulCoeffsMontgomery(A[i], p[i], tmp)
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

	// row 0: p[0] + eHat . zHat with eHat . zHat = sum(MulCoeffsMontgomeryAndAdd(eHat[j], zHat[j]))
	sum0 := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		tmpez := ringQ.NewPoly()
		ringQ.MulCoeffsMontgomery(eHat[j], zHat[j], tmpez)
		ringQ.Add(sum0, tmpez, sum0)
	}
	x[0] = ringQ.NewPoly()
	// In the PALISADE C++ implementation these terms are added
	// when forming the first two rows of x.
	ringQ.Add(p[0], sum0, x[0])

	// row 1: p[1] + zHat with rHat . zHat = sum(MulCoeffsMontgomeryAndAdd(rHat[j], zHat[j]))
	sum1 := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		tmprz := ringQ.NewPoly()
		ringQ.MulCoeffsMontgomery(rHat[j], zHat[j], tmprz)
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
	// Sanity-check block  (insert in GaussSamp right after sum0 / sum1 are built)
	// -----------------------------------------------------------------------------
	{
		// -------------------------------------------------------------------------
		// 1)  Verify   Σ ê·ẑ    and    Σ r̂·ẑ   against a fresh recomputation
		// -------------------------------------------------------------------------
		tmp0 := ringQ.NewPoly()
		tmp1 := ringQ.NewPoly()
		ringQ.InvNTT(sum0, tmp0) // tmp0 = ⟨ê , ẑ⟩
		ringQ.InvNTT(sum1, tmp1) // tmp1 = ⟨r̂ , ẑ⟩

		naiveE := ringQ.NewPoly()
		naiveR := ringQ.NewPoly()
		tmpMul := ringQ.NewPoly()

		for j := 0; j < k; j++ {
			ringQ.MulCoeffsMontgomery(eHat[j], zHat[j], tmpMul)
			ringQ.Add(naiveE, tmpMul, naiveE)

			ringQ.MulCoeffsMontgomery(rHat[j], zHat[j], tmpMul)
			ringQ.Add(naiveR, tmpMul, naiveR)
		}

		diffE := ringQ.NewPoly()
		diffR := ringQ.NewPoly()
		ringQ.Sub(naiveE, sum0, diffE) // diffE  should be 0
		ringQ.Sub(naiveR, sum1, diffR) // diffR  should be 0

		ringQ.InvNTT(diffE, diffE) // to coeff-domain for easy test
		ringQ.InvNTT(diffR, diffR)

		centre := func(v uint64) int64 { // helper already used elsewhere
			q := int64(ringQ.Modulus[0])
			x := int64(v)
			if x > q/2 {
				x -= q
			}
			return x
		}

		for t := 0; t < ringQ.N; t++ {
			if centre(diffE.Coeffs[0][t]) != 0 {
				log.Fatalf("[CHK-1] ⟨ê,ẑ⟩ mismatch at slot %d: %d",
					t, centre(diffE.Coeffs[0][t]))
			}
			if centre(diffR.Coeffs[0][t]) != 0 {
				log.Fatalf("[CHK-1] ⟨r̂,ẑ⟩ mismatch at slot %d: %d",
					t, centre(diffR.Coeffs[0][t]))
			}
		}
		fmt.Println("[CHK-1] inner products verified ✔")

		// -------------------------------------------------------------------------
		// 2)  Verify   ⟨ê,ẑ⟩ + a·⟨r̂,ẑ⟩   equals   Σ (a·r̂_j + ê_j)·ẑ_j
		// -------------------------------------------------------------------------
		// 2a) Left-hand side:  lhsEval = ⟨ê,ẑ⟩ + a·⟨r̂,ẑ⟩   (evaluation domain)
		lhsEval := ringQ.NewPoly()
		ring.Copy(sum0, lhsEval) // lhsEval ← ⟨ê,ẑ⟩
		tmpMul = ringQ.NewPoly()
		ringQ.MulCoeffsMontgomery(A[1], sum1, tmpMul) // tmpMul = a·⟨r̂,ẑ⟩
		ringQ.Add(lhsEval, tmpMul, lhsEval)

		// 2b) Right-hand side:  rhsEval = Σ (a·r̂_j + ê_j)·ẑ_j
		rhsEval := ringQ.NewPoly()
		tmp1 = ringQ.NewPoly()  // tmp1  = a·r̂_j + ê_j
		tmp2 := ringQ.NewPoly() // tmp2  = (a·r̂_j + ê_j)·ẑ_j
		for j := 0; j < k; j++ {
			ringQ.MulCoeffsMontgomery(A[1], rHat[j], tmp1) // a·r̂_j
			ringQ.Add(tmp1, eHat[j], tmp1)                 // a·r̂_j + ê_j
			ringQ.MulCoeffsMontgomery(tmp1, zHat[j], tmp2) // ⋅ ẑ_j
			ringQ.Add(rhsEval, tmp2, rhsEval)
		}

		// 2c) Compare in coefficient domain
		lhsCoeff := ringQ.NewPoly()
		rhsCoeff := ringQ.NewPoly()
		ringQ.InvNTT(lhsEval, lhsCoeff)
		ringQ.InvNTT(rhsEval, rhsCoeff)

		diff := ringQ.NewPoly()
		ringQ.Sub(lhsCoeff, rhsCoeff, diff)

		for t := 0; t < ringQ.N; t++ {
			if centre(diff.Coeffs[0][t]) != 0 {
				log.Fatalf("[CHK-2] mismatch at slot %d: %d",
					t, centre(diff.Coeffs[0][t]))
			}
		}
		fmt.Println("[CHK-2] ⟨ê,ẑ⟩ + a·⟨r̂,ẑ⟩  ≡  Σ (a·r̂_j+ê_j)·ẑ_j  ✔")
	}
	// -----------------------------------------------------------------------------
	// End of sanity-check block
	// -----------------------------------------------------------------------------

	//! ---------- DEBUG B : verify Ap + Gz = u ----------
	if true {
		// helper to centre a coeff into (−q/2, q/2]
		centre := func(v uint64) int64 {
			q := int64(ringQ.Modulus[0])
			x := int64(v)
			if x > q/2 {
				x -= q
			}
			return x
		}

		tmp := ringQ.NewPoly()

		// 1) Compute A·x in NTT, then back to coefficients
		AxEval := ringQ.NewPoly()
		for i := range A {
			ringQ.MulCoeffsMontgomery(A[i], x[i], tmp)
			ringQ.Add(AxEval, tmp, AxEval)
		}

		// 1b) also compute A·p and A·z separately for diagnostics
		ApEval := ringQ.NewPoly()
		for i := range p {
			ringQ.MulCoeffsMontgomery(A[i], p[i], tmp)
			ringQ.Add(ApEval, tmp, ApEval)
		}
		ApCoeff := ringQ.NewPoly()
		ringQ.InvNTT(ApEval, ApCoeff)

		AzEval := ringQ.NewPoly()
		// contributions from gadget columns
		for j := 0; j < k; j++ {
			ringQ.MulCoeffsMontgomery(A[j+2], zHat[j], tmp)
			ringQ.Add(AzEval, tmp, AzEval)
		}
		// contributions from x0,x1 parts
		eDotZEvalAcc := ringQ.NewPoly()
		rDotZEvalAcc := ringQ.NewPoly()
		for j := 0; j < k; j++ {
			ringQ.MulCoeffsMontgomery(eHat[j], zHat[j], tmp)
			ringQ.Add(eDotZEvalAcc, tmp, eDotZEvalAcc)
			ringQ.MulCoeffsMontgomery(rHat[j], zHat[j], tmp)
			ringQ.Add(rDotZEvalAcc, tmp, rDotZEvalAcc)
		}
		ringQ.MulCoeffsMontgomery(A[0], eDotZEvalAcc, tmp)
		ringQ.Add(AzEval, tmp, AzEval)
		ringQ.MulCoeffsMontgomery(A[1], rDotZEvalAcc, tmp)
		ringQ.Add(AzEval, tmp, AzEval)

		AzCoeff := ringQ.NewPoly()
		ringQ.InvNTT(AzEval, AzCoeff)

		// 1c) verify that A·p + G·z equals the target u
		Gz := CreateGadgetMatrix(ringQ, base, 1, k)
		GzEval := ringQ.NewPoly()
		for j := 0; j < k; j++ {
			ringQ.NTT(Gz[j], Gz[j])
			ringQ.MulCoeffsMontgomery(Gz[j], zHat[j], tmp)
			ringQ.Add(GzEval, tmp, GzEval)
		}

		GzCoeff_FromVector := ringQ.NewPoly()
		ringQ.InvNTT(GzEval, GzCoeff_FromVector)

		uCoeff := ringQ.NewPoly()
		ringQ.InvNTT(u, uCoeff)

		sumCheck := ringQ.NewPoly()
		ringQ.Add(ApCoeff, GzCoeff_FromVector, sumCheck)
		ringQ.Sub(sumCheck, uCoeff, sumCheck)

		diff0 := centre(sumCheck.Coeffs[0][0])
		if diff0 != 0 {
			log.Printf("[CHECK] A·p coeff[0]=%d", centre(ApCoeff.Coeffs[0][0]))
			log.Printf("[CHECK] G·z coeff[0]=%d", centre(GzCoeff_FromVector.Coeffs[0][0]))
			log.Printf("[CHECK] u coeff[0]=%d", centre(uCoeff.Coeffs[0][0]))
			log.Fatalf("[CHECK] A·p+G·z mismatch at slot 0: %d", diff0)
		}
		for t := 1; t < ringQ.N; t++ {
			if centre(sumCheck.Coeffs[0][t]) != 0 {
				log.Fatalf("[CHECK] A·p+G·z mismatch at slot %d: %d", t, centre(sumCheck.Coeffs[0][t]))
			}
		}
		fmt.Println("[CHECK] A·p + G·z ≡ u  ✔")

		// 1d) verify that the object per object assembled x indeed satisfies A·x = u
		AxEvalCheck := ringQ.NewPoly()
		ring.Copy(ApEval, AxEvalCheck)

		blockA2 := ringQ.NewPoly()
		ring.Copy(GzEval, blockA2)
		ringQ.MulCoeffsMontgomery(A[1], rDotZEvalAcc, tmp)
		ringQ.Sub(blockA2, tmp, blockA2)
		ringQ.Sub(blockA2, eDotZEvalAcc, blockA2)

		blockX01 := ringQ.NewPoly()
		ringQ.MulCoeffsMontgomery(A[1], rDotZEvalAcc, blockX01)
		ringQ.Add(blockX01, eDotZEvalAcc, blockX01)

		ringQ.Add(AxEvalCheck, blockA2, AxEvalCheck)
		ringQ.Add(AxEvalCheck, blockX01, AxEvalCheck)

		AxCoeff := ringQ.NewPoly()
		ringQ.InvNTT(AxEvalCheck, AxCoeff)

		diffAx := ringQ.NewPoly()
		ringQ.Sub(AxCoeff, uCoeff, diffAx)

		for t := 1; t < ringQ.N; t++ {
			if centre(diffAx.Coeffs[0][t]) != 0 {
				log.Fatalf("[CHECK] A·x mismatch at slot %d: %d", t, centre(diffAx.Coeffs[0][t]))
			}
		}
		fmt.Println("[CHECK] A·x ≡ u  ✔")

		// 6a) accumEval = A ⋅ x in evaluation domain
		accumEval := ringQ.NewPoly()
		tmpEval := ringQ.NewPoly()
		for i := 0; i < len(A); i++ {
			ringQ.MulCoeffsMontgomery(A[i], x[i], tmpEval)
			ringQ.Add(accumEval, tmpEval, accumEval)
		}

		// 6b) Verify that A·x ≡ u mod q
		for i := 0; i < ringQ.N; i++ {
			if accumEval.Coeffs[0][i] != u.Coeffs[0][i] {
				log.Fatalf("[GaussSamp] A·x mismatch at slot %d: got %d want %d",
					i, accumEval.Coeffs[0][i], u.Coeffs[0][i])
			}
		}
	}
	//! ---------- END DEBUG B ----------
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
		ringQ.MulCoeffsMontgomery(trap.A[i], x[i], tmpEval)
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
