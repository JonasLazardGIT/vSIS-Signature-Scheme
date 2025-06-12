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

	q := int64(ringQ.Modulus[0])
	out := make([]*ring.Poly, k)

	for row := 0; row < k; row++ {
		if len(Z[row]) != N {
			log.Fatalf("row %d length %d ≠ %d", row, len(Z[row]), N)
		}

		p := ringQ.NewPoly() // coefficient domain
		for t := 0; t < N; t++ {
			v := Z[row][t] % q // canonical 0…q-1
			if v < 0 {
				v += q
			}
			p.Coeffs[0][t] = uint64(v)
		}
		ringQ.NTT(p, p) // to evaluation domain
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

	if true {
		coeff := ringQ.NewPoly()
		ringQ.InvNTT(zHat[0], coeff)
		fmt.Printf("[DEBUG] zHat[0](0) = %d (expected %d)\n", UnsignedToSigned(coeff.Coeffs[0][0], ringQ.Modulus[0]), Zmat[0][0])
	}

	// Precompute the negacyclic transposes of the trapdoor polynomials for
	// assembling x. SamplePz uses the originals and applies the transpose
	// internally when computing its dot-products.
	rHatT := make([]*ring.Poly, k)
	eHatT := make([]*ring.Poly, k)
	for j := 0; j < k; j++ {
		rTmp := ringQ.NewPoly()
		ringQ.InvNTT(rHat[j], rTmp)
		rTmp = AutomorphismTranspose(ringQ, rTmp)
		ringQ.NTT(rTmp, rTmp)
		rHatT[j] = rTmp

		eTmp := ringQ.NewPoly()
		ringQ.InvNTT(eHat[j], eTmp)
		eTmp = AutomorphismTranspose(ringQ, eTmp)
		ringQ.NTT(eTmp, eTmp)
		eHatT[j] = eTmp
	}

	fmt.Println("GaussSamp: using transposed trapdoor for inner products")

	// 5) assemble x = [ p₀ + ê⋅zHat,  p₁ + r̂⋅zHat,  p₂+ẑ₀, …, p_{k+1}+ẑ_{k-1} ]
	x := make([]*ring.Poly, k+2)

	// row 0: p[0] + <eHatᵀ, zHat>
	sum0 := ringQ.NewPoly()
	tmpez := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		ringQ.MulCoeffsMontgomery(eHatT[j], zHat[j], tmpez)
		ringQ.Add(sum0, tmpez, sum0)
	}
	x[0] = ringQ.NewPoly()
	// In the PALISADE C++ implementation these terms are added
	// when forming the first two rows of x.
	ringQ.Add(p[0], sum0, x[0])

	// row 1: p[1] + <rHatᵀ, zHat>
	sum1 := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		ringQ.MulCoeffsMontgomery(rHatT[j], zHat[j], tmpez)
		ringQ.Add(sum1, tmpez, sum1)
	}
	x[1] = ringQ.NewPoly()
	// PALISADE also adds here
	ringQ.Add(p[1], sum1, x[1])

	// rows 2…k+1: just p[i] + zHat[i-2]
	for i := 2; i < k+2; i++ {
		x[i] = ringQ.NewPoly()
		ringQ.Add(p[i], zHat[i-2], x[i])
	}

	//! ---------- DEBUG B : verify x₀ formula and A·x = u ----------
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
		AxCoeff := ringQ.NewPoly()
		ringQ.InvNTT(AxEval, AxCoeff)

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
			ringQ.MulCoeffsMontgomery(eHatT[j], zHat[j], tmp)
			ringQ.Add(eDotZEvalAcc, tmp, eDotZEvalAcc)
			ringQ.MulCoeffsMontgomery(rHatT[j], zHat[j], tmp)
			ringQ.Add(rDotZEvalAcc, tmp, rDotZEvalAcc)
		}
		ringQ.MulCoeffsMontgomery(A[0], eDotZEvalAcc, tmp)
		ringQ.Add(AzEval, tmp, AzEval)
		ringQ.MulCoeffsMontgomery(A[1], rDotZEvalAcc, tmp)
		ringQ.Add(AzEval, tmp, AzEval)

		AzCoeff := ringQ.NewPoly()
		ringQ.InvNTT(AzEval, AzCoeff)

		fmt.Printf("[CHECK] A·p coeff[0]=%d  A·z coeff[0]=%d\n", UnsignedToSigned(ApCoeff.Coeffs[0][0], ringQ.Modulus[0]), UnsignedToSigned(AzCoeff.Coeffs[0][0], ringQ.Modulus[0]))

		// 1c) verify that A·p + G·z equals the target u
		Gz := CreateGadgetMatrix(ringQ, base, 1, k)
		GzEval := ringQ.NewPoly()
		for j := 0; j < k; j++ {
			ringQ.NTT(Gz[j], Gz[j])
			ringQ.MulCoeffsMontgomery(Gz[j], zHat[j], tmp)
			ringQ.Add(GzEval, tmp, GzEval)
		}
		GzCoeff_FromMat := ringQ.NewPoly()
		for j := 0; j < ringQ.N; j++ {
			recomb := int64(0)
			for idx := k - 1; idx >= 0; idx-- {
				recomb = (recomb*int64(base) + Zmat[idx][j]) % int64(ringQ.Modulus[0])
			}
			GzCoeff_FromMat.Coeffs[0][j] = SignedToUnsigned(recomb, ringQ.Modulus[0])
		}

		GzCoeff_FromVector := ringQ.NewPoly()
		ringQ.InvNTT(GzEval, GzCoeff_FromVector)
		for t := 0; t < 10; t++ {
			fmt.Printf(("G.Z fromMat coeff[%d]=%d  G.Z FromVector coeff[%d]=%d\n"), t, (GzCoeff_FromMat.Coeffs[0][t]), t, (GzCoeff_FromVector.Coeffs[0][t]))
		}

		uCoeff := ringQ.NewPoly()
		ringQ.InvNTT(u, uCoeff)

		sumCheck := ringQ.NewPoly()
		ringQ.Add(ApCoeff, GzCoeff_FromMat, sumCheck)
		ringQ.Sub(sumCheck, uCoeff, sumCheck)

		diff0 := centre(sumCheck.Coeffs[0][0])
		if diff0 != 0 {
			log.Printf("[CHECK] A·p coeff[0]=%d", centre(ApCoeff.Coeffs[0][0]))
			log.Printf("[CHECK] G·z coeff[0]=%d", centre(GzCoeff_FromMat.Coeffs[0][0]))
			log.Printf("[CHECK] u coeff[0]=%d", centre(uCoeff.Coeffs[0][0]))
			log.Fatalf("[CHECK] A·p+G·z mismatch at slot 0: %d", diff0)
		}
		for t := 1; t < ringQ.N; t++ {
			if centre(sumCheck.Coeffs[0][t]) != 0 {
				log.Fatalf("[CHECK] A·p+G·z mismatch at slot %d: %d", t, centre(sumCheck.Coeffs[0][t]))
			}
		}
		fmt.Println("[CHECK] A·p + G·z ≡ u  ✔")

		// 3) Inverse‐NTT p₀ and extract constant
		p0CoeffPoly := ringQ.NewPoly()
		ringQ.InvNTT(p[0], p0CoeffPoly)

		// 4) Inverse‐NTT x₀
		xCoeffPoly := ringQ.NewPoly()
		ringQ.InvNTT(x[0], xCoeffPoly)

		// 5) Compute e·z in NTT, then back to coeffs
		eDotZEval := ringQ.NewPoly()
		for j := 0; j < k; j++ {
			ringQ.MulCoeffsMontgomery(eHatT[j], zHat[j], tmp)
			ringQ.Add(eDotZEval, tmp, eDotZEval)
		}
		eDotZCoeff := ringQ.NewPoly()
		ringQ.InvNTT(eDotZEval, eDotZCoeff)
		// coefficient-wise cancel check for x0
		for t := 0; t < ringQ.N; t++ {
			p0t := centre(p0CoeffPoly.Coeffs[0][t])
			ez := centre(eDotZCoeff.Coeffs[0][t])
			xt := centre(xCoeffPoly.Coeffs[0][t])
			if xt < 0 {
				xt += int64(ringQ.Modulus[0])
			}
			p0tplusEz := p0t + ez
			if p0tplusEz < 0 && xt >= 0 {
				p0tplusEz += int64(ringQ.Modulus[0])
			}
			if p0tplusEz != xt {
				log.Fatalf("Cancel-check x0 mismatch slot %d: p0=%d e·z=%d sum=%d x0=%d", t, p0t, ez, p0tplusEz, xt)
			}
		}

		p1CoeffPoly := ringQ.NewPoly()
		ringQ.InvNTT(p[1], p1CoeffPoly)
		x1CoeffPoly := ringQ.NewPoly()
		ringQ.InvNTT(x[1], x1CoeffPoly)
		rDotZEval := ringQ.NewPoly()
		for j := 0; j < k; j++ {
			ringQ.MulCoeffsMontgomery(rHatT[j], zHat[j], tmp)
			ringQ.Add(rDotZEval, tmp, rDotZEval)
		}
		rDotZCoeff := ringQ.NewPoly()
		ringQ.InvNTT(rDotZEval, rDotZCoeff)
		for t := 0; t < ringQ.N; t++ {
			p1t := centre(p1CoeffPoly.Coeffs[0][t])
			rz := centre(rDotZCoeff.Coeffs[0][t])
			xt := centre(x1CoeffPoly.Coeffs[0][t])
			p1tplusrz := p1t + rz
			if p1tplusrz < 0 && xt >= 0 {
				p1tplusrz += int64(ringQ.Modulus[0])
			}
			if p1tplusrz > 0 && xt < 0 {
				p1tplusrz -= int64(ringQ.Modulus[0])
			}
			if p1tplusrz != xt {
				log.Fatalf("Cancel-check x1 mismatch slot %d: p1=%d r·z=%d sum=%d x1=%d", t, p1t, rz, p1tplusrz, xt)
			}
		}

		//! ---------- DEBUG C : verify gadget rows ----------
		G := CreateGadgetMatrix(ringQ, base, 1, k)
		for j := 0; j < k; j++ {
			ringQ.NTT(G[j], G[j])
		}

		for j := 0; j < k; j++ {
			evalGad := ringQ.NewPoly()
			evalGp := ringQ.NewPoly()
			evalGz := ringQ.NewPoly()

			ringQ.MulCoeffsMontgomery(A[j+2], x[j+2], evalGad)
			ringQ.MulCoeffsMontgomery(G[j], p[j+2], evalGp)
			ringQ.MulCoeffsMontgomery(G[j], zHat[j], evalGz)

			// (a*rHat_j + eHat_j) * x_{j+2}
			trapTerm := ringQ.NewPoly()
			ringQ.MulCoeffsMontgomery(A[1], rHat[j], trapTerm)
			ringQ.Add(trapTerm, eHat[j], trapTerm)
			ringQ.MulCoeffsMontgomery(trapTerm, x[j+2], trapTerm)

			coeffGad := ringQ.NewPoly()
			coeffGp := ringQ.NewPoly()
			coeffGz := ringQ.NewPoly()
			coeffTrap := ringQ.NewPoly()

			ringQ.InvNTT(evalGad, coeffGad)
			ringQ.InvNTT(evalGp, coeffGp)
			ringQ.InvNTT(evalGz, coeffGz)
			ringQ.InvNTT(trapTerm, coeffTrap)

			rhs := ringQ.NewPoly()
			ringQ.Add(coeffGp, coeffGz, rhs)
			ringQ.Sub(rhs, coeffTrap, rhs)

			diff := ringQ.NewPoly()
			ringQ.Sub(coeffGad, rhs, diff)

			for t := 0; t < ringQ.N; t++ {
				if centre(diff.Coeffs[0][t]) != 0 {
					log.Fatalf("Gad-row %d, slot %d: diff = %d ≠ 0", j, t, centre(diff.Coeffs[0][t]))
				}
			}
			fmt.Printf("Gad-row %d: all %d coefficients matched\n", j, ringQ.N)
		}
		//! ------- END DEBUG C -------
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

	// --- trapdoor self-check (Alg 1 correctness) ---
	G := CreateGadgetMatrix(ringQ, trap.Base, trap.Rows, trap.K)
	for j := 0; j < trap.K; j++ {
		// send G[j] to EVAL
		ringQ.NTT(G[j], G[j])

		// tmp = a * r̂_j + ê_j
		tmp := ringQ.NewPoly()
		a := trap.A1[1]    // the "a" sample
		rj := trap.R[0][j] // r̂_j
		ej := trap.R[1][j] // ê_j
		ringQ.MulCoeffsMontgomery(a, rj, tmp)
		ringQ.Add(tmp, ej, tmp)

		// sum = tmp + A₂[j]
		sum := ringQ.NewPoly()
		ringQ.Add(tmp, trap.A2[j], sum)

		// check sum == G[j] coefficient-wise
		for t := 0; t < ringQ.N; t++ {
			want := G[j].Coeffs[0][t]
			got := sum.Coeffs[0][t]
			if want != got {
				log.Fatalf("TrapGen self-check failed at row %d, slot %d: got %d want %d", j, t, got, want)
			}
		}
	}
	fmt.Println("✔ Trapdoor generation self-check passed")

	// ----------------------------------------------------------------------------
	// 4) Sample a random syndrome u ∈ R_q in EVAL
	// ----------------------------------------------------------------------------
	u := ringQ.NewPoly()
	for i := 0; i < ringQ.N; i++ {
		u.Coeffs[0][i] = uint64(i % int(qMod))
	}
	uEval := ringQ.NewPoly() // u in evaluation domain
	ringQ.NTT(u, uEval)

	// keep a copy of uEval before any mutation
	uEvalCopy := ringQ.NewPoly()
	ring.Copy(uEval, uEvalCopy)

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
	for i, xi := range x {
		ringQ.MulCoeffsMontgomery(trap.A[i], xi, tmpEval)
		ringQ.Add(accumEval, tmpEval, accumEval)
	}

	// 6b) bring both accumEval and uEvalCopy back to coefficient domain
	accumCoeff := ringQ.NewPoly()
	uCoeffOrig := ringQ.NewPoly()
	ringQ.InvNTT(accumEval, accumCoeff)
	ringQ.InvNTT(uEvalCopy, uCoeffOrig)
	fmt.Printf("accumCoeff.Coeffs[0][:8]=%v\n", accumCoeff.Coeffs[0][:8])
	fmt.Printf("uCoeffOrig.Coeffs[0][:8]=%v\n", uCoeffOrig.Coeffs[0][:8])

	// 6c) form difference diff = accumCoeff - uCoeffOrig in Z_q
	diff := ringQ.NewPoly()
	ringQ.Sub(accumCoeff, uCoeffOrig, diff)

	// 6d) map to signed integers and assert all zero
	for i, coeff := range diff.Coeffs[0] {
		signed := UnsignedToSigned(coeff, ringQ.Modulus[0])
		if signed != 0 {
			log.Fatalf("FAIL: verification mismatch at slot %d: ≡ %d mod q", i, signed)
		}
	}

	fmt.Println("SUCCESS: exact A·x ≡ u mod q")
}
