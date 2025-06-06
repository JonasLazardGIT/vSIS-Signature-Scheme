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

// calculateParams computes σₜ=(t+1)σ and
func CalculateParams(base uint64, n, k int) (sigmaT, s float64) {
	sigma := 3.19                    // smoothing parameter from Sec V-A1
	sigmaT = float64(base+1) * sigma // σₜ = (t+1)·σ
	s = SpectralBound(n, k, base)    // 1.3 * sigma * sigmaT * sqrt(n*k + 2*n + 4.7)
	return
}

// ZtoZhat converts an integer matrix Z ∈ ℤ^{κ×N} into polys in R_q^κ (full CRT).
func ZtoZhat(Z [][]int64, ringQ *ring.Ring) []*ring.Poly {
	κ, N := len(Z), ringQ.N
	out := make([]*ring.Poly, κ)
	for i := 0; i < κ; i++ {
		P := ringQ.NewPoly()
		for lvl, qi := range ringQ.Modulus {
			mod := int64(qi)
			for t := 0; t < N; t++ {
				c := SignedToUnsigned(Z[i][t], uint64(mod))
				P.Coeffs[lvl][t] = c
			}
		}
		ringQ.NTT(P, P)
		// fmt.Print("P = ")
		// for t := 0; t < N; t++ {
		// 	fmt.Print(P.Coeffs[0][t], " ")
		// }
		out[i] = P
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
	p := SamplePz(ringQ, s, (2+1)*sigma, [2][]*ring.Poly{rHat, eHat}, k+int(base), 256)

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
	dgg := NewDiscreteGaussian(sigma) // discrete Gaussian sampler for G-sampling
	// fmt.Printf("uCoeffs = %v\n", uCoeffs)
	Zmat := SampleGDiscrete(ringQ, (float64(base+1) * sigma), base, uCoeffs, k, dgg)

	// Zmat is a matrix of integers Z ∈ ℤ^{κ×N} where each row is a poly in R_q

	// 4) CRT+NTT each row of Z to get zHat ∈ R_q^κ
	zHat := ZtoZhat(Zmat, ringQ)

	// 5) assemble x = [ p₀ + ê⋅zHat,  p₁ + r̂⋅zHat,  p₂+ẑ₀, …, p_{k+1}+ẑ_{k-1} ]
	x := make([]*ring.Poly, k+2)

	// row 0: p[0] + <eHat, zHat>
	sum0 := ringQ.NewPoly()
	tmpez := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		ringQ.MulCoeffsMontgomery(eHat[j], zHat[j], tmpez)
		ringQ.Add(sum0, tmpez, sum0)
	}
	x[0] = ringQ.NewPoly()
	ringQ.Add(p[0], sum0, x[0])

	// row 1: p[1] + <rHat, zHat>
	sum1 := ringQ.NewPoly()
	for j := 0; j < k; j++ {
		ringQ.MulCoeffsMontgomery(rHat[j], zHat[j], tmpez)
		ringQ.Add(sum1, tmpez, sum1)
	}
	x[1] = ringQ.NewPoly()
	ringQ.Add(p[1], sum1, x[1])

	// rows 2…k+1: just p[i] + zHat[i-2]
	for i := 2; i < k+2; i++ {
		x[i] = ringQ.NewPoly()
		ringQ.Add(p[i], zHat[i-2], x[i])
	}

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

		// check sum == G[j]
		if !ringQ.Equal(sum, G[j]) {
			panic(fmt.Sprintf("TrapGen self-check failed at j=%d", j))
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
