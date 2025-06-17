package zkproof

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// QuadraticGate holds the (R2,R1,R0) polynomials defining a quadratic constraint.
type QuadraticGate struct {
	R2 [][]*ring.Poly
	R1 []*ring.Poly
	R0 *ring.Poly
}

// BuildQuadraticGate constructs a QuadraticGate following Fig. 6 of
// "vSIS-BBS: Blind Signatures from Lattices". All polynomials are kept
// in the NTT domain.
func BuildQuadraticGate(ringQ *ring.Ring, b1 []*ring.Poly, A [][]*ring.Poly, B0 [][]*ring.Poly, m, lm, lr int) *QuadraticGate {
	fmt.Println("BuildQuadraticGate: start")
	n := len(A)
	prng, _ := utils.NewPRNG()
	uni := ring.NewUniformSampler(prng, ringQ)

	mu := make([]*ring.Poly, n)
	for i := range mu {
		mu[i] = ringQ.NewPoly()
		uni.Read(mu[i])
		ringQ.NTT(mu[i], mu[i])
	}

	// α = μᵀ · (b1 ⊙ A)
	alpha := make([]*ring.Poly, m)
	tmp1 := ringQ.NewPoly()
	tmp2 := ringQ.NewPoly()
	for j := 0; j < m; j++ {
		acc := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			ringQ.MulCoeffs(b1[i], A[i][j], tmp1)
			ringQ.MulCoeffs(mu[i], tmp1, tmp2)
			ringQ.Add(acc, tmp2, acc)
		}
		alpha[j] = acc
	}

	// γ = μᵀ · A
	gamma := make([]*ring.Poly, m)
	for j := 0; j < m; j++ {
		acc := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			ringQ.MulCoeffs(mu[i], A[i][j], tmp1)
			ringQ.Add(acc, tmp1, acc)
		}
		gamma[j] = acc
	}

	// (δ0, δU) = μᵀ · B0
	colsB0 := len(B0[0])
	delta0 := ringQ.NewPoly()
	deltaU := make([]*ring.Poly, colsB0-1)
	for j := 0; j < colsB0; j++ {
		acc := ringQ.NewPoly()
		for i := 0; i < n; i++ {
			ringQ.MulCoeffs(mu[i], B0[i][j], tmp1)
			ringQ.Add(acc, tmp1, acc)
		}
		if j == 0 {
			delta0 = acc
		} else {
			deltaU[j-1] = ringQ.NewPoly()
			ring.Copy(acc, deltaU[j-1])
		}
	}

	L := m + lm + lr
	R2 := make([][]*ring.Poly, L+1)
	for i := range R2 {
		R2[i] = make([]*ring.Poly, L+1)
		for j := 0; j < L+1; j++ {
			R2[i][j] = ringQ.NewPoly()
		}
	}

	// set R2[j][L] and R2[L][j]
	for j := 0; j < m; j++ {
		tmp := ringQ.NewPoly()
		tmpJ := ringQ.NewPoly()
		for lvl, qi := range ringQ.Modulus {
			half := new(big.Int).SetUint64((qi + 1) >> 1)
			ringQ.MulScalarBigintLvl(lvl, gamma[j], half, tmp)
			copy(tmpJ.Coeffs[lvl], tmp.Coeffs[lvl])
		}
		ringQ.Neg(tmpJ, tmpJ)
		R2[j][L] = tmpJ
		R2[L][j] = tmpJ
	}

	// r1 = α ‖ (-δU) ‖ 0
	r1 := make([]*ring.Poly, L+1)
	idx := 0
	for _, p := range alpha {
		r1[idx] = p
		idx++
	}
	for _, p := range deltaU {
		neg := ringQ.NewPoly()
		ringQ.Neg(p, neg)
		r1[idx] = neg
		idx++
	}
	r1[idx] = ringQ.NewPoly() // zero

	// r0 = -δ0
	r0 := ringQ.NewPoly()
	ringQ.Neg(delta0, r0)
	fmt.Println("BuildQuadraticGate: done")
	return &QuadraticGate{R2: R2, R1: r1, R0: r0}
}
