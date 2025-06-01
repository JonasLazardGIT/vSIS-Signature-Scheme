package Preimage_Sampler

import (
	"log"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Trapdoor holds the public key A and the secret trapdoor (r̂,ê).
// Base = gadget base t, K = gadget dimension κ.
type Trapdoor struct {
	A    []*ring.Poly    // full public row: [1, a, g₁−(a·r̂₁+ê₁), …, g_κ−(a·r̂_κ+ê_κ)]
	A1   []*ring.Poly    // first block [a, ]
	A2   []*ring.Poly    // gadget block [gᵢ−(…)]
	R    [2][]*ring.Poly // secret trapdoor rows: R[0]=r̂, R[1]=ê
	Base uint64          // gadget base t
	K    int             // gadget length κ
	Rows int             // =1
	Cols int             // =2
}

// CreateGadgetMatrix returns the gadget matrix G of size (rows × k).
// In particular, for rows=1 it returns the row-vector
//
//	[ g_1, g_2, …, g_k ]  where  g_j = base^(j-1) mod q
//
// as constant polynomials in Rq.
// We allocate rows*k polys, in row-major order; unused for rows>1 in this trapdoor scheme.
//
// ringQ : the R_q polynomial ring
// base  : the gadget base t
// rows  : number of gadget rows (for our TrapGen, rows == 1)
// k     : gadget length κ
func CreateGadgetMatrix(ringQ *ring.Ring, base uint64, rows, k int) []*ring.Poly {
	N := ringQ.N
	// total number of polys = rows×k
	G := make([]*ring.Poly, rows*k)

	// For each row i and gadget index j, set G[i*k + j] = constant poly t^(j)
	for i := 0; i < rows; i++ {
		for j := 0; j < k; j++ {
			idx := i*k + j
			p := ringQ.NewPoly()

			// compute t^j mod each prime in the CRT chain
			for tIdx, qi := range ringQ.Modulus {
				mod := uint64(qi)
				// fast exponentiation of base^j mod mod
				var power uint64 = 1
				for e := 0; e < j; e++ {
					power = (power * base) % mod
				}
				// fill *all* N coefficients with this constant
				for coeffIdx := 0; coeffIdx < N; coeffIdx++ {
					p.Coeffs[tIdx][coeffIdx] = power
				}
			}

			G[idx] = p
		}
	}
	return G
}

// TrapGen implements Algorithm 1 (Ring-LWE trapdoor).
//   - ringQ: the R_q ring
//   - base t, sigmaT: Gaussian width σₜ for sampling r̂, ê
func TrapGen(ringQ *ring.Ring, base uint64, sigmaT float64) Trapdoor {
	// 0) Compute κ = ceil(log_t q)
	q := float64(ringQ.Modulus[0])
	k := int(math.Ceil(math.Log(q) / math.Log(float64(base))))

	// 1) sample a ← Uniform(R_q)
	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	uniform := ring.NewUniformSampler(prng, ringQ)
	a := ringQ.NewPoly()
	uniform.Read(a)
	ringQ.NTT(a, a)

	// 2) sample r̂, ê ∈ R_q^κ from D_{R,σₜ} (discrete Gaussian)
	Rhat := make([]*ring.Poly, k)
	Ehat := make([]*ring.Poly, k)
	dg := NewDiscreteGaussian(sigmaT)
	for j := 0; j < k; j++ {
		// draw each coefficient exactly via Karney
		Rhat[j] = ringQ.NewPoly()
		Ehat[j] = ringQ.NewPoly()
		for lvl, qi := range ringQ.Modulus {
			mod := int64(qi)
			for i := 0; i < ringQ.N; i++ {
				// sample one integer, reduce mod qi
				r := dg.Draw(0)
				e := dg.Draw(0)
				Rhat[j].Coeffs[lvl][i] = uint64((r%mod + mod) % mod)
				Ehat[j].Coeffs[lvl][i] = uint64((e%mod + mod) % mod)
			}
		}
		ringQ.NTT(Rhat[j], Rhat[j])
		ringQ.NTT(Ehat[j], Ehat[j])
	}

	// 3) build gadget row g₁,…,g_κ as constant polys
	G := CreateGadgetMatrix(ringQ, base, 1, k)
	for j := range G {
		ringQ.NTT(G[j], G[j])
	}
	// 4) build A₂[j] = g_j − (a⋅r̂_j + ê_j)
	A2 := make([]*ring.Poly, k)
	for j := 0; j < k; j++ {
		// assume a, Rhat[j], Ehat[j], G[j] are all in EVAL
		tmp := ringQ.NewPoly()
		ringQ.MulCoeffsMontgomery(a, Rhat[j], tmp)
		ringQ.Add(tmp, Ehat[j], tmp)
		A2[j] = ringQ.NewPoly()
		ringQ.Sub(G[j], tmp, A2[j])
	}

	// 5) build A₁ = [1, a]
	one := ringQ.NewPoly()
	for lvl := range ringQ.Modulus {
		one.Coeffs[lvl][0] = 1
	}
	ringQ.NTT(one, one) // send to EVALUATION
	A1 := []*ring.Poly{one, a}

	// 6) flatten A = [ A₁ ∥ A₂ ]
	A := append(A1, A2...)

	return Trapdoor{
		A:    A,
		A1:   A1,
		A2:   A2,
		R:    [2][]*ring.Poly{Rhat, Ehat},
		Base: base,
		K:    k,
		Rows: 1,
		Cols: 2,
	}
}

// Perturb implements PERTURB(σ,ℓ,h,base) from Alg 3 of the paper.
// It uses the new discrete–Gaussian sampler so that each zi ← Dℤ(mean,σi).
func Perturb(sigma float64, ell, h []float64, base uint64) []int64 {
	k := len(ell)
	// 1) sample z-vector via exact Dℤ
	z := make([]int64, k)
	beta := 0.0
	for i := 0; i < k; i++ {
		mean := beta / ell[i]
		sigmaI := sigma / ell[i]
		dg := NewDiscreteGaussian(sigmaI)
		z[i] = dg.Draw(mean)
		beta = -float64(z[i]) * h[i]
	}

	// 2) build the gadget–base combination p
	p := make([]int64, k)

	// first coordinate: (2·base+1)*z0 + base·z1
	p[0] = int64(2*base+1)*z[0] + int64(base)*z[1]

	// middle coordinates: base*(z[i−1] + 2·z[i] + z[i+1])
	for i := 1; i < k-1; i++ {
		sum := z[i-1] + 2*z[i] + z[i+1]
		p[i] = int64(base) * sum
	}

	// last coordinate: base*(z[k-2] + 2·z[k-1])
	p[k-1] = int64(base) * (z[k-2] + 2*z[k-1])

	return p
}

// SampleD implements the discrete‐Sample_D of Alg.3.
//
//	sigma: the stddev σ
//	  a   : the accumulator vector a ∈ ℝ^k (we modify it in place)
//	  c   : the conditioning vector c ∈ ℝ^k
//
// returns z ∈ ℤ^k
func SampleD(sigma float64, a, c []float64) []int {
	k := len(a)
	z := make([]int, k)

	// 1) last coord
	mean := -a[k-1] / c[k-1]
	stdp := sigma / c[k-1]
	dg0 := NewDiscreteGaussian(stdp)
	z[k-1] = int(dg0.Draw(mean))

	// update a ← a + zₖ₋₁·c
	for i := range a {
		a[i] += float64(z[k-1]) * c[i]
	}

	// 2) remaining coords
	for i := 0; i < k-1; i++ {
		dg := NewDiscreteGaussian(sigma)
		z[i] = int(dg.Draw(-a[i]))
	}

	return z
}

// SampleGDiscrete runs the *discrete* G-sampling (Alg 3) exactly.
//
//	ringQ   : R_q ring
//	sigma   : continuous Gaussian width σₜ
//	base    : gadget base t
//	uCoeff  : coefficient form of u(x) in [0,q), len=N=ringQ.N
//	k       : gadget length (digits)
//
// returns a k×N matrix Z of ints.
func SampleGDiscrete(
	ringQ *ring.Ring,
	sigma float64,
	base uint64,
	uCoeff []uint64,
	k int,
) [][]int64 {
	const bigPrec = 256 // précision binaire pour tous les big.Float
	N := ringQ.N

	// 1) build ell[], h[] en float64 (pas critique pour la précision des centres)
	ell := make([]float64, k)
	h := make([]float64, k)
	ell[0] = math.Sqrt(float64(base)*(1+1/float64(k)) + 1)
	for i := 1; i < k; i++ {
		ell[i] = math.Sqrt(float64(base) * (1 + 1/float64(k-i)))
	}
	h[0] = 0
	for i := 1; i < k; i++ {
		h[i] = math.Sqrt(float64(base) * (1 - 1/float64(k-(i-1))))
	}

	// 2) build dBig[] en big.Float à partir des digits de q
	q := ringQ.Modulus[0]
	modDigits := baseDigits(int64(q), int64(base), k)

	// dBig[i] = (modDigits[0] + modDigits[1] + ... + modDigits[i]) / 2^i
	dBig := make([]*big.Float, k)
	// On commence avec dBig[0] = modDigits[0] / 2
	twoB := new(big.Float).SetPrec(bigPrec).SetFloat64(2.0)
	dBig[0] = new(big.Float).SetPrec(bigPrec).
		Quo(
			new(big.Float).SetPrec(bigPrec).SetFloat64(float64(modDigits[0])),
			twoB,
		)
	// puis dBig[i] = (dBig[i-1] + modDigits[i]) / 2
	for i := 1; i < k; i++ {
		tmp := new(big.Float).SetPrec(bigPrec).
			Add(
				dBig[i-1],
				new(big.Float).SetPrec(bigPrec).SetFloat64(float64(modDigits[i])),
			)
		dBig[i] = new(big.Float).SetPrec(bigPrec).Quo(tmp, twoB)
	}

	// 3) sigma′ = σₜ/(base+1)
	sigmaP := sigma / float64(base+1)

	// 4) préparer la sortie Z[k][N]
	Z := make([][]int64, k)
	for i := range Z {
		Z[i] = make([]int64, N)
	}

	// 5) boucle par coefficient j = 0..N-1
	for j := 0; j < N; j++ {
		v := uCoeff[j]                                  // u_j ∈ [0,q)
		vDigits := baseDigits(int64(v), int64(base), k) // décomposition base‐t en ints

		// 5a) perturb → p[0..k-1] (on ignore la précision ici, c’est du perturbe exact)
		p := Perturb(sigmaP, ell, h, base)

		// 5b) construire les centres cBig[0..k-1] en big.Float
		cBig := make([]*big.Float, k)

		// cBig[0] = (vDigits[0] - p[0]) / base
		numer0 := new(big.Float).SetPrec(bigPrec).SetFloat64(
			float64(vDigits[0] - p[0]),
		)
		cBig[0] = new(big.Float).SetPrec(bigPrec).Quo(numer0, twoB)

		// Pour i=1..k-1 : cBig[i] = (cBig[i-1] + vDigits[i] - p[i]) / base
		for i := 1; i < k; i++ {
			tmp := new(big.Float).SetPrec(bigPrec).
				Add(
					cBig[i-1],
					new(big.Float).SetPrec(bigPrec).SetFloat64(float64(vDigits[i]-p[i])),
				)
			cBig[i] = new(big.Float).SetPrec(bigPrec).Quo(tmp, twoB)

			// Vérification de l’invariant carry en multiprécision :
			//   base*c[i] ?= c[i-1] + vDigits[i] - p[i]
			lhs := new(big.Float).SetPrec(bigPrec).
				Mul(new(big.Float).SetPrec(bigPrec).SetFloat64(float64(base)), cBig[i])
			rhs := new(big.Float).SetPrec(bigPrec).
				Add(
					cBig[i-1],
					new(big.Float).SetPrec(bigPrec).SetFloat64(float64(vDigits[i]-p[i])),
				)
			diff := new(big.Float).SetPrec(bigPrec).Sub(lhs, rhs)
			absDiff := new(big.Float).SetPrec(bigPrec).Abs(diff)
			tolBig := new(big.Float).SetPrec(bigPrec).SetFloat64(1e-30)
			if absDiff.Cmp(tolBig) > 0 {
				lhsF, _ := lhs.Float64()
				rhsF, _ := rhs.Float64()
				log.Printf(
					"CARRY_ERR (j=%d,i=%d): base*c=%.12g, (c_prev+v-p)=%.12g\n",
					j, i, lhsF, rhsF,
				)
			}
		}

		// 5c) construire d[] en float64 à partir de dBig (pour le SampleD)
		d := make([]float64, k)
		for i := 0; i < k; i++ {
			d[i], _ = dBig[i].Float64() // conversion unique ici
		}

		// 5d) bâtir le slice “centers” en float64 pour SampleD :
		c := make([]float64, k)
		for i := 0; i < k; i++ {
			c[i], _ = cBig[i].Float64() // conversion unique ici
		}

		// 5e) discrete‐Gaussian sample z[0..k-1]
		z := SampleD(sigmaP, c, d)

		// 5f) reconstituer t = (t₀,...,t_{k-1}) en int64
		//    t₀ = 2 z₀ + q₀ z_{k-1} + v₀
		Z[0][j] = 2*int64(z[0]) +
			int64(modDigits[0])*int64(z[k-1]) +
			int64(vDigits[0])

		//    t_i = 2 z_i − z_{i-1} + q_i z_{k-1} + v_i  (pour i=1..k-2)
		for i := 1; i < k-1; i++ {
			Z[i][j] = 2*int64(z[i]) -
				int64(z[i-1]) +
				int64(modDigits[i])*int64(z[k-1]) +
				int64(vDigits[i])
		}

		//    t_{k-1} = q_{k-1} z_{k-1} − z_{k-2} + v_{k-1}
		Z[k-1][j] = int64(modDigits[k-1])*int64(z[k-1]) -
			int64(z[k-2]) +
			int64(vDigits[k-1])
	}

	return Z
}
