package Preimage_Sampler

import (
	"math"

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
//	[ g_1, g_2, …, g_k ]  where  g_j = base^(j) mod q
//
// as constant polynomials in Rq.
// We allocate rows*k polys, in row-major order; unused for rows>1 in this trapdoor scheme.
//
// ringQ : the R_q polynomial ring
// base  : the gadget base t
// rows  : number of gadget rows (for our TrapGen, rows == 1)
// k     : gadget length κ
func CreateGadgetMatrix(ringQ *ring.Ring, base uint64, rows, k int) []*ring.Poly {
	// total number of polys = rows×k
	G := make([]*ring.Poly, rows*k)

	// For each row i and gadget index j, set G[i*k + j] = constant poly t^(j-1)
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
				// set only the constant coefficient
				p.Coeffs[tIdx][0] = power
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
		ringQ.MulCoeffs(a, Rhat[j], tmp)
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

// Perturb mirrors PALISADE's discrete Perturb (Fig. 2, 2017/308).
//
//	l, h : Cholesky scalars (length k)
//	base : gadget radix  t  (≥2)
//
// RETURNS p ∈ ℤᵏ  and keeps all intermediate computation in int64.
func Perturb(
	sigma float64,
	l, h []float64,
	base uint64,
) []int64 {

	k := len(l)
	z := make([]int64, k)

	var d float64
	for i := 0; i < k; i++ {
		dgg := NewDiscreteGaussian(sigma / l[i]) // reinitialize dgg for each i
		z[i] = dgg.Draw(d / l[i])
		d = -float64(z[i]) * h[i]
	}

	p := make([]int64, k)
	if k == 1 {
		// degenerate case: p₀ = (2·base+1) z₀
		p[0] = int64(2*base+1) * z[0]
		return p
	}

	p[0] = int64(2*base+1)*z[0] + int64(base)*z[1]
	for i := 1; i < k-1; i++ {
		p[i] = int64(base) * (z[i-1] + 2*z[i] + z[i+1])
	}
	p[k-1] = int64(base) * (z[k-2] + 2*z[k-1])
	return p
}

// SampleC reproduces PALISADE's SampleC (Figure 2, 2017/308).
//
//	c     :  constant divisor vector  (length k)
//	a     :  *mutable* accumulator  (length k) – will be updated in place
//	sigma :  continuous σ'
//
// It returns the lattice vector z (length k)   INT64.
func SampleC(
	c []float64,
	sigma float64,
	a []float64,
) []int64 {

	k := len(c)
	z := make([]int64, k)

	// 1) draw last coordinate with conditional parameters
	dgg := NewDiscreteGaussian(sigma / c[k-1]) // corrected to initialize dgg
	z[k-1] = dgg.Draw(
		-a[k-1] / c[k-1])

	// 2) propagate carry:  a ← a + z_{k-1}·c   (note the **plus**)
	zLastFloat := float64(z[k-1])
	for i := 0; i < k; i++ {
		a[i] += zLastFloat * c[i]
	}

	// 3) draw remaining coordinates (independent, mean = -a_i)
	for i := 0; i < k-1; i++ {
		dgg = NewDiscreteGaussian(sigma) // reinitialize dgg for each i
		z[i] = dgg.Draw(-a[i])
	}
	return z
}

// SampleGDiscrete implements the *discrete* G–sampling exactly like the
// PALISADE C++ reference GaussSampGqArbBase (Fig. 2, ePrint 2017/308).
//
//	ringQ   : R_q ring
//	sigma   : continuous Gaussian width σ_t
//	base    : gadget base t  (may be >2)
//	vCoeff  : centered syndrome coefficients, len=N=ringQ.N
//	k       : gadget length (digits)
//
// RETURNS  a k×N matrix Z of int64.
func SampleGDiscrete(
	ringQ *ring.Ring,
	sigma float64,
	base uint64,
	vCoeff []int64,
	k int,
) [][]int64 {

	N := ringQ.N
	q := ringQ.Modulus[0]                             // single-precision ring
	modDigits := baseDigits(int64(q), int64(base), k) // q_i

	//--------------------------------------------------------------------
	// 1)  L = diag(l)  +  diag(h,1)  (only σ-scaling constants)
	//--------------------------------------------------------------------
	l := make([]float64, k)
	h := make([]float64, k)

	l[0] = math.Sqrt(float64(base)*(1+1/float64(k)) + 1)
	for i := 1; i < k; i++ {
		l[i] = math.Sqrt(float64(base) * (1 + 1/float64(k-i)))
	}
	h[0] = 0
	for i := 1; i < k; i++ {
		h[i] = math.Sqrt(float64(base) * (1 - 1/float64(k-(i-1))))
	}

	//--------------------------------------------------------------------
	// 2)  cConst = constant divisor vector  (depends ONLY on modulus)
	//--------------------------------------------------------------------
	cConst := make([]float64, k) // acts as "c" in C++
	cConst[0] = float64(modDigits[0]) / float64(base)
	for i := 1; i < k; i++ {
		cConst[i] = (cConst[i-1] + float64(modDigits[i])) / float64(base)
	}

	//--------------------------------------------------------------------
	// 3)  σ'  =  σ_t /(t+1)
	//--------------------------------------------------------------------
	sigmaP := sigma / float64(base+1)

	//--------------------------------------------------------------------
	// 4)  allocate result  Z[k][N]
	//--------------------------------------------------------------------
	Z := make([][]int64, k)
	for i := range Z {
		Z[i] = make([]int64, N)
	}

	//--------------------------------------------------------------------
	// 5)  main loop over polynomial coefficients
	//--------------------------------------------------------------------
	for j := 0; j < N; j++ {
		// 5a)  base-t digits of current syndrome coeff
		vDigits := baseDigits(vCoeff[j], int64(base), k)

		// 5b)  perturbation  p   (discrete variant)
		p := Perturb(sigmaP, l, h, base) // []int64 length k

		// 5c)  build *accumulator*  a  (running carry)      [! changed]
		a := make([]float64, k)
		a[0] = float64(int64(vDigits[0])-p[0]) / float64(base)
		for i := 1; i < k; i++ {
			a[i] = (a[i-1] + float64(int64(vDigits[i])-p[i])) / float64(base)
		}

		// 5d)  sample  z  from the sparse lattice           [! changed]
		z := SampleC(cConst, sigmaP, a)
		// 5e)  recombine t-vector and store into Z          [! changed]
		//      t₀  =  base·z₀  +  q₀·z_{k−1}  +  v₀
		Z[0][j] = int64(base)*int64(z[0]) +
			int64(modDigits[0])*int64(z[k-1]) +
			int64(vDigits[0])

		//      t_i =  base·z_i − z_{i−1} + q_i·z_{k−1} + v_i   for 1≤i≤k−2
		for i := 1; i < k-1; i++ {
			Z[i][j] = int64(base)*int64(z[i]) -
				int64(z[i-1]) +
				int64(modDigits[i])*int64(z[k-1]) +
				int64(vDigits[i])
		}

		//      t_{k−1} = q_{k−1}·z_{k−1} − z_{k−2} + v_{k−1}
		Z[k-1][j] = int64(modDigits[k-1])*int64(z[k-1]) -
			int64(z[k-2]) +
			int64(vDigits[k-1])

		//! ---------- DEBUG A : verify Σ t_i·t^i = sub[j] ----------
		// if j < 8 { // limit noise; increase if needed
		// 	q := int64(ringQ.Modulus[0]) // import ringQ via closure or param
		// 	recomb := int64(0)
		// 	for idx := k - 1; idx >= 0; idx-- {
		// 		recomb = (recomb*int64(base) + Z[idx][j]) % q
		// 	}
		// 	diff := (recomb - int64(uCoeff[j])) % q
		// 	if diff < 0 {
		// 		diff += q
		// 	}
		// 	fmt.Printf("G-recomb j=%d : v0=%d  zLast=%d  t0=%d  recomb=%d  want=%d  diff=%d\n",
		// 		j, uCoeff[j], Z[k-1][j], Z[0][j], recomb, uCoeff[j], diff)
		// }
		//! ----------------------------------------------------------

	}

	return Z
}
