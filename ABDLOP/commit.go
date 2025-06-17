package ABDLOP

import (
	"math"

	ps "vSIS-Signature/Preimage_Sampler"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type Params struct {
	Q     uint64
	D     int
	N     int
	M1    int
	M2    int
	Beta1 float64
	Beta2 float64
}

type PublicKey struct {
	Ring   *ring.Ring
	Params Params
	A1     [][]*ring.Poly
	A2     [][]*ring.Poly
	B      [][]*ring.Poly
	b      []*ring.Poly
}

type Commitment struct {
	TA []*ring.Poly
	TB []*ring.Poly
	T  *ring.Poly
}

type Opening struct {
	S1 []*ring.Poly
	S2 []*ring.Poly
	M  []*ring.Poly
}

func KeyGen(p Params) *PublicKey {
	ringQ, err := ring.NewRing(p.D, []uint64{p.Q})
	if err != nil {
		panic(err)
	}
	prng, _ := utils.NewPRNG()
	uni := ring.NewUniformSampler(prng, ringQ)

	A1 := make([][]*ring.Poly, p.N)
	for i := 0; i < p.N; i++ {
		A1[i] = make([]*ring.Poly, p.M1)
		for j := 0; j < p.M1; j++ {
			poly := ringQ.NewPoly()
			uni.Read(poly)
			ringQ.NTT(poly, poly)
			A1[i][j] = poly
		}
	}

	A2 := make([][]*ring.Poly, p.N)
	for i := 0; i < p.N; i++ {
		A2[i] = make([]*ring.Poly, p.M2)
		for j := 0; j < p.M2; j++ {
			poly := ringQ.NewPoly()
			uni.Read(poly)
			ringQ.NTT(poly, poly)
			A2[i][j] = poly
		}
	}

	B := make([][]*ring.Poly, p.N)
	for i := 0; i < p.N; i++ {
		B[i] = make([]*ring.Poly, p.M2)
		for j := 0; j < p.M2; j++ {
			poly := ringQ.NewPoly()
			uni.Read(poly)
			ringQ.NTT(poly, poly)
			B[i][j] = poly
		}
	}

	bvec := make([]*ring.Poly, p.M2)
	for i := 0; i < p.M2; i++ {
		poly := ringQ.NewPoly()
		uni.Read(poly)
		ringQ.NTT(poly, poly)
		bvec[i] = poly
	}

	return &PublicKey{Ring: ringQ, Params: p, A1: A1, A2: A2, B: B, b: bvec}
}

func sampleGaussianVector(ringQ *ring.Ring, dim int, stddev float64) []*ring.Poly {
	out := make([]*ring.Poly, dim)
	dg := ps.NewDiscreteGaussian(stddev)
	for j := 0; j < dim; j++ {
		poly := ringQ.NewPoly()
		for lvl, qi := range ringQ.Modulus {
			mod := int64(qi)
			for i := 0; i < ringQ.N; i++ {
				x := dg.Draw(0)
				poly.Coeffs[lvl][i] = uint64((x%mod + mod) % mod)
			}
		}
		ringQ.NTT(poly, poly)
		out[j] = poly
	}
	return out
}

func Commit(pk *PublicKey, s1 []*ring.Poly, m []*ring.Poly) (*Commitment, *Opening) {
	ringQ := pk.Ring
	s2 := sampleGaussianVector(ringQ, pk.Params.M2, pk.Params.Beta2)

	tA := make([]*ring.Poly, pk.Params.N)
	tmp := ringQ.NewPoly()
	for i := 0; i < pk.Params.N; i++ {
		acc := ringQ.NewPoly()
		for j := 0; j < pk.Params.M1; j++ {
			ringQ.MulCoeffs(pk.A1[i][j], s1[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		for j := 0; j < pk.Params.M2; j++ {
			ringQ.MulCoeffs(pk.A2[i][j], s2[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		tA[i] = acc
	}

	tB := make([]*ring.Poly, pk.Params.N)
	for i := 0; i < pk.Params.N; i++ {
		acc := ringQ.NewPoly()
		for j := 0; j < pk.Params.M2; j++ {
			ringQ.MulCoeffs(pk.B[i][j], s2[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		mNTT := ringQ.NewPoly()
		ringQ.NTT(m[i], mNTT)
		ringQ.Add(acc, mNTT, acc)
		tB[i] = acc
	}

	tag := ringQ.NewPoly()
	for j := 0; j < pk.Params.M2; j++ {
		ringQ.MulCoeffs(pk.b[j], s2[j], tmp)
		ringQ.Add(tag, tmp, tag)
	}

	com := &Commitment{TA: tA, TB: tB, T: tag}
	open := &Opening{S1: s1, S2: s2, M: m}
	return com, open
}

// CommitWithRand is identical to Commit but uses the caller-provided s2
// randomness. It additionally returns the commitment tag b^T·s2 in NTT form.
// Panics if len(s2) != pk.Params.M2.
func CommitWithRand(pk *PublicKey, s1, s2 []*ring.Poly, m []*ring.Poly) (*Commitment, *Opening, *ring.Poly) {
	if len(s2) != pk.Params.M2 {
		panic("CommitWithRand: wrong s2 length")
	}

	ringQ := pk.Ring

	tA := make([]*ring.Poly, pk.Params.N)
	tmp := ringQ.NewPoly()
	for i := 0; i < pk.Params.N; i++ {
		acc := ringQ.NewPoly()
		for j := 0; j < pk.Params.M1; j++ {
			ringQ.MulCoeffs(pk.A1[i][j], s1[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		for j := 0; j < pk.Params.M2; j++ {
			ringQ.MulCoeffs(pk.A2[i][j], s2[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		tA[i] = acc
	}

	tB := make([]*ring.Poly, pk.Params.N)
	for i := 0; i < pk.Params.N; i++ {
		acc := ringQ.NewPoly()
		for j := 0; j < pk.Params.M2; j++ {
			ringQ.MulCoeffs(pk.B[i][j], s2[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		mNTT := ringQ.NewPoly()
		ringQ.NTT(m[i], mNTT)
		ringQ.Add(acc, mNTT, acc)
		tB[i] = acc
	}

	tag := ringQ.NewPoly()
	for j := 0; j < pk.Params.M2; j++ {
		ringQ.MulCoeffs(pk.b[j], s2[j], tmp)
		ringQ.Add(tag, tmp, tag)
	}

	com := &Commitment{TA: tA, TB: tB, T: tag}
	open := &Opening{S1: s1, S2: s2, M: m}
	return com, open, tag
}

func polyNormInf(r *ring.Ring, p *ring.Poly, q uint64) int64 {
	coeff := r.NewPoly()
	r.InvNTT(p, coeff)
	var max int64
	for _, c := range coeff.Coeffs[0] {
		v := ps.UnsignedToSigned(c, q)
		if v < 0 {
			v = -v
		}
		if v > max {
			max = v
		}
	}
	return max
}

func Open(pk *PublicKey, com *Commitment, op *Opening) bool {
	ringQ := pk.Ring

	bound1 := pk.Params.Beta1 * math.Sqrt(float64(pk.Params.D))
	bound2 := pk.Params.Beta2 * math.Sqrt(float64(pk.Params.D))

	for _, p := range op.S1 {
		if float64(polyNormInf(ringQ, p, pk.Params.Q)) > bound1 {
			return false
		}
	}
	for _, p := range op.S2 {
		if float64(polyNormInf(ringQ, p, pk.Params.Q)) > bound2 {
			return false
		}
	}

	tmp := ringQ.NewPoly()
	// recompute tA
	for i := 0; i < pk.Params.N; i++ {
		acc := ringQ.NewPoly()
		for j := 0; j < pk.Params.M1; j++ {
			ringQ.MulCoeffs(pk.A1[i][j], op.S1[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		for j := 0; j < pk.Params.M2; j++ {
			ringQ.MulCoeffs(pk.A2[i][j], op.S2[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		if !ringQ.Equal(acc, com.TA[i]) {
			return false
		}
	}

	for i := 0; i < pk.Params.N; i++ {
		acc := ringQ.NewPoly()
		for j := 0; j < pk.Params.M2; j++ {
			ringQ.MulCoeffs(pk.B[i][j], op.S2[j], tmp)
			ringQ.Add(acc, tmp, acc)
		}
		mNTT := ringQ.NewPoly()
		ringQ.NTT(op.M[i], mNTT)
		ringQ.Add(acc, mNTT, acc)
		if !ringQ.Equal(acc, com.TB[i]) {
			return false
		}
	}

	tag := ringQ.NewPoly()
	for j := 0; j < pk.Params.M2; j++ {
		ringQ.MulCoeffs(pk.b[j], op.S2[j], tmp)
		ringQ.Add(tag, tmp, tag)
	}
	if !ringQ.Equal(tag, com.T) {
		return false
	}

	return true
}
