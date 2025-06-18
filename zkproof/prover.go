package zkproof

import (
	"fmt"
	"math"
	"math/rand"

	abd "vSIS-Signature/ABDLOP"
	ps "vSIS-Signature/Preimage_Sampler"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// Witness represents (s,u,x0,x1) with all polynomials in NTT form.
type Witness struct {
	S  []*ring.Poly
	U  []*ring.Poly
	X0 []*ring.Poly
	X1 []*ring.Poly
	S2 []*ring.Poly
}

// Transcript captures a single execution of the Σ-protocol.
type Transcript struct {
	W  []*ring.Poly
	T  *ring.Poly
	V  *ring.Poly
	Z1 []*ring.Poly
	Z2 []*ring.Poly
	C  int
}

// sampleGaussianVector is a helper that samples dim polynomials from D_β.
func sampleGaussianVector(r *ring.Ring, dim int, beta float64) []*ring.Poly {
	out := make([]*ring.Poly, dim)
	dg := ps.NewDiscreteGaussian(beta)
	for j := 0; j < dim; j++ {
		poly := r.NewPoly()
		for lvl, qi := range r.Modulus {
			mod := int64(qi)
			for i := 0; i < r.N; i++ {
				x := dg.Draw(0)
				poly.Coeffs[lvl][i] = uint64((x%mod + mod) % mod)
			}
		}
		out[j] = poly
	}
	return out
}

func sampleGaussianScalar(r *ring.Ring, beta float64) *ring.Poly {
	return sampleGaussianVector(r, 1, beta)[0]
}

func makeNTTslice(r *ring.Ring, in []*ring.Poly) []*ring.Poly {
	out := make([]*ring.Poly, len(in))
	for i, p := range in {
		out[i] = r.NewPoly()
		r.NTT(p, out[i])
	}
	return out
}

// makeInvNTTslice converts each polynomial of in with InvNTT and returns the new slice.
func makeInvNTTslice(r *ring.Ring, in []*ring.Poly) []*ring.Poly {
	out := make([]*ring.Poly, len(in))
	for i, p := range in {
		out[i] = r.NewPoly()
		r.InvNTT(p, out[i])
	}
	return out
}

func toNTT(r *ring.Ring, p *ring.Poly) *ring.Poly {
	out := r.NewPoly()
	r.NTT(p, out)
	return out
}

func matVecMul(mat [][]*ring.Poly, vec []*ring.Poly, r *ring.Ring) []*ring.Poly {
	out := make([]*ring.Poly, len(mat))
	tmp := r.NewPoly()
	for i := 0; i < len(mat); i++ {
		acc := r.NewPoly()
		for j := 0; j < len(mat[i]); j++ {
			r.MulCoeffs(mat[i][j], vec[j], tmp)
			r.Add(acc, tmp, acc)
		}
		out[i] = acc
	}
	return out
}

func innerProd(a, b []*ring.Poly, r *ring.Ring) *ring.Poly {
	res := r.NewPoly()
	tmp := r.NewPoly()
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		r.MulCoeffs(a[i], b[i], tmp)
		r.Add(res, tmp, res)
	}
	return res
}

func applyAutoNTT(vec []*ring.Poly, k int, r *ring.Ring) []*ring.Poly {
	out := make([]*ring.Poly, len(vec))
	for i, p := range vec {
		out[i] = r.NewPoly()
		r.Shift(p, k, out[i])
	}
	return out
}

func negateSlice(vec []*ring.Poly, r *ring.Ring) []*ring.Poly {
	out := make([]*ring.Poly, len(vec))
	for i, p := range vec {
		out[i] = r.NewPoly()
		r.Neg(p, out[i])
	}
	return out
}

func scalarMul(vec []*ring.Poly, c int, r *ring.Ring) []*ring.Poly {
	out := make([]*ring.Poly, len(vec))
	for i, p := range vec {
		out[i] = r.NewPoly()
		switch c {
		case 1:
			ring.Copy(p, out[i])
		case -1:
			r.Neg(p, out[i])
		default:
		}
	}
	return out
}

// Prove executes the prover algorithm. This is a simplified placeholder
// implementation matching the published interface.
func Prove(pk *abd.PublicKey, gate *QuadraticGate, witness *Witness) *Transcript {
	fmt.Println("Prover: start")
	ringQ := pk.Ring

	// 1. witness split
	s1 := append(append(witness.S, witness.U...), witness.X0...)
	mPoly := witness.X1

	var (
		y1NTT []*ring.Poly
		y2NTT []*ring.Poly
		s1NTT []*ring.Poly
		s2NTT []*ring.Poly
		com   *abd.Commitment
		open  *abd.Opening
		v     *ring.Poly
		t     *ring.Poly
		w     []*ring.Poly
		z1    []*ring.Poly
		z2    []*ring.Poly
		c     int
	)

	bound1 := pk.Params.Beta1 * math.Sqrt(float64(pk.Params.D*2))
	bound2 := pk.Params.Beta2 * math.Sqrt(float64(pk.Params.D*2))

	for {
		// 2. masks
		y1 := sampleGaussianVector(ringQ, len(s1), pk.Params.Beta1)
		y2Vec := sampleGaussianVector(ringQ, pk.Params.M2, pk.Params.Beta2)
		y1NTT = makeNTTslice(ringQ, y1)
		y2NTT = makeNTTslice(ringQ, y2Vec)

		// 3. commit using real witness s2
		com, open, _ = abd.CommitWithRand(pk, s1, witness.S2, mPoly)
		s1NTT = open.S1
		s2NTT = open.S2

		// compute w = A1*y1 + A2*y2
		part1 := matVecMul(pk.A1, y1NTT, ringQ)
		part2 := matVecMul(pk.A2, y2NTT, ringQ)
		w = make([]*ring.Poly, len(part1))
		for i := 0; i < len(part1); i++ {
			w[i] = ringQ.NewPoly()
			ringQ.Add(part1[i], part2[i], w[i])
		}

		// 4. automorphisms (only i=0 supported)
		yLeft := applyAutoNTT(y1NTT, 0, ringQ)
		by2 := matVecMul(pk.B, y2NTT, ringQ)
		yRight := negateSlice(applyAutoNTT(by2, 0, ringQ), ringQ)
		y := append(yLeft, yRight...)

		// helper vector s = s1||s2
		sVec := append(s1NTT, s2NTT...)

		// 5. scalars
		g1 := innerProd(sVec, matVecMul(gate.R2, y, ringQ), ringQ)
		ringQ.Add(g1, innerProd(y, matVecMul(gate.R2, sVec, ringQ), ringQ), g1)
		ringQ.Add(g1, innerProd(gate.R1, y, ringQ), g1)

		t := ring.NewPoly(ringQ.N, len(ringQ.Modulus)-1)
		ring.Copy(com.T, t)
		ringQ.Add(t, g1, t)

		v := innerProd(y, matVecMul(gate.R2, y, ringQ), ringQ)
		ringQ.Add(v, innerProd(pk.Bvec, y2NTT, ringQ), v)

		// 6. send first msg - compute challenge
		c = rand.Intn(3) - 1

		// 7. response
		z1 = make([]*ring.Poly, len(s1NTT))
		for i := range z1 {
			switch c {
			case 1:
				z1[i] = ringQ.NewPoly()
				ringQ.Add(s1NTT[i], y1NTT[i], z1[i])
			case -1:
				tmp := ringQ.NewPoly()
				ringQ.Neg(s1NTT[i], tmp)
				z1[i] = ringQ.NewPoly()
				ringQ.Add(tmp, y1NTT[i], z1[i])
			default:
				z1[i] = ring.NewPoly(ringQ.N, len(ringQ.Modulus)-1)
				ring.Copy(y1NTT[i], z1[i])
			}
		}

		z2 = make([]*ring.Poly, len(s2NTT))
		for i := range z2 {
			switch c {
			case 1:
				z2[i] = ringQ.NewPoly()
				ringQ.Add(s2NTT[i], y2NTT[i], z2[i])
			case -1:
				tmp := ringQ.NewPoly()
				ringQ.Neg(s2NTT[i], tmp)
				z2[i] = ringQ.NewPoly()
				ringQ.Add(tmp, y2NTT[i], z2[i])
			default:
				z2[i] = ring.NewPoly(ringQ.N, len(ringQ.Modulus)-1)
				ring.Copy(y2NTT[i], z2[i])
			}
		}

		// rejection sampling
		z1Coeff := makeInvNTTslice(ringQ, z1)
		z2Coeff := makeInvNTTslice(ringQ, z2)
		reject := false
		for _, p := range z1Coeff {
			if float64(abd.PolyNormInf(ringQ, p, pk.Params.Q)) > bound1 {
				reject = true
				break
			}
		}
		if !reject {
			for _, p := range z2Coeff {
				if float64(abd.PolyNormInf(ringQ, p, pk.Params.Q)) > bound2 {
					reject = true
					break
				}
			}
		}
		if !reject {
			break
		}
	}

	fmt.Println("Prover: end")
	return &Transcript{W: w, T: t, V: v, Z1: z1, Z2: z2, C: c}
}
