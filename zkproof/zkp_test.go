package zkproof

import (
	"testing"

	abd "vSIS-Signature/ABDLOP"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func makeUniformMatrix(r *ring.Ring, rows, cols int) [][]*ring.Poly {
	prng, _ := utils.NewPRNG()
	uni := ring.NewUniformSampler(prng, r)
	mat := make([][]*ring.Poly, rows)
	for i := 0; i < rows; i++ {
		mat[i] = make([]*ring.Poly, cols)
		for j := 0; j < cols; j++ {
			mat[i][j] = r.NewPoly()
			uni.Read(mat[i][j])
			r.NTT(mat[i][j], mat[i][j])
		}
	}
	return mat
}

func TestSigmaProtocol(t *testing.T) {
	p := abd.Params{Q: 257, D: 16, N: 2, M1: 2, M2: 2, Beta1: 2, Beta2: 2}
	pk := abd.KeyGen(p)

	// Build random gate parameters
	b1 := makeUniformMatrix(pk.Ring, p.N, 1)
	// flatten b1 to slice length N
	b1vec := make([]*ring.Poly, p.N)
	for i := range b1vec {
		b1vec[i] = b1[i][0]
	}
	A := pk.A1
	B0 := makeUniformMatrix(pk.Ring, p.N, 1)

	gate := BuildQuadraticGate(pk.Ring, b1vec, A, B0, p.M1, 0, 0)

	witness := &Witness{
		S:  sampleGaussianVector(pk.Ring, p.M1, p.Beta1),
		U:  []*ring.Poly{},
		X0: sampleGaussianVector(pk.Ring, 1, p.Beta1),
		X1: sampleGaussianVector(pk.Ring, p.N, p.Beta1),
		S2: sampleGaussianVector(pk.Ring, p.M2, p.Beta2),
	}

	transcript := Prove(pk, gate, witness)
	s1 := append(append(witness.S, witness.U...), witness.X0...)
	com, _, _ := abd.CommitWithRand(pk, s1, witness.S2, witness.X1)

	if !Verify(pk, gate, com, transcript) {
		t.Fatalf("Sigma protocol verification failed")
	}
}
