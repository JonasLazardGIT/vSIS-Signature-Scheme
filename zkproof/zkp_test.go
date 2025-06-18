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
	p := abd.Params{Q: 257, D: 16, N: 2, M1: 2, M2: 2, Beta1: 3, Beta2: 3}
	pk := abd.KeyGen(p)

	witness := &Witness{
		S:  sampleGaussianVector(pk.Ring, p.M1, p.Beta1),
		U:  []*ring.Poly{},
		X0: sampleGaussianVector(pk.Ring, 1, p.Beta1),
		X1: sampleGaussianVector(pk.Ring, 1, p.Beta1),
		S2: sampleGaussianVector(pk.Ring, p.M2, p.Beta2),
	}

	// Build random gate parameters
	// Use a zero gate to focus on the linear Σ-protocol checks
	L := p.M1 + len(witness.U) + len(witness.X0)
	R2 := make([][]*ring.Poly, L+1)
	for i := range R2 {
		R2[i] = make([]*ring.Poly, L+1)
		for j := range R2[i] {
			R2[i][j] = pk.Ring.NewPoly()
		}
	}
	R1 := make([]*ring.Poly, L+1)
	for i := range R1 {
		R1[i] = pk.Ring.NewPoly()
	}
	gate := &QuadraticGate{R2: R2, R1: R1, R0: pk.Ring.NewPoly()}

	var transcript *Transcript
	for {
		transcript = Prove(pk, gate, witness)
		if transcript.C == 1 {
			break
		}
	}
	t.Logf("challenge=%d", transcript.C)

	if !Verify(pk, gate, transcript.Com, transcript) {
		t.Fatalf("Sigma protocol verification failed")
	}
}
