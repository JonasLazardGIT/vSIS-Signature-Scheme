package ABDLOP

import (
	"testing"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func TestCommitOpen(t *testing.T) {
	p := Params{Q: 257, D: 16, N: 2, M1: 2, M2: 2, Beta1: 2, Beta2: 2}
	pk := KeyGen(p)

	s1 := sampleGaussianVector(pk.Ring, p.M1, p.Beta1)

	prng, _ := utils.NewPRNG()
	uni := ring.NewUniformSampler(prng, pk.Ring)
	m := make([]*ring.Poly, p.N)
	for i := 0; i < p.N; i++ {
		m[i] = pk.Ring.NewPoly()
		uni.Read(m[i])
	}

	com, open := Commit(pk, s1, m)
	if com.T == nil {
		t.Fatalf("commitment tag is nil")
	}
	if !Open(pk, com, open) {
		t.Fatalf("open failed on valid commitment")
	}
}

func TestOpenBoundFail(t *testing.T) {
	p := Params{Q: 257, D: 16, N: 2, M1: 2, M2: 2, Beta1: 1, Beta2: 1}
	pk := KeyGen(p)
	s1 := sampleGaussianVector(pk.Ring, p.M1, p.Beta1)

	prng, _ := utils.NewPRNG()
	uni := ring.NewUniformSampler(prng, pk.Ring)
	m := make([]*ring.Poly, p.N)
	for i := 0; i < p.N; i++ {
		m[i] = pk.Ring.NewPoly()
		uni.Read(m[i])
	}

	com, open := Commit(pk, s1, m)

	// tamper with s2 to exceed bound
	bad := &Opening{S1: open.S1, S2: make([]*ring.Poly, len(open.S2)), M: open.M}
	for i, poly := range open.S2 {
		coeff := pk.Ring.NewPoly()
		pk.Ring.InvNTT(poly, coeff)
		if i == 0 {
			coeff.Coeffs[0][0] = 5
		}
		pk.Ring.NTT(coeff, coeff)
		bad.S2[i] = coeff
	}

	if Open(pk, com, bad) {
		t.Fatalf("open succeeded with malformed witness")
	}
}
