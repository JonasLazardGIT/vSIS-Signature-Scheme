package zkproof

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ring"

	abd "vSIS-Signature/ABDLOP"
)

// Verify executes the verifier algorithm. This is a placeholder that
// always returns false until fully implemented.
func Verify(pk *abd.PublicKey, gate *QuadraticGate, com *abd.Commitment, transcript *Transcript) bool {
	fmt.Println("Verifier: start")
	ringQ := pk.Ring
	c := transcript.C

	bound1 := pk.Params.Beta1 * math.Sqrt(float64(pk.Params.D*2))
	bound2 := pk.Params.Beta2 * math.Sqrt(float64(pk.Params.D*2))

	for _, p := range transcript.Z1 {
		if float64(abd.PolyNormInf(ringQ, p, pk.Params.Q)) > bound1 {
			return false
		}
	}
	for _, p := range transcript.Z2 {
		if float64(abd.PolyNormInf(ringQ, p, pk.Params.Q)) > bound2 {
			return false
		}
	}

	leftA := matVecMul(pk.A1, transcript.Z1, ringQ)
	rightA := matVecMul(pk.A2, transcript.Z2, ringQ)
	left := make([]*ring.Poly, len(leftA))
	for i := range left {
		left[i] = ringQ.NewPoly()
		ringQ.Add(leftA[i], rightA[i], left[i])
	}

	right := make([]*ring.Poly, len(com.TA))
	for i := range right {
		right[i] = ringQ.NewPoly()
		switch c {
		case 1:
			ringQ.Add(com.TA[i], transcript.W[i], right[i])
		case -1:
			tmp := ringQ.NewPoly()
			ringQ.Neg(com.TA[i], tmp)
			ringQ.Add(transcript.W[i], tmp, right[i])
		default:
			ring.Copy(transcript.W[i], right[i])
		}
	}
	for i := range left {
		if !ringQ.Equal(left[i], right[i]) {
			return false
		}
	}

	ctB := scalarMul(com.TB, c, ringQ)
	diffB := make([]*ring.Poly, len(ctB))
	mulB := matVecMul(pk.B, transcript.Z2, ringQ)
	for i := range diffB {
		diffB[i] = ringQ.NewPoly()
		ringQ.Sub(ctB[i], mulB[i], diffB[i])
	}
	zLeft := applyAutoNTT(transcript.Z1, 0, ringQ)
	zRight := applyAutoNTT(diffB, 0, ringQ)
	z := append(zLeft, zRight...)

	f := innerProd(z, matVecMul(gate.R2, z, ringQ), ringQ)
	r1z := innerProd(gate.R1, z, ringQ)
	if c != 0 {
		tmp := ringQ.NewPoly()
		if c == 1 {
			ring.Copy(r1z, tmp)
		} else {
			ringQ.Neg(r1z, tmp)
		}
		ringQ.Add(f, tmp, f)
		if c*c != 0 {
			tmp0 := ringQ.NewPoly()
			ringQ.MulScalar(gate.R0, uint64(c*c), tmp0)
			ringQ.Add(f, tmp0, f)
		}
	}
	bZ2 := innerProd(pk.Bvec, transcript.Z2, ringQ)
	ringQ.Sub(transcript.T, bZ2, bZ2)
	ringQ.Sub(f, bZ2, f)

	if !ringQ.Equal(f, transcript.V) {
		return false
	}
	fmt.Println("Verifier: accept")
	return true
}
