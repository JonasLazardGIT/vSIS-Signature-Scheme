package vsishash

import (
	"fmt"
	ps "vSIS-Signature/Preimage_Sampler"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// GenerateB samples a common random string B ∈ Rq^(1×4), returns
// 4 field elements in evaluation domain (CyclotomicFieldElem.Domain=Eval).
func GenerateB(ringQ *ring.Ring, prec uint, prng utils.PRNG) ([]*ps.CyclotomicFieldElem, error) {
	uni := ring.NewUniformSampler(prng, ringQ)
	B := make([]*ps.CyclotomicFieldElem, 4)
	for i := 0; i < 4; i++ {
		// 1) sample p ∈ Rq (coefficient)
		p := ringQ.NewPoly()
		uni.Read(p)
		fmt.Printf("[GenerateB] sampled p[%d][0] = %d\n", i, p.Coeffs[0][0])

		// 2) embed → CyclotomicFieldElem in Eval domain
		B[i] = ps.ConvertFromPolyBig(ringQ, p, prec)
		// Print first evaluation coefficient to trace collapse
		fmt.Printf("[GenerateB] B[%d].Coeffs[0] = %v + %vi\n", i,
			B[i].Coeffs[0].Real, B[i].Coeffs[0].Imag)
	}
	return B, nil
}

// ComputeBBSHash computes the BBS hash in BigComplex land and prints intermediate values.
func ComputeBBSHash(
	ringQ *ring.Ring,
	B []*ps.CyclotomicFieldElem,
	m, x0, x1 *ring.Poly,
	prec uint,
) (*ps.CyclotomicFieldElem, error) {
	// 1) lift m, x0, x1 into Eval domain
	mEval := ps.ConvertFromPolyBig(ringQ, m, prec)
	x0Eval := ps.ConvertFromPolyBig(ringQ, x0, prec)
	x1Eval := ps.ConvertFromPolyBig(ringQ, x1, prec)

	// 2) build "one" as a field element [1,1,1,...]
	oneC := ps.NewBigComplex(1, 0, prec)
	oneF := ps.NewFieldElemBig(ringQ.N, prec)
	oneF.Domain = ps.Eval
	for j := range oneF.Coeffs {
		oneF.Coeffs[j] = oneC.Copy()
	}

	// 3) r = B0*1 + B1*m + B2*x0
	r0 := ps.FieldMulBig(B[0], oneF)
	r1 := ps.FieldMulBig(B[1], mEval)
	r2 := ps.FieldMulBig(B[2], x0Eval)
	tmp := ps.FieldAddBig(r0, r1)
	r := ps.FieldAddBig(tmp, r2)

	// 4) denom = B3 – x1
	denom := ps.FieldSubBig(B[3], x1Eval)

	// 5) invert diag: invD * norms = diag(d)^{-1}
	invD, norms := ps.FieldInverseDiagWithNorm(denom)

	// 6) denomInv = invD ⊘ norms  (slot-wise div by real norm)
	denomInv := ps.FieldScalarDiv(invD, norms)

	// 7) hadamard division: tEval = r ⊙ denomInv
	tEval := ps.FieldMulBig(r, denomInv)

	return tEval, nil
}

// ToPolyNTT converts an Eval-domain CyclotomicFieldElem back into an
// .*ring.Poly in NTT form (ready for GaussSamp).
func ToPolyNTT(
	elem *ps.CyclotomicFieldElem,
	ringQ *ring.Ring,
) (*ring.Poly, error) {
	// 1) bring back to coefficient domain
	elem.SwitchToCoeff(ringQ)
	fmt.Printf("[ToPolyNTT] after SwitchToCoeff: elem.Coeffs[0] = %v + %vi\n", elem.Coeffs[0].Real, elem.Coeffs[0].Imag)

	// 2) interpolate via IFFT → ring.Poly (coeff), then NTT
	P := ps.ConvertToPolyBig(elem, ringQ)
	ringQ.NTT(P, P)
	fmt.Printf("[ToPolyNTT] output P.Coeffs[0][0] = %d\n", P.Coeffs[0][0])

	return P, nil
}
