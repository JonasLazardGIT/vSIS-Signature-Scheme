package Preimage_Sampler

import "testing"

// TestZtoZhatRoundTrip ensures that NTT + InvNTT round-trips each gadget row.
func TestZtoZhatRoundTrip(t *testing.T) {
	const n = 16
	const q = uint64(97)
	ringQ := makeSmallRing(n, q)

	for j := 0; j < n; j++ {
		row := make([]int64, n)
		row[j] = 1
		Z := [][]int64{row}
		polys := ZtoZhat(Z, ringQ)
		if len(polys) != 1 {
			t.Fatalf("expected 1 poly, got %d", len(polys))
		}
		coeffs := ringQ.NewPoly()
		ringQ.InvNTT(polys[0], coeffs)
		for i := 0; i < n; i++ {
			var want uint64
			if i == j {
				want = 1
			} else {
				want = 0
			}
			if coeffs.Coeffs[0][i] != want {
				t.Fatalf("delta at %d mismatch index %d: got %d want %d", j, i, coeffs.Coeffs[0][i], want)
			}
		}
	}
}
