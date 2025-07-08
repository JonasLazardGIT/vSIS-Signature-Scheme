// PACS_Test.go
package PIOP

import "testing"

// TestPACSEquation7 loads the JSON fixtures, reconstructs the witness,
// builds Γ′, γ′, the masking rows Mᵢ, all F_j / F′_j polynomials, then
// computes Q = (Q₁,…,Q_ρ) as in Eq.(4) and finally checks
//
//	∀i  Σ_{ω∈Ω}  Q_i(ω)  = 0      (Equation 7)
//
// Everything else (Merkle-oracle consistency, degree tests, etc.) is
// already enforced inside the helpers we call, so this test focuses
// exclusively on the last step the verifier performs.
func TestPACSEquation7(t *testing.T) {

	// Build the Q-polynomials with fresh Fiat–Shamir randomness
	Q, omega, ringQ := BuildQFromDisk()

	// Verifier’s final check
	if !VerifyQ(ringQ, Q, omega) {
		t.Fatalf("Equation 7 failed: Σ_{ω∈Ω} Q_i(ω) ≠ 0 for some i")
	}
}
