// cmd/app/main.go
package main

import (
	"fmt"
	"log"

	// "vSIS-Signature/Preimage_Sampler"
	"vSIS-Signature/PIOP"
	signer "vSIS-Signature/Signer"
	Parameters "vSIS-Signature/System"
	verifier "vSIS-Signature/Verifier"
)

func main() {
	// Preimage_Sampler.Main()
	// 1) Generate (or load) public system parameters
	fmt.Println("🔧 Generating public parameters...")
	Parameters.Generate()

	// 2) Run the signer: key-gen + sign one syndrome
	fmt.Println("✍️  Generating keypair and signature...")
	signer.Sign()

	// 3) Verify the signature as before
	fmt.Println("🔍 Verifying signature...")
	if ok := verifier.Verify(); !ok {
		log.Fatal("❌ Signature verification failed")
	}

	fmt.Println("⚙️  Building Witnesses …")
	w1, w2, w3 := PIOP.BuildWitnessFromDisk()
	fmt.Printf("   Gate built, lengths are : w1: %d, w2: %d, w3: %d\n",
		len(w1), len(w2.Coeffs[0]), len(w3))

	fmt.Println("✅ All done.")
}
