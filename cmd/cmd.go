// cmd/app/main.go
package main

import (
	"fmt"
	"log"

	"vSIS-Signature/Preimage_Sampler"
	signer "vSIS-Signature/Signer"
	Parameters "vSIS-Signature/System"
	verifier "vSIS-Signature/Verifier"
)

func main() {
	Preimage_Sampler.Main()
	// 1) Generate (or load) public system parameters
	fmt.Println("🔧 Generating public parameters...")
	Parameters.Generate()

	// 2) Run the signer: key-gen + sign one syndrome
	fmt.Println("✍️  Generating keypair and signature...")
	signer.Sign()

	// 3) Verify the freshly produced signature
	fmt.Println("🔍 Verifying signature...")
	if ok := verifier.Verify(); !ok {
		log.Fatal("❌ Signature verification failed")
	}

	fmt.Println("✅ All done: signature is valid.")
}
