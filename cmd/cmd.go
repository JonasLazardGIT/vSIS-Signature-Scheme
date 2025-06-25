// cmd/app/main.go
package main

import (
	"fmt"
	"log"

	// "vSIS-Signature/Preimage_Sampler"
	proverzkp "vSIS-Signature/ProverZKP"
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

	// 4) Build the quadratic gate and check ZKP prerequisites
	fmt.Println("⚙️  Building quadratic gate …")
	pub, priv, err := proverzkp.BuildGateFromDisk()
	if err != nil {
		log.Fatalf("gate build: %v", err)
	}
	fmt.Printf("   Gate built: |R1|=%d  witnessDim=%d\n",
		len(pub.R1), len(priv.X))

	fmt.Println("✅ All done.")
}
