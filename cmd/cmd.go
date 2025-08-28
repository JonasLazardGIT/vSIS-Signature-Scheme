// cmd/app/main.go
package main

import (
    "fmt"
    "log"
    "time"

    // "vSIS-Signature/Preimage_Sampler"

    signer "vSIS-Signature/Signer"
    Parameters "vSIS-Signature/System"
    verifier "vSIS-Signature/Verifier"
    measure "vSIS-Signature/measure"
    prof "vSIS-Signature/prof"
)

func main() {
	defer prof.Track(time.Now(), "main")
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

	// Emit a single consolidated measure report for the whole run
	measure.Global.Dump()

	fmt.Println("✅ All done.")
}
