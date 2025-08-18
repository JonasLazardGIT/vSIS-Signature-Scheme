// Parameters/system.go
package Parameters

import (
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"

	"vSIS-Signature/Preimage_Sampler"
	prof "vSIS-Signature/prof"
)

// SystemParams holds all public parameters for vSIS
type SystemParams struct {
	N      int     `json:"n"`       // ring degree
	Q      uint64  `json:"q"`       // modulus
	Base   uint64  `json:"base"`    // gadget base t
	K      int     `json:"k"`       // gadget decomposition length
	SigmaT float64 `json:"sigma_t"` // gadget‐scaled Gaussian width
	Sigma  float64 `json:"sigma"`   // raw Gaussian width
	Bound  float64 `json:"bound"`   // spectral (norm) bound s
	Beta   uint64  `json:"beta"`    // public L2 bound β
}

// BMatrix holds four Rq‐polynomials in coefficient form
type BMatrix struct {
	// B[i][j] is the j-th coefficient of the i-th polynomial
	B [][]uint64 `json:"B"`
}

// Generate computes the public parameters, writes them to
// ./Parameters/Parameters.json, and also generates the common‐random
// B‐matrix and saves it to ./Parameters/Bmatrix.json.
func Generate() {
	defer prof.Track(time.Now(), "GenerateParameters")
	// 1) Parameter arithmetic
	n := 512
	q := uint64(8399873)
	base := uint64(2)
	// k = ⌈log_base(q)⌉
	k := int(math.Ceil(math.Log(float64(q)) / math.Log(float64(base))))

	// Compute σₜ and bound s
	sigmaT, bound := Preimage_Sampler.CalculateParams(base, n, k)
	sigma := sigmaT / float64(base+1)
	beta := uint64(math.Ceil(bound))

	// 2) Instantiate R_q
	ringQ, err := ring.NewRing(n, []uint64{q})
	if err != nil {
		log.Fatalf("ring.NewRing failed: %v", err)
	}

	// 3) Bundle into our struct
	params := SystemParams{
		N:      n,
		Q:      q,
		Base:   base,
		K:      k,
		SigmaT: sigmaT,
		Sigma:  sigma,
		Bound:  bound,
		Beta:   beta,
	}

	// 4) Serialize parameters to JSON
	if err := os.MkdirAll("Parameters", 0755); err != nil {
		log.Fatalf("failed to create Parameters folder: %v", err)
	}
	paramData, err := json.MarshalIndent(params, "", "  ")
	if err != nil {
		log.Fatalf("failed to marshal parameters: %v", err)
	}
	if err := os.WriteFile("Parameters/Parameters.json", paramData, 0644); err != nil {
		log.Fatalf("failed to write Parameters.json: %v", err)
	}
	log.Println("✔ Parameters written to ./Parameters/Parameters.json")

	// 5) Generate B‐matrix (4 uniform polynomials)
	prng, err := utils.NewPRNG()
	if err != nil {
		log.Fatalf("failed to initialize PRNG: %v", err)
	}
	sampler := ring.NewUniformSampler(prng, ringQ)

	Bcoeffs := make([][]uint64, 4)
	for i := 0; i < 4; i++ {
		p := ringQ.NewPoly()
		sampler.Read(p)
		// single‐modulus ring, so level 0 contains all coefficients
		coeffCopy := make([]uint64, ringQ.N)
		copy(coeffCopy, p.Coeffs[0])
		Bcoeffs[i] = coeffCopy
	}

	bmat := BMatrix{B: Bcoeffs}
	f, err := os.Create("Parameters/Bmatrix.json")
	if err != nil {
		log.Fatalf("failed to open Bmatrix.json: %v", err)
	}
	defer f.Close()
	if _, err := f.WriteString("{\n  \"B\": [\n"); err != nil {
		log.Fatalf("write failed: %v", err)
	}
	for i, poly := range bmat.B {
		line, _ := json.Marshal(poly)
		if i < len(bmat.B)-1 {
			fmt.Fprintf(f, "    %s,\n", line)
		} else {
			fmt.Fprintf(f, "    %s\n", line)
		}
	}
	if _, err := f.WriteString("  ]\n}\n"); err != nil {
		log.Fatalf("write failed: %v", err)
	}
	log.Println("✔ B‐matrix written to ./Parameters/Bmatrix.json")
}
