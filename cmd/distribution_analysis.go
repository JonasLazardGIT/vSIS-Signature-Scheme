//go:build analysis
// +build analysis

package main

import (
	"encoding/json"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"time"

	signer "vSIS-Signature/Signer"
	Parameters "vSIS-Signature/System"

	"github.com/tuneinsight/lattigo/v4/ring"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

// saveCoeffSignature reads the last signature and stores its centered coefficients in dstDir.
func saveCoeffSignature(dstDir string) error {
	params, err := signer.LoadParams("Parameters/Parameters.json")
	if err != nil {
		return err
	}
	ringQ, err := ring.NewRing(params.N, []uint64{params.Q})
	if err != nil {
		return err
	}
	raw, err := os.ReadFile("Signature/Signature.json")
	if err != nil {
		return err
	}
	var sig struct {
		Signature [][]uint64 `json:"signature"`
	}
	if err := json.Unmarshal(raw, &sig); err != nil {
		return err
	}
	halfQ := params.Q / 2
	coeffSigs := make([][]int64, len(sig.Signature))
	for i, row := range sig.Signature {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], row)
		ringQ.InvNTT(p, p)
		coeffRow := make([]int64, params.N)
		for j, cu := range p.Coeffs[0] {
			if cu > halfQ {
				coeffRow[j] = int64(cu) - int64(params.Q)
			} else {
				coeffRow[j] = int64(cu)
			}
		}
		coeffSigs[i] = coeffRow
	}
	out := struct {
		Timestamp string    `json:"timestamp"`
		Coeffs    [][]int64 `json:"coeffs"`
	}{
		Timestamp: time.Now().Format("20060102_150405"),
		Coeffs:    coeffSigs,
	}
	if err := os.MkdirAll(dstDir, 0755); err != nil {
		return err
	}
	fname := filepath.Join(dstDir, fmt.Sprintf("coeffsig_%s.json", out.Timestamp))
	f, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer f.Close()
	enc := json.NewEncoder(f)
	enc.SetIndent("", "  ")
	if err := enc.Encode(out); err != nil {
		return err
	}
	log.Printf("saved coefficient signature to %s", fname)
	return nil
}

// collectCoeffs reads all signature files in dir and returns flattened coefficients.
func collectCoeffs(dir string) ([]float64, error) {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return nil, err
	}
	var values []float64
	for _, e := range entries {
		if e.IsDir() {
			continue
		}
		data, err := os.ReadFile(filepath.Join(dir, e.Name()))
		if err != nil {
			return nil, err
		}
		var sig struct {
			Coeffs [][]int64 `json:"coeffs"`
		}
		if err := json.Unmarshal(data, &sig); err != nil {
			return nil, err
		}
		for _, row := range sig.Coeffs {
			for _, c := range row {
				values = append(values, float64(c))
			}
		}
	}
	return values, nil
}

// plotHistogram plots the histogram of values and saves it to path.
func plotHistogram(values []float64, path string) error {
	p := plot.New()
	p.Title.Text = "Coefficient Distribution"
	h, err := plotter.NewHist(plotter.Values(values), 50)
	if err != nil {
		return err
	}
	p.Add(h)
	if err := p.Save(6*vg.Inch, 4*vg.Inch, path); err != nil {
		return err
	}
	return nil
}

func main() {
	const runs = 10

	Parameters.Generate()
	for i := 0; i < runs; i++ {
		log.Printf("Run %d/%d", i+1, runs)
		signer.Sign()
		if err := saveCoeffSignature("Read_signatures"); err != nil {
			log.Fatalf("saveCoeffSignature: %v", err)
		}
	}
	values, err := collectCoeffs("Read_signatures")
	if err != nil {
		log.Fatalf("collectCoeffs: %v", err)
	}
	if err := plotHistogram(values, filepath.Join("Read_signatures", "coefficient_distribution.png")); err != nil {
		log.Fatalf("plotHistogram: %v", err)
	}
	fmt.Println("Histogram saved to Read_signatures/coefficient_distribution.png")
}
