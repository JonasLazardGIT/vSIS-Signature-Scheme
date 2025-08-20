//go:build analysis
// +build analysis

package main

import (
	"encoding/json"
	"fmt"
	"image/color"
	"log"
	"math"
	"os"
	"path/filepath"
	"sort"
	"time"

	signer "vSIS-Signature/Signer"
	Parameters "vSIS-Signature/System"

	"github.com/tuneinsight/lattigo/v4/ring"

	// Static plots
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"

	// Interactive HTML histogram
	"github.com/go-echarts/go-echarts/v2/charts"
	"github.com/go-echarts/go-echarts/v2/components"
	"github.com/go-echarts/go-echarts/v2/opts"
)

// ----------------------------- I/O: save centered coeffs -----------------------------

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
	if err := os.MkdirAll(dstDir, 0o755); err != nil {
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
			continue // skip non-matching files gracefully
		}
		for _, row := range sig.Coeffs {
			for _, c := range row {
				values = append(values, float64(c))
			}
		}
	}
	return values, nil
}

// ----------------------------------- Statistics -------------------------------------

type summaryStats struct {
	Count    int     `json:"count"`
	Mean     float64 `json:"mean"`
	Std      float64 `json:"std"`
	Min      float64 `json:"min"`
	Q1       float64 `json:"q1"`
	Median   float64 `json:"median"`
	Q3       float64 `json:"q3"`
	Max      float64 `json:"max"`
	IQR      float64 `json:"iqr"`
	Skewness float64 `json:"skewness"`
	Kurtosis float64 `json:"kurtosis_excess"`
}

func computeStats(x []float64) summaryStats {
	n := len(x)
	if n == 0 {
		return summaryStats{}
	}
	cp := append([]float64(nil), x...)
	sort.Float64s(cp)

	min, max := cp[0], cp[n-1]
	median := quantileSorted(cp, 0.5)
	q1 := quantileSorted(cp, 0.25)
	q3 := quantileSorted(cp, 0.75)
	iqr := q3 - q1

	var m float64
	for _, v := range x {
		m += v
	}
	m /= float64(n)

	var m2, m3, m4 float64
	for _, v := range x {
		d := v - m
		d2 := d * d
		m2 += d2
		m3 += d2 * d
		m4 += d2 * d2
	}
	varVar := m2 / float64(n-1)
	std := math.Sqrt(varVar)

	var skew, kurtEx float64
	if std > 0 {
		m2n := m2 / float64(n)
		m3n := m3 / float64(n)
		m4n := m4 / float64(n)
		skew = m3n / math.Pow(m2n, 1.5)
		kurtEx = m4n/m2n/m2n - 3.0
	}

	return summaryStats{
		Count:    n,
		Mean:     m,
		Std:      std,
		Min:      min,
		Q1:       q1,
		Median:   median,
		Q3:       q3,
		Max:      max,
		IQR:      iqr,
		Skewness: skew,
		Kurtosis: kurtEx,
	}
}

// quantileSorted returns the p-quantile of a sorted slice using linear interpolation.
func quantileSorted(sorted []float64, p float64) float64 {
	if p <= 0 {
		return sorted[0]
	}
	if p >= 1 {
		return sorted[len(sorted)-1]
	}
	pos := p * float64(len(sorted)-1)
	l := int(math.Floor(pos))
	r := int(math.Ceil(pos))
	if l == r {
		return sorted[l]
	}
	weight := pos - float64(l)
	return sorted[l]*(1-weight) + sorted[r]*weight
}

// Freedman–Diaconis recommended bin count; clamped to [50, 2000].
func freedmanDiaconisBins(x []float64) int {
	n := len(x)
	if n < 2 {
		return 1
	}
	cp := append([]float64(nil), x...)
	sort.Float64s(cp)
	iqr := quantileSorted(cp, 0.75) - quantileSorted(cp, 0.25)
	if iqr == 0 {
		return minInt(200, n)
	}
	binWidth := 2 * iqr * math.Pow(float64(n), -1.0/3.0)
	if binWidth <= 0 {
		return minInt(200, n)
	}
	r := cp[n-1] - cp[0]
	k := int(math.Ceil(r / binWidth))
	if k < 50 {
		k = 50
	}
	if k > 2000 {
		k = 2000
	}
	return k
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// computeHistogram returns equally spaced edges and counts for values.
func computeHistogram(values []float64, nbins int) (edges []float64, counts []int) {
	cp := append([]float64(nil), values...)
	sort.Float64s(cp)
	minv, maxv := cp[0], cp[len(cp)-1]
	if nbins < 1 {
		nbins = 1
	}
	width := (maxv - minv) / float64(nbins)
	if width <= 0 {
		width = 1
	}
	edges = make([]float64, nbins+1)
	for i := 0; i <= nbins; i++ {
		edges[i] = minv + float64(i)*width
	}
	counts = make([]int, nbins)
	for _, v := range values {
		idx := int(math.Floor((v - minv) / width))
		if idx < 0 {
			idx = 0
		}
		if idx >= nbins {
			idx = nbins - 1
		}
		counts[idx]++
	}
	return
}

// ---------------------------- Static PNG: fine histogram -----------------------------

// plotHistogramPNG plots a high-resolution histogram with mean/±1σ lines.
func plotHistogramPNG(values []float64, outPNG string, stats summaryStats) error {
	if len(values) == 0 {
		return fmt.Errorf("no values to plot")
	}

	p := plot.New()
	p.Title.Text = "Coefficient Distribution (Fine Binning)"
	p.X.Label.Text = "coefficient value"
	p.Y.Label.Text = "count"

	nbins := freedmanDiaconisBins(values)
	h, err := plotter.NewHist(plotter.Values(values), nbins)
	if err != nil {
		return err
	}
	h.Normalize(0) // keep counts
	p.Add(h)

	// We need an explicit y-extent for vertical lines; compute from our own binning.
	_, counts := computeHistogram(values, nbins)
	var ymax float64
	for _, c := range counts {
		if float64(c) > ymax {
			ymax = float64(c)
		}
	}
	if ymax <= 0 {
		ymax = 1
	}

	addVLine := func(x float64, width vg.Length, col color.Color, label string) error {
		pts := plotter.XYs{{X: x, Y: 0}, {X: x, Y: ymax}}
		l, err := plotter.NewLine(pts)
		if err != nil {
			return err
		}
		l.Width = width
		l.Color = col
		p.Add(l)
		p.Legend.Add(label, l)
		return nil
	}

	// Mean (thicker), std bounds (thinner)
	if err := addVLine(stats.Mean, 1.5, color.RGBA{0, 0, 0, 255}, fmt.Sprintf("mean=%.3f", stats.Mean)); err != nil {
		return err
	}
	if stats.Std > 0 {
		if err := addVLine(stats.Mean-stats.Std, 0.8, color.RGBA{80, 80, 80, 220}, fmt.Sprintf("mean-σ=%.3f", stats.Mean-stats.Std)); err != nil {
			return err
		}
		if err := addVLine(stats.Mean+stats.Std, 0.8, color.RGBA{80, 80, 80, 220}, fmt.Sprintf("mean+σ=%.3f", stats.Mean+stats.Std)); err != nil {
			return err
		}
	}

	// Save very high-res for crisp zoom
	if err := p.Save(12*vg.Inch, 8*vg.Inch, outPNG); err != nil {
		return err
	}
	return nil
}

// ------------------------------ Interactive HTML plot --------------------------------

// buildInteractiveHistogram writes a zoomable HTML histogram with tooltips.
func buildInteractiveHistogram(values []float64, outHTML string, stats summaryStats) error {
	if len(values) == 0 {
		return fmt.Errorf("no values to plot")
	}

	// Bucketization (use same FD bins so HTML ≈ PNG)
	nbins := freedmanDiaconisBins(values)
	edges, counts := computeHistogram(values, nbins)

	// label each bin by its center
	xLabels := make([]string, nbins)
	for i := 0; i < nbins; i++ {
		center := 0.5 * (edges[i] + edges[i+1])
		xLabels[i] = fmt.Sprintf("%.2f", center)
	}

	bar := charts.NewBar()
	bar.SetGlobalOptions(
		charts.WithTitleOpts(opts.Title{
			Title: "Coefficient Distribution (Interactive)",
			Subtitle: fmt.Sprintf("n=%d, mean=%.3f, std=%.3f, median=%.3f, IQR=%.3f",
				stats.Count, stats.Mean, stats.Std, stats.Median, stats.IQR),
		}),
		charts.WithInitializationOpts(opts.Initialization{
			PageTitle: "Coefficient Histogram",
			Width:     "1200px",
			Height:    "700px",
		}),
		// zoom/pan
		charts.WithDataZoomOpts(
			opts.DataZoom{Type: "inside"},
			opts.DataZoom{Type: "slider"},
		),
		// default tooltip on
	)
	bar.SetXAxis(xLabels).
		AddSeries("count", toBarItems(counts)).
		SetSeriesOptions(
			charts.WithLabelOpts(opts.Label{Show: opts.Bool(false)}),
		)

	page := components.NewPage()
	page.AddCharts(bar)

	f, err := os.Create(outHTML)
	if err != nil {
		return err
	}
	defer f.Close()
	return page.Render(f)
}

func toBarItems(vals []int) []opts.BarData {
	out := make([]opts.BarData, len(vals))
	for i, v := range vals {
		out[i] = opts.BarData{Value: v}
	}
	return out
}

// ----------------------------- Stats export and helpers ------------------------------

func saveStatsJSON(stats summaryStats, path string) error {
	b, err := json.MarshalIndent(stats, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0o644)
}

// --------------------------------------- main ---------------------------------------

func main() {
	const runs = 10
	const outDir = "Read_signatures"

	Parameters.Generate()
	for i := 0; i < runs; i++ {
		log.Printf("Run %d/%d", i+1, runs)
		signer.Sign()
		if err := saveCoeffSignature(outDir); err != nil {
			log.Fatalf("saveCoeffSignature: %v", err)
		}
	}

	values, err := collectCoeffs(outDir)
	if err != nil {
		log.Fatalf("collectCoeffs: %v", err)
	}

	stats := computeStats(values)
	// Save stats JSON and print concise line
	if err := saveStatsJSON(stats, filepath.Join(outDir, "coefficient_stats.json")); err != nil {
		log.Printf("warn: save stats json: %v", err)
	}
	log.Printf("n=%d mean=%.6f std=%.6f min=%.0f median=%.0f max=%.0f IQR=%.0f skew=%.3f kurt(excess)=%.3f",
		stats.Count, stats.Mean, stats.Std, stats.Min, stats.Median, stats.Max, stats.IQR, stats.Skewness, stats.Kurtosis)

	// Static high-res PNG (fine bins + mean/σ lines)
	pngPath := filepath.Join(outDir, "coefficient_distribution.png")
	if err := plotHistogramPNG(values, pngPath, stats); err != nil {
		log.Fatalf("plotHistogramPNG: %v", err)
	}
	fmt.Println("High-res histogram saved to", pngPath)

	// Interactive HTML (zoomable)
	htmlPath := filepath.Join(outDir, "coefficient_distribution.html")
	if err := buildInteractiveHistogram(values, htmlPath, stats); err != nil {
		log.Fatalf("buildInteractiveHistogram: %v", err)
	}
	fmt.Println("Interactive histogram saved to", htmlPath)
}
