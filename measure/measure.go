package measure

import (
	"fmt"
	"math/big"
	"os"
	"sync"
)

var Enabled bool
var Global Counter

func init() {
	Enabled = os.Getenv("MEASURE_SIZES") == "1"
	Global = Counter{M: make(map[string]int64)}
}

func BytesField(q *big.Int) int {
	// ceil(bitlen(q)/8)
	if q == nil {
		return 0
	}
	b := q.BitLen()
	return (b + 7) / 8
}

func BytesRing(phi int, q *big.Int) int {
	return phi * BytesField(q)
}

func BytesPolyDegree(d int, bytesF int) int {
	return (d + 1) * bytesF
}

func Human(n int64) string {
	const (
		KiB = 1024
		MiB = 1024 * KiB
	)
	switch {
	case n >= MiB:
		return fmt.Sprintf("%.1f MiB", float64(n)/float64(MiB))
	case n >= KiB:
		return fmt.Sprintf("%.1f KiB", float64(n)/float64(KiB))
	default:
		return fmt.Sprintf("%d B", n)
	}
}

type Counter struct {
	mu sync.Mutex
	M  map[string]int64
}

func (c *Counter) Add(key string, n int64) {
	if !Enabled {
		return
	}
	c.mu.Lock()
	c.M[key] += n
	c.mu.Unlock()
}

func (c *Counter) Dump() {
	if !Enabled {
		return
	}
	fmt.Println("[measure] Size report:")
	for k, v := range c.M {
		fmt.Printf("[measure] %s = %s\n", k, Human(v))
	}
}

func Section(name string, f func()) {
	if !Enabled {
		f()
		return
	}
	fmt.Printf("[measure] Begin %s\n", name)
	f()
	fmt.Printf("[measure] End %s\n", name)
}
