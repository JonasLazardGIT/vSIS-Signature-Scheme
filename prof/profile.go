package prof

import (
	"log"
	"time"
)

// Track logs the duration since start with the given name.
func Track(start time.Time, name string) {
	elapsed := time.Since(start)
	log.Printf("%s took %s", name, elapsed)
}
