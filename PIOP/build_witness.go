// quadraticgate/quadratic_gate.go
package PIOP

import (
	"crypto/rand"
	"encoding/binary"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"time"
	signer "vSIS-Signature/Signer"
	prof "vSIS-Signature/prof"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/utils"
)

type signatureJSON struct {
	Message   []uint64   `json:"message"`
	X0        []uint64   `json:"x0"`
	X1        []uint64   `json:"x1"`
	Target    []uint64   `json:"target"`
	Signature [][]uint64 `json:"signature"`
}

func loadSignature(path string) (*signatureJSON, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var s signatureJSON
	return &s, json.Unmarshal(raw, &s)
}

func loadBmatrixCoeffs(path string) ([][]uint64, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var tmp struct {
		B [][]uint64 `json:"B"`
	}
	return tmp.B, json.Unmarshal(raw, &tmp)
}

// zeroPoly allocates an all-zero polynomial in NTT form.
func zeroPoly(r *ring.Ring) *ring.Poly { p := r.NewPoly(); return p }

// randomVector returns a slice of n fresh elements in [0,q).
// Panics on entropy failure.
func randomVector(n int, q uint64) []uint64 {
	if q == 0 {
		panic("randomVector: modulus q == 0")
	}

	// largest multiple of q that fits in 64 bits – values ≥bound are rejected
	var bound uint64 = ^uint64(0) - (^uint64(0) % q)

	vec := make([]uint64, n)
	for i := 0; i < n; {
		var buf [8]byte
		if _, err := rand.Read(buf[:]); err != nil {
			panic("randomVector: entropy read failed: " + err.Error())
		}
		v := binary.LittleEndian.Uint64(buf[:])
		if v >= bound { // rejection step: avoid modulo bias
			continue
		}
		vec[i] = v % q
		i++
	}
	return vec
}

// copyPoly returns a deep copy of p.
func copyPoly(r *ring.Ring, p *ring.Poly) *ring.Poly {
	out := r.NewPoly()
	ring.Copy(p, out)
	return out
}

// addInto  :  dst += src
func addInto(r *ring.Ring, dst, src *ring.Poly) { r.Add(dst, src, dst) }

// mulScalarNTT multiplies every coefficient of p by the scalar c (mod q)
// and writes the result into dst.  Both p and dst must be in NTT form;
// dst may alias p for in-place updates.
func mulScalarNTT(r *ring.Ring, p *ring.Poly, c uint64, dst *ring.Poly) {
	if c == 0 {
		// quick-zero
		for i := range dst.Coeffs[0] {
			dst.Coeffs[0][i] = 0
		}
		return
	}

	q := r.Modulus[0]
	c %= q

	// Montgomery representation: coefficient-wise multiplication
	for i := range p.Coeffs[0] {
		dst.Coeffs[0][i] = (p.Coeffs[0][i] * c) % q
	}
}

func BuildWitness(
	ringQ *ring.Ring,
	A [][]*ring.Poly,
	b1 []*ring.Poly,
	B0Const []*ring.Poly,
	B0Msg [][]*ring.Poly,
	B0Rnd [][]*ring.Poly,
	/* private */
	s []*ring.Poly,
	x1 *ring.Poly,
	u []*ring.Poly,
	x0 []*ring.Poly,
) (w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly) {
	defer prof.Track(time.Now(), "BuildWitness")

	n := len(A)    // #rows in A
	m := len(s)    // len(signature vector)
	lu := len(u)   // len(message block
	lx0 := len(x0) // len(mask block
	// witness vector (s, u, x0) has length m + lu + lx0
	k := m + lu + lx0

	// -------------------------------------------------------------------------
	// 0)  Sanity-check dimensions
	// -------------------------------------------------------------------------
	if len(b1) != n || len(B0Const) != n {

		log.Fatal("dimension mismatch in public vectors")
	}
	for _, row := range A {
		if len(row) != m {
			log.Fatal("wrong A row length ≠ m")

		}
	}

	// -------------------------------------------------------------------------
	// 1)  Verify the proof-friendly equation
	// -------------------------------------------------------------------------
	// (b1 ⊙ A)·s
	left1 := make([]*ring.Poly, n)
	for j := 0; j < n; j++ {
		left1[j] = zeroPoly(ringQ)
		for t := 0; t < m; t++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(b1[j], A[j][t], tmp) // b₁ⱼ * Aⱼ,t
			ringQ.MulCoeffs(tmp, s[t], tmp)
			addInto(ringQ, left1[j], tmp)
		}
	}

	// (A·s) * x1
	left2 := make([]*ring.Poly, n)
	for j := 0; j < n; j++ {
		left2[j] = zeroPoly(ringQ)
		for t := 0; t < m; t++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(A[j][t], s[t], tmp) // Aⱼ,t * s_t
			ringQ.MulCoeffs(tmp, x1, tmp)
			addInto(ringQ, left2[j], tmp)
		}
	}

	// B0(1;u;x0)
	right := make([]*ring.Poly, n)
	one := zeroPoly(ringQ)
	one.Coeffs[0][0] = 1 // constant 1
	for j := 0; j < n; j++ {
		right[j] = copyPoly(ringQ, B0Const[j]) // 1 · B0const

		// + message part
		for i := 0; i < lu; i++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(B0Msg[i][j], u[i], tmp)
			addInto(ringQ, right[j], tmp)
		}
		// + randomness part
		for i := 0; i < lx0; i++ {
			tmp := ringQ.NewPoly()
			ringQ.MulCoeffs(B0Rnd[i][j], x0[i], tmp)
			addInto(ringQ, right[j], tmp)
		}
	}

	// Check equality row by row
	for j := 0; j < n; j++ {
		tmp := ringQ.NewPoly()
		ringQ.Sub(left1[j], left2[j], tmp) // (b⊙A)s − (A s)x1
		ringQ.Sub(tmp, right[j], tmp)      // − B0(...)
		if !ringQ.Equal(tmp, ringQ.NewPoly()) {
			fmt.Printf("Want 0 got %d\n", tmp.Coeffs[0][j])
			log.Fatal("proof-friendly eq. fails on row", j)
		}
	}
	// -------------------------------------------------------------------------
	// Build the witnesses such as : w1 = (s, u, x0), w2 = x1, w3 = w1.w2= (w_{1,i}*x1)_i
	// -------------------------------------------------------------------------
	w1 = make([]*ring.Poly, k)
	for i := 0; i < m; i++ {
		w1[i] = copyPoly(ringQ, s[i]) // w1[i] = s_i
	}
	for i := 0; i < lu; i++ {
		w1[m+i] = copyPoly(ringQ, u[i]) // w1[m+i] = u_i
	}
	for i := 0; i < lx0; i++ {
		w1[m+lu+i] = copyPoly(ringQ, x0[i]) // w1[m+lu+i] = x0_i
	}
	w2 = copyPoly(ringQ, x1) // w2 = x1
	w3 = make([]*ring.Poly, k)
	for i := 0; i < k; i++ {
		w3[i] = ringQ.NewPoly()
		ringQ.MulCoeffs(w1[i], w2, w3[i]) // w3[i] = w1[i] * x1
	}
	// -------------------------------------------------------------------------
	// Return the witnesses
	return w1, w2, w3
}

func BuildWitnessFromDisk() (w1 []*ring.Poly, w2 *ring.Poly, w3 []*ring.Poly) {
	defer prof.Track(time.Now(), "BuildWitnessFromDisk")

	// ‣ 0. parameters ----------------------------------------------------------
	par, _ := signer.LoadParams("../Parameters/Parameters.json")

	ringQ, _ := ring.NewRing(par.N, []uint64{par.Q})

	// convenience: explicit in-place lift
	toNTT := func(p *ring.Poly) { ringQ.NTT(p, p) }

	// ‣ 1. matrix A  (already stored in NTT, **do not** touch) -----------------
	pk, _ := signer.LoadPublicKey("../public_key/public_key.json")

	A := make([][]*ring.Poly, 1)
	A[0] = make([]*ring.Poly, len(pk.A))
	for i, coeffs := range pk.A {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs) // pk.A is already NTT
		A[0][i] = p               // no toNTT
	}

	// ‣ 2. B-matrix columns  (stored in coefficient domain → lift) -------------
	Bcoeffs, _ := loadBmatrixCoeffs("../Parameters/Bmatrix.json")

	B0Const := make([]*ring.Poly, 1)
	B0Msg := [][]*ring.Poly{make([]*ring.Poly, 1)}
	B0Rnd := [][]*ring.Poly{make([]*ring.Poly, 1)}
	b1 := make([]*ring.Poly, 1)

	makePolyNTT := func(raw []uint64) *ring.Poly {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], raw)
		toNTT(p) // lift to NTT **once**
		return p
	}

	B0Const[0] = makePolyNTT(Bcoeffs[0])  //   B₀,const
	B0Msg[0][0] = makePolyNTT(Bcoeffs[1]) //   B₀,msg
	B0Rnd[0][0] = makePolyNTT(Bcoeffs[2]) //   B₀,rnd
	b1[0] = makePolyNTT(Bcoeffs[3])       //   b₁

	// ‣ 3. ρ  (compression vector) --------------------------------------------
	prng, _ := utils.NewPRNG()
	rho := make([]uint64, 1)
	rbuf := make([]byte, 8)
	prng.Read(rbuf)
	rho[0] = uint64(rbuf[0]) % par.Q // small is fine

	// ‣ 4. signature blobs -----------------------------------------------------
	sig, _ := loadSignature("../Signature/Signature.json")

	// message  u   and masks  x₀, x₁  – stored in coeff domain → lift
	m := makePolyNTT(sig.Message)
	x0 := makePolyNTT(sig.X0)
	x1 := makePolyNTT(sig.X1)

	// signature vector s  – already NTT (output of GaussSamp) ------------------
	s := make([]*ring.Poly, len(sig.Signature))
	for i, coeffs := range sig.Signature {
		p := ringQ.NewPoly()
		copy(p.Coeffs[0], coeffs) // keep as is (NTT)
		s[i] = p
	}

	// ‣ 5. build quadratic gate -----------------------------------------------
	w1, w2, w3 = BuildWitness(
		ringQ,
		A, b1,
		B0Const, B0Msg, B0Rnd,
		/*private*/ s, x1, []*ring.Poly{m}, []*ring.Poly{x0})

	return w1, w2, w3
}
