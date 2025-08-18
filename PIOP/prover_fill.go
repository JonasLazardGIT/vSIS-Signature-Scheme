package PIOP

import (
	"math/big"
	"time"

	"github.com/tuneinsight/lattigo/v4/ring"
	prof "vSIS-Signature/prof"
)

// GlobSlack holds global slack digits and their bit decomposition.
type GlobSlack struct {
	D     []*ring.Poly   // Δℓ const polys (Δℓ*S0inv)
	DBits [][]*ring.Poly // Δℓ bits
}

// appendGlobalSlack allocates slack digit columns.
func appendGlobalSlack(r *ring.Ring, LS, W int) GlobSlack {
	defer prof.Track(time.Now(), "appendGlobalSlack")
	mk := func() *ring.Poly { p := r.NewPoly(); r.NTT(p, p); return p }
	s := GlobSlack{
		D:     make([]*ring.Poly, LS),
		DBits: make([][]*ring.Poly, LS),
	}
	for l := 0; l < LS; l++ {
		s.D[l] = mk()
		s.DBits[l] = make([]*ring.Poly, W)
		for u := 0; u < W; u++ {
			s.DBits[l][u] = mk()
		}
	}
	return s
}

// ProverFillIntegerL2 populates all limb digits, remainder chains, slack and carries.
func ProverFillIntegerL2(
	r *ring.Ring, w1 []*ring.Poly, mSig int,
	spec BoundSpec,
	cols DecompCols, glob GlobCarry, slack GlobSlack,
	S0, S0inv uint64,
) error {
	defer prof.Track(time.Now(), "ProverFillIntegerL2")
	q := r.Modulus[0]
	N := r.N
	R := spec.R
	LS := spec.LS
	W := spec.W

	coeffW1 := make([]*ring.Poly, mSig)
	for k := 0; k < mSig; k++ {
		coeffW1[k] = r.NewPoly()
		r.InvNTT(w1[k], coeffW1[k])
	}

	SRow := make([]*big.Int, N)
	Drows := make([][]uint64, LS)
	for l := 0; l < LS; l++ {
		Drows[l] = make([]uint64, N)
	}

	bq := new(big.Int).SetUint64(q)
	halfq := new(big.Int).Rsh(bq, 1)
	Rbig := new(big.Int).SetUint64(R)

	for j := 0; j < N; j++ {
		sum := new(big.Int)
		for k := 0; k < mSig; k++ {
			a := new(big.Int).SetUint64(coeffW1[k].Coeffs[0][j])
			if a.Cmp(halfq) > 0 {
				a.Sub(a, bq)
			}
			a.Mul(a, a)
			sum.Add(sum, a)
		}
		SRow[j] = sum
		tmp := new(big.Int).Set(sum)
		for l := 0; l < LS; l++ {
			if tmp.Sign() == 0 {
				Drows[l][j] = 0
				continue
			}
			rem := new(big.Int)
			tmp.DivMod(tmp, Rbig, rem)
			Drows[l][j] = rem.Uint64()
		}
	}

	for l := 0; l < LS; l++ {
		fillRowUint64(r, cols.D[l], Drows[l], q)
		for u := 0; u < W; u++ {
			Bu := make([]uint64, N)
			for j := 0; j < N; j++ {
				Bu[j] = (Drows[l][j] >> uint(u)) & 1
			}
			fillRowUint64(r, cols.Bit[l][u], Bu, q)
		}
	}

	Trows := make([][]uint64, LS+1)
	for l := 0; l <= LS; l++ {
		Trows[l] = make([]uint64, N)
	}
	for j := 0; j < N; j++ {
		t := new(big.Int).Set(SRow[j])
		Trows[0][j] = modU64(t, q)
		for l := 0; l < LS; l++ {
			t.Sub(t, new(big.Int).SetUint64(Drows[l][j]))
			t.Div(t, Rbig)
			Trows[l+1][j] = modU64(t, q)
		}
	}
	for l := 0; l <= LS; l++ {
		fillRowUint64(r, cols.T[l], Trows[l], q)
	}

	SumDigits := make([]*big.Int, LS)
	totalSum := new(big.Int)
	for l := 0; l < LS; l++ {
		SumDigits[l] = new(big.Int)
		for j := 0; j < N; j++ {
			SumDigits[l].Add(SumDigits[l], new(big.Int).SetUint64(Drows[l][j]))
		}
	}
	for j := 0; j < N; j++ {
		totalSum.Add(totalSum, SRow[j])
	}

	B := limbsToBig(spec.Beta2, R)
	Delta := new(big.Int).Sub(B, totalSum)
	if Delta.Sign() < 0 {
		Delta.SetInt64(0)
	}

	DeltaLimbs := splitBaseR(Delta, R, LS)

	carries := make([]*big.Int, LS+1)
	carries[0] = new(big.Int)
	for l := 0; l < LS; l++ {
		lhs := new(big.Int).Add(SumDigits[l], carries[l])
		lhs.Add(lhs, new(big.Int).SetUint64(DeltaLimbs[l]))
		lhs.Sub(lhs, new(big.Int).SetUint64(spec.Beta2[l]))
		cnext := new(big.Int).Div(lhs, Rbig)
		carries[l+1] = cnext
	}

	fillSlackDigits(r, slack, DeltaLimbs, S0inv, q, R, W)
	fillCarries(r, glob, carries, S0inv, q, R)
	return nil
}

func fillRowUint64(r *ring.Ring, p *ring.Poly, vals []uint64, q uint64) {
	for j := 0; j < len(vals); j++ {
		p.Coeffs[0][j] = vals[j] % q
	}
	for j := len(vals); j < len(p.Coeffs[0]); j++ {
		p.Coeffs[0][j] = 0
	}
	r.NTT(p, p)
}

func fillConstUint64(r *ring.Ring, p *ring.Poly, val uint64, q uint64) {
	for j := range p.Coeffs[0] {
		p.Coeffs[0][j] = val % q
	}
	r.NTT(p, p)
}

func modU64(x *big.Int, q uint64) uint64 {
	mod := new(big.Int).Mod(x, new(big.Int).SetUint64(q))
	return mod.Uint64()
}

func splitBaseR(x *big.Int, R uint64, L int) []uint64 {
	out := make([]uint64, L)
	tmp := new(big.Int).Set(x)
	Rb := new(big.Int).SetUint64(R)
	zero := new(big.Int)
	for i := 0; i < L; i++ {
		if tmp.Cmp(zero) == 0 {
			out[i] = 0
			continue
		}
		rem := new(big.Int)
		tmp.DivMod(tmp, Rb, rem)
		out[i] = rem.Uint64()
	}
	return out
}

func limbsToBig(limbs []uint64, R uint64) *big.Int {
	acc := new(big.Int)
	Rb := new(big.Int).SetUint64(R)
	pow := new(big.Int).SetUint64(1)
	for _, d := range limbs {
		term := new(big.Int).Mul(new(big.Int).SetUint64(d), pow)
		acc.Add(acc, term)
		pow.Mul(pow, Rb)
	}
	return acc
}

func fillSlackDigits(r *ring.Ring, slack GlobSlack, deltaLimbs []uint64, S0inv, q, R uint64, W int) {
	for l := 0; l < len(deltaLimbs); l++ {
		val := (deltaLimbs[l] % q) * (S0inv % q) % q
		fillConstUint64(r, slack.D[l], val, q)
		for u := 0; u < W; u++ {
			bit := (deltaLimbs[l] >> uint(u)) & 1
			fillConstUint64(r, slack.DBits[l][u], bit, q)
		}
	}
}

func fillCarries(r *ring.Ring, glob GlobCarry, carries []*big.Int, S0inv, q, R uint64) {
	for l := 0; l < len(carries); l++ {
		val := new(big.Int).Mul(carries[l], new(big.Int).SetUint64(S0inv))
		val.Mod(val, new(big.Int).SetUint64(q))
		fillConstUint64(r, glob.C[l], val.Uint64(), q)
		for u := 0; u < len(glob.CBits[l]); u++ {
			mask := new(big.Int).Lsh(big.NewInt(1), uint(u))
			b := new(big.Int).And(carries[l], mask)
			bit := uint64(0)
			if b.Cmp(big.NewInt(0)) != 0 {
				bit = 1
			}
			fillConstUint64(r, glob.CBits[l][u], bit, q)
		}
	}
}
