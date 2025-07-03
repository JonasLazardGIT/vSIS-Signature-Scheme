package decs

import (
	"bytes"
	"crypto/sha256"
)

// MerkleTree is a full binary Merkle tree of 32-byte hashes.
type MerkleTree struct {
	layers [][][32]byte
}

// BuildMerkleTree builds a balanced tree from leaves.
func BuildMerkleTree(leaves [][]byte) *MerkleTree {
	n := len(leaves)
	size := 1
	for size < n {
		size <<= 1
	}
	layer := make([][32]byte, size)
	for i := 0; i < n; i++ {
		layer[i] = sha256.Sum256(leaves[i])
	}
	for i := n; i < size; i++ {
		layer[i] = sha256.Sum256(nil)
	}
	layers := [][][32]byte{layer}

	for sz := size; sz > 1; sz >>= 1 {
		prev := layers[len(layers)-1]
		next := make([][32]byte, sz/2)
		for i := 0; i < sz; i += 2 {
			h := sha256.New()
			h.Write(prev[i][:])
			h.Write(prev[i+1][:])
			next[i/2] = sha256.Sum256(h.Sum(nil))
		}
		layers = append(layers, next)
	}

	return &MerkleTree{layers: layers}
}

// Root returns the root hash.
func (mt *MerkleTree) Root() [32]byte {
	return mt.layers[len(mt.layers)-1][0]
}

// Path returns the sibling path for leaf idx.
func (mt *MerkleTree) Path(idx int) [][]byte {
	path := make([][]byte, len(mt.layers)-1)
	for lvl := 0; lvl < len(mt.layers)-1; lvl++ {
		sib := idx ^ 1
		h := mt.layers[lvl][sib]
		path[lvl] = h[:]
		idx >>= 1
	}
	return path
}

// VerifyPath checks leaf→root via path.
func VerifyPath(leaf []byte, path [][]byte, root [32]byte, idx int) bool {
	h := sha256.Sum256(leaf)
	for _, sib := range path {
		tmp := sha256.New()
		if idx&1 == 0 {
			tmp.Write(h[:])
			tmp.Write(sib)
		} else {
			tmp.Write(sib)
			tmp.Write(h[:])
		}
		h = sha256.Sum256(tmp.Sum(nil))
		idx >>= 1
	}
	return bytes.Equal(h[:], root[:])
}
