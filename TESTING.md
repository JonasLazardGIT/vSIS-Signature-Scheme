# PACS and LVCS Test Harness

This repository contains a deterministic simulation of the PACS argument with
LVCS row commitments. The test suite exercises both the happy path and a wide
set of negative cases that each trigger a specific verifier check.

## What is covered

| Tamper | Rejecting layer |
| ------ | --------------- |
| modify `bar` or `C` | LVCS linear map |
| head index in `E`, mismatched opening indices | LVCS binding / tail-only |
| flip opened `Pvals` or Merkle path bytes | DECS/Merkle integrity |
| change `Q`, `Γ'`, `γ'`, or mask `M` | Eq.(4) consistency |
| bump constant term of `Q` | ΣΩ invariant |
| exceed degree in `R_k` | LVCS degree bound |
| alter witness so `w3 != w1*w2` | `VerifyFullPACS` witness check |
| unbound `Sqs` row | commitment coupling (expected failure) |

## Running the tests

```bash
# run unit tests (including tamper cases)
go test ./PIOP -run TestPACS -count=1 -v

# full repo tests
go test ./... -race -count=1
```

## Parameters

Tests use `SimOpts` to vary the number of masked points `Ell`, the grid size
`Ncols`, batching factor `Rho`, and toggles such as `FSforC`, `TailOnly` and
`BindSqs`. See `PACS_Simulation_test.go` for examples (e.g., `TestPACSParamGrid`).

Deterministic challenges are sampled from a Fiat–Shamir RNG bound to
`{root, hash(R), C, Ω}` to ensure reproducible transcripts.
