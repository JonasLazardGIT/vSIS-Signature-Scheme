# Running and Testing vSIS-Signature

This guide explains how to execute the command-line tools and tests in the vSIS-Signature repository. It covers available flags, environment variables, and package-specific options, with a focus on the PIOP proof system.

## Prerequisites

* Go 1.23 or newer (tested with Go 1.24).
* Required Go modules are pulled automatically by `go test` or `go run`.
* The repository expects a POSIX-like shell.

---

## Running the Main Program

### Basic execution

```
go run ./cmd
```

This command performs three steps sequentially:
1. Generate system parameters and common random data.
2. Generate a trapdoor keypair and produce one signature.
3. Verify the signature.

The main entry point is `cmd/cmd.go`【F:cmd/cmd.go†L17-L33】.

### Generated files

Running the program creates or reuses several folders at the repository root:

| Folder / File | Purpose |
|---------------|---------|
| `Parameters/Parameters.json` | Serialized public parameters such as `N`, `Q`, and gadget data【F:System/System.go†L72-L83】 |
| `Parameters/Bmatrix.json` | Common random B-matrix used by the hash【F:System/System.go†L103-L119】 |
| `public_key/public_key.json` & `private_key/private_key.json` | Trapdoor keypair saved by `Sign()` if not already present【F:Signer/signer.go†L66-L87】【F:Signer/signer.go†L220-L258】 |
| `Signature/Signature.json` | Message, randomness, target syndrome, and signature bundle【F:Signer/signer.go†L160-L166】 |

`Verify()` consumes these files during the final check【F:Verifier/verifier.go†L36-L70】.

### Measuring object sizes

Setting the environment variable `MEASURE_SIZES=1` enables additional size accounting logs:

```
MEASURE_SIZES=1 go run ./cmd
```

The `measure` package activates when this variable equals `1`【F:measure/measure.go†L10-L17】.

### Distribution analysis build

The `cmd` directory also contains an instrumented program for coefficient distribution analysis guarded by the `analysis` build tag. Run it as:

```
go run -tags analysis ./cmd
```

This version generates 50 signatures, stores the centred coefficients, and writes
`coefficient_stats.json`, `coefficient_distribution.png`, and
`coefficient_distribution.html` to `Read_signatures/`【F:cmd/distribution_analysis.go†L418-L458】.
Adjust the number of runs or the destination folder by editing the `runs` and
`outDir` constants in `distribution_analysis.go`.

### Preimage sampler configuration

The pre-image sampling routines react to several environment variables:

| Variable | Effect |
|----------|--------|
| `S_CONST` | Override the spectral-bound constant used in `SpectralBound`【F:Preimage_Sampler/Preimage_Sampling.go†L18-L28】 |
| `SIGMA_CONST` | Override the smoothing parameter in `CalculateParams`【F:Preimage_Sampler/Preimage_Sampling.go†L38-L47】 |
| `ALPHA_MODE` | Set to `GEN` to use α=(t+1)σ; otherwise α=2σ in `GaussSamp`【F:Preimage_Sampler/Preimage_Sampling.go†L106-L113】 |
| `CENTER_SYNDROME` | Set to `1` to center the syndrome before G-sampling【F:Preimage_Sampler/Preimage_Sampling.go†L132-L135】 |
| `RANDOMIZED_ROUNDING` | Set to `0` to disable randomized rounding in perturbation sampling【F:Preimage_Sampler/Perturbation_Sampling.go†L426-L431】 |

Additional debug flags enable verbose logs in other packages:

* `DEBUG_DECS` – verbose DECS verifier output【F:DECS/decs_verifier.go†L11-L12】
* `DEBUG_LVCS` – verbose LVCS verifier output【F:LVCS/lvcs_verifier.go†L14-L15】

---

## Testing

### Common `go test` flags

* `-run <regex>` – run only tests whose names match the regular expression.
* `-count=1` – disable the test cache to force fresh execution.
* `-v` – verbose output, listing each test name.
* `-race` – enable the race detector (slower, but finds data races).
* `-bench <regex>` – run benchmarks matching the expression.
* `-benchmem` – print benchmark memory allocations.

Example of a full repository test run:

```
go test ./... -race -count=1
```

### PIOP package

The PIOP package contains the most comprehensive test harness, simulating and tampering with the PACS argument. Typical invocation:

```
go test ./PIOP -run TestPACS -count=1 -v
```

This runs the main simulation plus tamper cases, printing detailed timings and checks【F:TESTING.md†L22-L27】. The simulator is configurable through the `SimOpts` structure defined in `PACS_Simulation_test.go`:

| Field    | Meaning                                            |
|----------|----------------------------------------------------|
| `Ell`    | Number of masked points per row (ℓ). Defaults to 1. |
| `Ncols`  | Size of the evaluation grid Ω. Defaults to 8.       |
| `Rho`    | Batching factor. Defaults to 1.                     |
| `FSforC` | Derive coefficient matrix `C` via Fiat–Shamir if `true`.|
| `TailOnly` | Restrict openings to the masked tail.             |
| `BindSqs` | Bind `Sqs` into the same commitment.               |

Definition excerpt【F:PIOP/PACS_Simulation_test.go†L73-L80】.

To sweep parameters or inject faults, edit the test cases in `TestPACSParamGrid` or add new cases using these options.

### Other packages

Run tests package by package as needed:

```
# Differential equation commitment scheme
go test ./DECS

# Linear verifiable commitment scheme
go test ./LVCS

# Polynomial commitment scheme
go test ./PCS

# Preimage sampler
go test ./Preimage_Sampler

# Helper utilities for size measurement
go test ./measure
```

Packages without `_test.go` files (e.g. `Signer`, `Verifier`, `System`, `vSIS-HASH`) report "[no test files]" when tested.

### Full suite

Running `go test ./...` executes all available package tests. The current suite completes successfully【676805†L1-L8】.

---

## Additional Notes

* The test harness writes intermediate files (e.g., parameters, keys, signatures) under the repository root. Clean them manually if needed.
* Use `GOMAXPROCS` to limit CPU parallelism when running heavy tests on constrained machines:
  `GOMAXPROCS=2 go test ./PIOP`.
* Benchmarks can be run with `go test -bench .` in any package containing benchmark functions.

This documentation should provide a comprehensive reference for executing and customizing the main program and test suites within the repository.
