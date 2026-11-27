# FHE-C

**A C implementation and experimental optimization framework for Fully Homomorphic Encryption (FHE).**

FHE-C is a C-based project exploring efficient implementations of FHE-related computations, with a focus on low-level optimization, architecture-specific compilation, and performance evaluation across x86 and ARM platforms.

The project includes a modular source structure, configurable build system, testing infrastructure, and tools for profiling and simulation.

---

## Features

- **C implementation** — Low-level implementation designed for performance and experimentation.
- **Multi-architecture support** — Build configurations for x86 and ARM.
- **Architecture-specific optimizations** — Support for SIMD-oriented compilation and hardware-specific instruction sets.
- **Configurable tests** — Select different test programs through the Makefile.
- **Performance analysis** — Integration with Intel VTune and Valgrind.
- **ARM simulation** — Support for running experiments with gem5.
- **Cloud workflow** — Helper targets for transferring files to and connecting to cloud instances.

---

## Repository Structure

```text
FHE-C/
├── gem5/       # gem5 simulation configurations
├── src/        # Core implementation
├── test/       # Test programs
├── utils/      # Utility tools
├── Makefile    # Build, execution and profiling workflows
├── ifft.py    # Python utility
├── ifft_tp.py # Python utility
├── omega5h.py # Python utility
└── roots.py   # Python utility
```

---

## Requirements

### General

- GCC
- GNU Make
- Linux environment recommended

### x86

The default configuration uses:

- AVX2
- FMA
- Native architecture optimization (`-march=native`)

Your processor must support the instruction sets enabled by the compiler flags.

### ARM

The ARM configuration uses:

- AArch64 cross-compilation toolchain (`aarch64-linux-gnu-gcc`)
- ARMv8.2-A
- SIMD and cryptographic instruction extensions

---

## Getting Started

### 1. Clone the repository

```bash
git clone https://github.com/Joao-Pedro-Cabral/FHE-C.git
cd FHE-C
```

### 2. Build and run

The default configuration targets x86 and runs the default test:

```bash
make
```

By default, the Makefile selects:

- Architecture: `x86`
- Test: `ckks_test`

### 3. Build without running

```bash
make build
```

### 4. Select a test

Tests can be selected through the `TEST` variable:

```bash
make TEST=ckks_test
```

Replace `ckks_test` with the name of another available test in the `test/` directory.

---

## Build Configurations

### x86

The default target:

```bash
make ARCH=x86
```

Uses GCC with optimization flags including:

```text
-O3
-mavx2
-mfma
-march=native
```

### ARM

To build for ARM:

```bash
make ARCH=arm
```

The build uses the AArch64 cross-compiler by default:

```text
aarch64-linux-gnu-gcc
```

You can override the compiler:

```bash
make ARCH=arm ARM_CC=<your-cross-compiler>
```

---

## Debugging

Compile and run with debug configuration:

```bash
make debug
```

This enables debug-related compiler options and symbols.

---

## Memory Analysis

Run the selected test under Valgrind:

```bash
make valgrind
```

The target enables debug compilation and performs memory leak analysis.

---

## Performance Profiling

### Intel VTune

To profile hotspots with Intel VTune:

```bash
make vtune
```

This target collects hotspot profiling data using VTune and generates a result directory.

---

## ARM Simulation with gem5

FHE-C includes a workflow for running experiments using gem5.

The Makefile expects a gem5 ARM binary and configuration script.

Configure the gem5 path:

```bash
export GEM5_PATH=/path/to/gem5
```

Then run:

```bash
make gem5
```

The gem5 executable and configuration can also be overridden:

```bash
make gem5 GEM5=/path/to/gem5.opt GEM5_CONFIG=/path/to/config.py
```

---

## Cloud Workflow

The Makefile includes helper targets for interacting with a configured Google Cloud VM.

### SSH

```bash
make cloud_ssh
```

### Copy project files

```bash
make cloud_scp
```

These targets assume that Google Cloud CLI is installed and configured, and that the target instance is available.

---

## Cleaning

Remove generated build directories:

```bash
make clean
```

---

## Contributing

Contributions, suggestions, bug reports, and improvements are welcome.

If you'd like to contribute:

1. Fork the repository.
2. Create a feature branch.
3. Implement and test your changes.
4. Submit a pull request.

---

## License

Please refer to the repository for licensing information.