# Homomorphic Linear Transformation
This is an example implementation of homomorphic linear transformation over
RLWE-based ciphertexts based on the SEAL library.

## Requirements
Note that on linux, a cmake of version >= 3.16.3 ,a SEAL library version >= 3.6.0 and C++ compiler with C++17 features

## Quick Start
Please also refer to READSEAL.md for building the SEAL library. Here, we only
provide a simplified version of the build process on Ubuntu 20.04.2 LTS.

First, build SEAL by running

```sh
cmake -S . -B build
cmake --build build
```

The binary for running the linear transformation will be in build/bin
("direct_conv" for direct convolution, "packed_conv" for packed convolution,
and "general_lt" for general linear transformation). To run the benchmark,
simply execute the binary with proper arguments, e.g.,
```sh
./direct_conv 1024 9 2048
```
for running a direct convolution with input dimension 1024, filter dimension
 3x3=9, and lattice dimension 2048.

 You can also run 
```sh
./hlt_bench.sh
```
to benchmark a set of pre-defined computations, where the benchmarked results
will be in the folder hlt_result/


## Directory Layout
- Benchmarks
  - bench_hlt: Benchmark primitives
- Custom library routines
  - src/util: Direct and packed convolution
