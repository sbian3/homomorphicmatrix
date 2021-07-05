# Homomorphic Linear Transformation
This is an example implementation of homomorphic linear transformation over
RLWE-based ciphertexts based on the [SEAL library](https://github.com/microsoft/SEAL "SEAL").

## Requirements
- cmake (>= 3.16.3) 
- C++ compiler with C++17 compatibility

## Quick Start
Please also refer to README.md for building the SEAL library. Here, we only
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

## CMake Options

| CMake option | Values | Information |
| --- | --- | --- |
| HLT_FETCH_THIRDPARTY | **ON**/OFF | Automatically download and build dependecies
| HLT_BUILD_TEST | ON/**OFF** | Build test to make sure our routines works

## Directory Layout
- Benchmarks
  - bench_hlt: Benchmark primitives
- Custom library routines
  - src/util: Direct and packed convolution implementation
