# Experimental SIMD implementation

Here we provide versions of the scripts which implement finite field arithmetic using C language intrinsics which exploit the 
SSE/AVX2/AVX512/AVX-IFMA and 
NEON SIMD extensions supported by many Intel AMD and ARM processors. The SSE, AVX2/AVX512/AVX-IFMA and NEON instructions typically perform 
operations in parallel on a number of 64-bit lanes implemented in 128, 256-bit or 512-bit registers. This experimental implementation was 
built and tested on a standard x64 PC/laptop and a Raspberry Pi 400 computer which uses the Cortex-A72 processor.

SSE AVX-IFMA and NEON support two lanes. AVX2 supports four, and AVX512 supports 8. Support for AVX-IFMA and AVX512
is currently quite patchy.

Recall that the script generated code is fully constant time. Therefore more than one simultaneous but independent calculation can be 
performed, one in each lane.

# RFC7748

RFC7748 defines a simple but powerful elliptic curve key exchange algorithm. For example it lies at the heart of the TLS and Signal
protocols. Here we provide a simple implementation.

# Quick Start on an Intel/AMD desktop/laptop PC

The setup is the same as for the original scripts. Make sure all dependencies are available. 
Generate the finite field arithmetic for an X25519 elliptic curve

	python3 pseudo_sse.py X25519

First run the *time* application which has been created by the script to get benchmark timings. Next drop the generated file *field.c* 
into *rfc7748_simd.c*. Compile as

	gcc -O3 -march=native rfc7748_simd.c -lcpucycles -o rfc7748_simd

When the program is run, a pair of test vectors from RFC7748 are calculated simultaneously, timings are taken, and a double key exchange 
is performed. Next generate the finite field arithmetic for the larger X448 curve 

	python3 monty_sse.py X448

Replace the *field.c* code in *rfc7748_simd.c* and compile and run the program again.

# Application

Two, four or eight protocol executions can be carried out simultaneously, which can be somewhat faster that executing them serially, one 
after another. For example a TLS server can calculate both its own public key and a shared secret at the same time.

