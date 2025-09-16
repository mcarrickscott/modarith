# Experimental SIMD implementation

Here we provide versions of the scripts which implement finite field arithmetic using C language intrinsics which exploit the SSE/AVX and 
NEON SIMD extensions supported by many Intel AMD and ARM processors. The SSE and NEON instructions typically perform operations in parallel 
on a pair of 64-bit lanes implemented in 128-bit registers. This experimental implementation was built and tested on a standard x64 
PC/laptop and a Rasberry Pi 400 computer which uses the Cortex-A72 processor.

Recall that the script generated code is fully constant time. Therefore two simultaneous but independent calculations can be performed,
one in each lane.

# RFC7748

RFC7748 defines a simple but powerful elliptic curve key exchange algorithm. For example it lies at the heart of the TLS and Signal
protocols. Here we provide a simple implementation.

# Quick Start on an Intel/AMD desktop/laptop PC

The setup is the same as for the original scripts. Make sure all dependencies are available. 
Generate the finite field arithmetic for an X25519 elliptic curve

	python3 pseudo_sse.py X25519

Now drop the generated file *field.c* into *rfc7748_simd.c*. Compile as

	gcc -O3 -march=native rfc7748_simd.c -lcpucycles -o rfc7748_simd

When the program is run, a pair of test vectors from RFC7748 are calculated simultaneously, timings are taken, and a double key exchange 
is performed. Next generate the finite field arithmetic for the larger X448 curve 

	python3 monty_sse.py X448

Replace the *field.c* code in *rfc7748_simd.c* and compile and run the program again.

# Application

Two algorithm executions can be carried out simultaneously, which can be somewhat faster that executing them serially, one after another.
For example a TLS server can calculate both its own public key and a shared secret at the same time.

