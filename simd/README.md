
# Exploiting SIMDs for faster finite field arithmetic

Consider an SIMD environment where SIMD registers are divided into *n* 64-bit lanes, For the Intel AVX512 architecture there would be 8 such lanes. 
Consider each lane as being the 64-bit data registers available to an independent virtual machine. Now recall that the script generated 
finite field code is fully constant time. Therefore an identical stream of data-independent instructions can be executed by each lane, 
each processing different data.

Clearly all of these virtual machines can share the same program counter, stack register, and memory.

In this way more than one simultaneous but independent calculation can be performed, one in each lane. The speed up achievable in terms 
of throughput will be very close to the number of lanes.

Modern processors contain multiple cores, which are independent processors each equipped with their own SIMD capabilities, and their own 
independent memory resource. So multiple threads each capable of processing *n* virtual finite field calculators could be executing 
simultaneously.

Sounds good, but there are some limitations to this approach.

First while it can massively increase throughput, it does nothing for latency. Using SIMDs to improve latency is possible in settings where 
parallelism is plentiful and inherent to the algorithm being executed. And in our experience algorithms using finite field arithmetic 
like RFC7748 do have a limited amount of such parallelism, but not enough to come close to achieving a times *n* improvement in latency.

Second, a modern 64-bit architecture supports an integer multiplication instruction that generates a 128-bit product. SIMDs typically 
do not, and must be content to perform 32-bit multiplications that produce a 64-bit product that will fit inside of a lane. So basically 
our virtual machines are 32-bit processors and not 64-bit processors. This is a limitation that is not inherent to SIMDs, a better SIMD 
might in future remove this restriction (and the Intel AVX-IFMA architecture is already a step in that direction).

Third, the proposed method depends on the algorithm being constant time. In cryptography this is often regarded as being an essential 
requirement anyway, to avoid timing based side channel attackers. But in some contexts, like signature verification with a public key, 
it is not and non-constant time implementations are commonly used, and are faster. One simple response would be to implement in constant 
time anyway, even if not necessary from a security viewpoint.


# Experimental SIMD implementation

Here we provide versions of the scripts which implement finite field arithmetic using C language intrinsics which exploit the 
SSE/AVX2/AVX512/AVX-IFMA and 
NEON SIMD extensions supported by many Intel AMD and ARM processors. The SSE, AVX2/AVX512/AVX-IFMA and NEON instructions typically perform 
operations in parallel on a number of 64-bit lanes implemented in 128, 256-bit or 512-bit registers. This experimental implementation was 
built and tested on a standard x64 PC/laptop and a Raspberry Pi 400 computer which uses the Cortex-A72 processor.

SSE AVX-IFMA and NEON support two lanes. AVX2 supports four, and AVX512 supports 8. Hardware support for AVX-IFMA and AVX512
is currently quite patchy.

Recall that the script generated code is fully constant time. Therefore more than one simultaneous but independent calculation can be 
performed, one in each lane.

The main support is provided by the scripts *pseudo_simd.py* and *monty_simd.py*. Edit these scripts to choose NEON, SSE, AVX2 or AVX-512.
Dedicated scripts are provided for AVX-IFMA.

## RFC7748

RFC7748 defines a simple but powerful elliptic curve key exchange algorithm. For example it lies at the heart of the TLS and Signal
protocols. Here we provide a simple implementation.

## Quick Start on an Intel/AMD desktop/laptop PC

The setup is the same as for the original scripts. Make sure all dependencies are available. 
Generate the finite field arithmetic for an X25519 elliptic curve

	python3 pseudo_simd.py X25519

First run the *time* application which has been created by the script to get benchmark timings. Next drop the generated file *field.c* 
into *rfc7748_simd.c*. Compile as

	gcc -O3 -march=native rfc7748_simd.c -lcpucycles -o rfc7748_simd

When the program is run, a pair of test vectors from RFC7748 are calculated simultaneously, timings are taken, and a double key exchange 
is performed. Next generate the finite field arithmetic for the larger X448 curve 

	python3 monty_simd.py X448

Replace the *field.c* code in *rfc7748_simd.c* and compile and run the program again.

# Experimental SIMT implementation

Versions of the scripts are now provided that generate code suitable for Nvidia GPUs. Code is generated in CUDA which is
basically C/C++ for GPUs. GPUs allow massive throughput at the cost of latency. Each core implements a simple and rather slow 32-bit 
architecture -- but there are typically thousands of them. To use these scripts an Nvidia graphics card is required, plus the Nvidia 
CUDA software development toolkit -- in particular the nvcc compiler.  

## Quick Start

The setup is the same as for the original scripts. Make sure all dependencies are available. 
Generate the finite field arithmetic for an X25519 elliptic curve

	python3 pseudo_cuda.py X25519

First run the *time* application which has been created by the script to get benchmark timings. Next drop the generated file *field.cu* 
into *rfc7748_simt.cu*. Compile as

	nvcc -O2 rfc7748_simt.cu -lcpucycles -o rfc7748_simt

When the program is run, a pair of test vectors from RFC7748 are calculated simultaneously, timings are taken, and a double key exchange is 
performed. Next generate the finite field arithmetic for the larger X448 curve 

	python3 monty_cuda.py X448

Replace the *field.cu* code in *rfc7748_simt.cu* and compile and run the program again.

# Application

Two, four or eight or more protocol executions can be carried out simultaneously, which can be somewhat faster that executing them serially, 
one after another. For example a TLS server can calculate both its own public key and a shared secret at the same time.

