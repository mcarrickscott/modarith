# Modarith - python scripts for generation of code for multi-precision modular arithmetic

Script generated finite field arithmetic for elliptic curve cryptography

This repository contains Python 3 scripts for the automatic generation of efficent code for multi-precision modular arithmetic, for 16, 32 and 64-bit architectures, in C and in Rust. 

The code uses a multi-limb unsaturated-radix representation for big numbers.

For both languages there are two scripts, one specialised for pseudo-Mersenne moduli of the form $2^m-c$ (*pseudo.py* and *pseudo_rust.py*), and one that generates efficient and tailored code for other prime shapes, in particular generalised Mersennes like those required for the standardised NIST256, NIST384, and C448 curves, and Montgomery friendly moduli, as sometimes recommended for use with isogenies (*monty.py* and *monty_rust.py*).

The code generated includes functions for modular addition, subtraction, multiplication, inversion, quadratic residuosity and square roots. In other words all of the requirements to implement field arithmetic in the context of elliptic curve cryptography.

# Dependencies

As dependencies it is required that an appropriate C compiler (*gcc/clang/icx*) and/or the Rust compiler (*rustc*) is included in the path. For the C code also ensure that the *clang-format* and *cppcheck* utilities are installed.

Also in the path must be the very useful utility addchain 
	
	https://github.com/mmcloughlin/addchain

For accurate timings across a range of architectures for the C code, install Dan Bernstein's libcpucycles utility from https://cpucycles.cr.yp.to/ . It may be necessary to run *ldconfig*.

# Quick start

For a quick start copy the files from here into a working directory, and try

	python3 pseudo.py 64 2**255-19
	./time

Then 64-bit code for the suggested modulus is generated and tested. An executable that times important functions will be created if the platform allows it. The standalone C timing code is output to *time.c*, and code for production use is output to *code.c* and *header.h*

For Rust 

	python3 pseudo_rust.py 64 2**255-19
	./time

In this case the output is directed to files *time.rs* and *code.rs*

For more details read the comments in the provided scripts.

As a Proof of Concept, elliptic curve code for RFC7748 is provided in the files *rfc7748.c* and *rfc7748.rs*

Read comments in these files for simple build instructions.

(RFC7748 describes an implementation of Diffie-Hellman key exchange on the Montgomery elliptic curves C25519 and C448. You can create your own Montgomery curve using the sagemath script provided in the file bowe.sage)

# Modular arithmetic API

Assume the modulus is $p$. The provided functions are

*nres()*   -- Convert a big number to internal format  \
*redc()*   -- Convert back from internal format, result $\lt p$  \
*modfsb()* -- Perform final subtraction to reduce from $\lt 2p$ to $\lt p$  \
*modadd()* -- Modular addition, result $\lt 2p$  \
*modsub()* -- Modular subtraction, result $\lt 2p$  \
*modneg()* -- Modular negation, result $\lt 2p$  \
*modmul()* -- Modular multiplication, result $\lt 2p$  \
*modsqr()* -- Modular squaring, result $\lt 2p$  \
*modmli()* -- Modular multiplication by a small integer, result $\lt 2p$ \
*modcpy()* -- Copy a big number \
*modpro()* -- Calculate progenitor, for subsequent use for modular inverses and square roots  \
*modinv()* -- Modular inversion  \
*modsqrt()* -- Modular square root \
*modis1()* -- Test for equal to unity  \
*modis0()* -- Test for equal to zero  \
*modone()* -- Set equal to unity  \
*modzer()* -- Set equal to zero  \
*modint()* -- Convert an integer to internal format  \
*modqr()*  -- Test for quadratic residue  \
*modcmv()* -- Conditional constant time move  \
*modcsw()* -- Conditional constant time swap  \
*modshl()* -- shift left by bits  \
*modshr()* -- shift right by bits  \
*modexp()* -- export from internal format to byte array \
*modimp()* -- import to internal format from byte array  \
*modsign()* -- Extract sign (parity bit) \
*modcmp()* -- Test for equality \

## Using the scripts

The scripts can be used out-of-the-box using simple command-line arguments, as in the examples above.

The scripts can also be tailored in various ways at the top of the script, and in the ``user editable area''.

There are default settings for the choice of compiler, choice of using a clock cycle counter, encouragement for inlining certain functions, code formatting, the use of Karatsuba for (maybe) faster multiplication, and, for the C code,  an option to ``decorate'' function names to avoid name clashes.

New named moduli can also be provided in the user editable area, and some settings (like radix choice) adapted individually.

Function name decoration may be required to avoid name clashes in C. If using C++ namespaces can be used to avoid this necessity. It is not an issue for Rust.

## Constant time

All generated functions are written with the expectation that they will execute in constant time. But high level code is nevertheless at the mercy of both the compiler and the architecture.

It is strongly recommended that the generated assembly language be closely studied to ensure that there are no compiler introduced timing leaks. In particular code generated for the functions *modcmv* and *modcsw* should be checked, bearing in mind that they may be inlined by the compiler. If necessary compiler-specific measures should be taken to prevent inlining, and/or place these functions into a separately compiled module.   


# New - automatic generation of Edwards or Weierstrass curve API code in C

1. Decide on wordlength (32 or 64 bit)
2. Choose the field prime.
3. Automatically generate field code to files *code.c* and *header.h*. If (pseudo)-Mersenne prime, use script *pseudo.py*, else use *monty.py*
4. Decide on Edwards or Weierstrass curve. The number of points on the curve is assumed to be a prime for Weierstrass, or 4 or 8 times a prime for Edwards. This prime is the order of the group.
5. Drop *code.c* into *edwards.c* or *weierstrass.c* and provide some curve constants (a constant *B* or *d* and a group generator point *x, y*) where indicated (some are provided already). If larger constants are required, a tool *make.py* is provided to generate them
6. Insert the prime group order into *testcurve.c* where indicated (several examples are there already)
7. Make sure fixed API header file *curve.h* is in the path. Note that this API is independent of the curve, its associated field and its parameters.
8. Compile and link *testcurve.c* with *edwards.c* or *weierstrass.c*
9. Run testcurve to test the arithmetic and perform some timings.

The API interface is as indicated in *curve.h*. The API is completely implemented in *edwards.c* or *weierstrass.c*

## Quickstart 1:-

	python pseudo.py 64 Ed25519
Drop code.c into edwards.c where indicated

	gcc -O2 testcurve.c edwards.c -lcpucycles -o testcurve
	./testcurve

Note that this intermediate API only provides the elliptic curve functionality. A higher level algorithm API (like that provided for Ed448 signature) would use this API while itself providing additional algorithm specific random number and hashing functionality. It may also use the *monty.py* script to generate code to perform arithmetic modulo the prime group order, if so required by the algorithm.

## Quickstart 2:-

	python monty.py 64 0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab5844f3
Drop code.c into Ed448.c where indicated

	python monty.py 64 Ed448
Drop code.c into edwards.c where indicated

	gcc -O2 Ed448.c edwards.c -o Ed448
	./Ed448
