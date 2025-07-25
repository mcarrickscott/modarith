# Modarith - python scripts for generation of code for finite field arithmetic

Script generated finite field arithmetic for "fast enough" elliptic curve cryptography. See https://eprint.iacr.org/2024/779

This repository contains Python 3 scripts for the automatic generation of efficent code for multi-precision modular arithmetic, for 16, 32 and 64-bit architectures, in C and in Rust. 

The code uses a multi-limb unsaturated-radix representation for big numbers.

For both languages there are two scripts, one specialised for pseudo-Mersenne moduli of the form $2^m-c$ (*pseudo.py* and *pseudo_rust.py*), and one that generates efficient and tailored code for other prime shapes, in particular generalised Mersennes like those required for the standardised NIST256, NIST384, and X448 curves, and Montgomery friendly moduli, as sometimes recommended for use with isogenies (*monty.py* and *monty_rust.py*).

The code generated includes functions for modular addition, subtraction, multiplication, inversion, quadratic residuosity and square roots. In other words all of the requirements to implement field arithmetic in the context of elliptic curve cryptography.

The code produced by these scripts was used to create the TLSECC library of C and Rust functions which implement all of the elliptic curve cryptography required by the TLS1.3 protocol. See https://github.com/mcarrickscott/TLSECC

# Dependencies

As dependencies it is required that an appropriate C compiler (*gcc/clang/icx*) and/or the Rust compiler (*rustc*) is included in the path. For the C code also ensure that the *clang-format* and *cppcheck* utilities are installed.

Also in the path must be the very useful utility addchain 
	
	https://github.com/mmcloughlin/addchain

For accurate timings across a range of architectures for the C code, install Dan Bernstein's libcpucycles utility from https://cpucycles.cr.yp.to/ . It may be necessary to run *ldconfig*.

# Quick start

For a quick start copy the files from here into a working directory, and try

	python3 pseudo.py 64 2**255-19
	./time

Then 64-bit code for the suggested modulus is generated and tested. An executable that times important functions will be created if the platform allows it. The standalone C timing code is output to *time.c*, and code for production use is output to *field.c*

For Rust 

	python3 pseudo_rust.py 64 2**255-19
	./time

In this case the output is directed to files *time.rs* and *field.rs*

For more details read the comments in the provided scripts.

As a Proof of Concept, elliptic curve code for RFC7748 is provided in the files *rfc7748.c* and *rfc7748.rs*

Read comments in these files for simple build instructions.

The Rust generated code uses wrapping arithmetic, so will panic due to integer overflow in Debug mode. Always run in Release mode.

(RFC7748 describes an implementation of Diffie-Hellman key exchange on the Montgomery elliptic curves X25519 and X448. You can create your own Montgomery and Edwards curve using this utility https://github.com/mcarrickscott/findmontgomery.git )

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
*mod2r()*  -- set to 2^r \
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

## Testing (new)

Some testing of code correctness is carried out by the scripts themselves. However field arithmetic can be further tested against a range of corner/edge cases using the Python script *edge.py*. See the comments at the start of that file for details.


# New - automatic generation of Edwards or Weierstrass curve API code in C

1. Copy all code from this directory to a working directory.
2. Run the python script curve.py, selecting curve name and wordlength (32 or 64 bit). The curve will be in either Edwards or Weierstrass form. Edit this script to add your own curve
3. The script generates 4 files *field.c* *group.c* *point.h* and *curve.c*
4. The script automatically injects this code into *edwards.c* or *weierstrass.c* and *curve.h*, and into EdDSA or ECDSA code if it exists, as for example into *ed448.c* and *nist256.c*

The API interface is as indicated in *curve.h*. The curve API is completely implemented in *edwards.c* or *weierstrass.c*

MAKE SURE to copy fresh copies of all files from source after each test, as some are modified by the python script

Support for standard SHA2 and SHA3 hashing algorithms is provided in *hash.c* and *hash.h*

## Quickstart 1:-

	python3 curve.py 64 ED25519
	gcc -O2 testcurve.c edwards.c -lcpucycles -o testcurve
	./testcurve

Note that this intermediate API only provides the elliptic curve functionality. A higher level algorithm API (like that provided for ed448 EdDSA signature) would use this API while itself providing additional algorithm specific random number and hashing functionality. It may also use the *monty.py* script to generate code to perform arithmetic modulo the prime group order, if so required by the algorithm.

## Quickstart 2 (EdDSA):-

	python3 curve.py 64 ED448
	gcc -O2 ed448.c edwards.c hash.c -o ed448
	./ed448

## Quickstart 3 (ECDSA):-

	python3 curve.py 64 NIST256
	gcc -O2 nist256.c weierstrass.c hash.c -o nist256
	./nist256

# Rust version

Copy all code from this directory to a working directory, and move to that directory. Create a standard rust project, named "ecc".

	cargo new ecc

Replace the default *Cargo.toml* file with the one provided. 

Copy *testcurve.rs* and *hash.rs* into the rust project src subdirectory. Delete the default *main.rs*

In the working directory

	python3 curve_rust.py 64 ED25519
Copy *edwards.rs* into the rust project src subdirectory. 

	cd ecc
	cargo run --release --features ED25519 --bin testcurve

Make sure to copy fresh copies of *edwards.rs* and *weierstrass.rs* from source to the working directory after each test.

	python3 curve_rust.py 64 NIST256

Copy *nist256.rs* and *weierstrass.rs* to the rust src subdirectory
	
	cd ecc
	cargo run --release --bin nist256

Make sure to copy fresh copies of *edwards.rs* and *weierstrass.rs* from source to the working directory after each test.

	python3 curve_rust.py 64 ED448

Copy *ed448.rs* and *edwards.rs* to the rust src subdirectory
	
	cd ecc
	cargo run --release --bin ed448

# A lazy reduction optimization

For the *ed448.c* example, modular additions and subtractions can be speeded up by first editing *monty.py* and *monty_rust.py* and setting *generic=False*.

This activates a curve-specific optimization described in https://eprint.iacr.org/2017/437

# Create your own ECC scheme

Here we describe the steps involved to create your own ECC based signature scheme in C. The process for Rust is very similar.

Copy the *signature.c* template to *<curve_name>.c* for the implementation of an API for a digital signature scheme. For working examples see *ed448.c* and *nist256.c* 

Start with a clean download of all of the files in this directory to a working directory.

Next ensure that the elliptic curve is named and supported in either *pseudo.py* or *monty.py*. 

Make sure the various elliptic curve constants are given in *curve.py* under the curve name. Then for a 64-bit build

	python3 curve.py 64 <Curve Name>

This will automatically create or modify the files *curve.h*, and either *edwards.c* or *weierstrass.c*. 

Now edit the *<curve_name>.c* file. Access the document that outlines the details of the signature scheme, and complete the API implementation. For working examples refer to *ed448.c* (EdDSA) and *nist256.c* (ECDSA)

Provide a main program in *<curve_name>.c* which implements some test vectors and compile as

	gcc -O2 <curve_name>.c <weierstrass.c or edwards.c> hash.c -o <curve_name>

# Testing

Any implementation should be tested using test vectors. A good source are the test vectors provided by the Wycheproof project - see https://github.com/C2SP/wycheproof

A short python script *parse.py* is provided which converts the JSON formatted test vectors for ECDSA and EdDSA to a form more easily digestible by a C or Rust program.


# Microsoft C 64-bit compiler support

Inexcusably the Microsoft C 64 bit compiler MSVC does not support 128-bit integers (BTW you really should be using *clang-cl*). However relatively efficient implementation is still possible with MSVC using special intrinsic functions.

The scripts *pseudoms64.py* and *montyms64.py* generate code which uses these intrinsics. By default x86-64 intrinsics are used. But generated code can also be configured for ARM64, or indeed any architecture which at a minimum allows access to the top 64 bits of a 128-bit product.

A minor change is also required in the *curve.py* script where indicated.  


