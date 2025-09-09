// Program to implement RFC7748 - https://datatracker.ietf.org/doc/html/rfc7748
// Montgomery curve key exchange code, as used by TLS
// Use associated python scripts to generate code for X25519 or X448, but easily modified for other Montgomery curves
//
// A good Montgomery curve can be found by running the sagemath script bowe.sage
//
// Mike Scott 23rd November 2023
// TII
//
// code for 16/32/64-bit processor for X25519 curve can be generated  by 
//
// python pseudo.py 16/32/64 X25519
// or
// python monty.py 16/32/64 X25519
//
// code for 16/32/64-bit processor for X448 curve can be generated  by
//
// python monty.py 16/32/64 X448

// make sure decoration and generic are both set to False
// Seems to prefer clang compiler and karatsuba set to False for X25519 and True for X448
// clang -O3 -march=native -mtune=native rfc7748.c -lcpucycles -o rfc7748

/*** Insert automatically generated code for modulus field.c here ***/



// Automatically generated modular arithmetic C code
// Command line : python monty_neon.py X448
// Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)

#include <arm_neon.h>
#include <stdint.h>
#include <stdio.h>
#define sspint int32x2_t
#define spint uint32x2_t
#define udpint uint64x2_t
#define dpint uint64x2_t

#define Wordlength 32
#define Nlimbs 16
#define Radix 28
#define Nbits 448
#define Nbytes 56

#define MONTGOMERY
#define X448

#define MULBYINT
// propagate carries
static spint inline prop(spint *n) {
  int i;
  spint mask = vdup_n_u32((1 << 28u) - 1);
  sspint carry = vreinterpret_s32_u32(n[0]);
  carry = vshr_n_s32(carry, 28u);
  n[0] = vand_u32(n[0], mask);
  for (i = 1; i < 15; i++) {
    carry = vadd_s32(carry, vreinterpret_s32_u32(n[i]));
    n[i] = vand_u32(vreinterpret_u32_s32(carry), mask);
    carry = vshr_n_s32(carry, 28u);
  }
  n[15] = vadd_u32(n[15], vreinterpret_u32_s32(carry));
  return vreinterpret_u32_s32(vshr_n_s32(vreinterpret_s32_u32(n[15]), 31u));
}

// propagate carries and add p if negative, propagate carries again
static spint flatten(spint *n) {
  int q = (1 << 28u);
  spint one = vdup_n_u32(1);
  spint carry = prop(n);
  spint mpy;
  mpy = vdup_n_u32(1u);
  n[0] = vsub_u32(n[0], vand_u32(mpy, carry));
  mpy = vdup_n_u32(1u);
  n[8] = vsub_u32(n[8], vand_u32(mpy, carry));
  mpy = vdup_n_u32(1u * q);
  n[15] = vadd_u32(n[15], vand_u32(mpy, carry));
  (void)prop(n);
  return vand_u32(carry, one);
}

// Montgomery final subtract
static spint modfsb(spint *n) {
  int q = (1 << 28u);
  spint mpy;
  mpy = vdup_n_u32(1u);
  n[0] = vadd_u32(n[0], mpy);
  mpy = vdup_n_u32(1u);
  n[8] = vadd_u32(n[8], mpy);
  mpy = vdup_n_u32(1u * q);
  n[15] = vsub_u32(n[15], mpy);
  return flatten(n);
}

// Modular addition - reduce less than 2p
static void modadd(const spint *a, const spint *b, spint *n) {
  spint mpy;
  int q = (1 << 28u);
  spint carry;
  n[0] = vadd_u32(a[0], b[0]);
  n[1] = vadd_u32(a[1], b[1]);
  n[2] = vadd_u32(a[2], b[2]);
  n[3] = vadd_u32(a[3], b[3]);
  n[4] = vadd_u32(a[4], b[4]);
  n[5] = vadd_u32(a[5], b[5]);
  n[6] = vadd_u32(a[6], b[6]);
  n[7] = vadd_u32(a[7], b[7]);
  n[8] = vadd_u32(a[8], b[8]);
  n[9] = vadd_u32(a[9], b[9]);
  n[10] = vadd_u32(a[10], b[10]);
  n[11] = vadd_u32(a[11], b[11]);
  n[12] = vadd_u32(a[12], b[12]);
  n[13] = vadd_u32(a[13], b[13]);
  n[14] = vadd_u32(a[14], b[14]);
  n[15] = vadd_u32(a[15], b[15]);
  mpy = vdup_n_u32(2u);
  n[0] = vadd_u32(n[0], mpy);
  mpy = vdup_n_u32(2u);
  n[8] = vadd_u32(n[8], mpy);
  mpy = vdup_n_u32(2u * q);
  n[15] = vsub_u32(n[15], mpy);
  carry = prop(n);
  mpy = vdup_n_u32(2u);
  n[0] = vsub_u32(n[0], vand_u32(mpy, carry));
  mpy = vdup_n_u32(2u);
  n[8] = vsub_u32(n[8], vand_u32(mpy, carry));
  mpy = vdup_n_u32(2u * q);
  n[15] = vadd_u32(n[15], vand_u32(mpy, carry));
  (void)prop(n);
}

// Modular subtraction - reduce less than 2p
static void modsub(const spint *a, const spint *b, spint *n) {
  spint mpy;
  int q = (1 << 28u);
  spint carry;
  n[0] = vsub_u32(a[0], b[0]);
  n[1] = vsub_u32(a[1], b[1]);
  n[2] = vsub_u32(a[2], b[2]);
  n[3] = vsub_u32(a[3], b[3]);
  n[4] = vsub_u32(a[4], b[4]);
  n[5] = vsub_u32(a[5], b[5]);
  n[6] = vsub_u32(a[6], b[6]);
  n[7] = vsub_u32(a[7], b[7]);
  n[8] = vsub_u32(a[8], b[8]);
  n[9] = vsub_u32(a[9], b[9]);
  n[10] = vsub_u32(a[10], b[10]);
  n[11] = vsub_u32(a[11], b[11]);
  n[12] = vsub_u32(a[12], b[12]);
  n[13] = vsub_u32(a[13], b[13]);
  n[14] = vsub_u32(a[14], b[14]);
  n[15] = vsub_u32(a[15], b[15]);
  carry = prop(n);
  mpy = vdup_n_u32(2u);
  n[0] = vsub_u32(n[0], vand_u32(mpy, carry));
  mpy = vdup_n_u32(2u);
  n[8] = vsub_u32(n[8], vand_u32(mpy, carry));
  mpy = vdup_n_u32(2u * q);
  n[15] = vadd_u32(n[15], vand_u32(mpy, carry));
  (void)prop(n);
}

// Modular negation
static void modneg(const spint *b, spint *n) {
  spint mpy;
  spint zero = vdup_n_u32(0);
  int q = ((int)1 << 28u);
  spint carry;
  n[0] = vsub_u32(zero, b[0]);
  n[1] = vsub_u32(zero, b[1]);
  n[2] = vsub_u32(zero, b[2]);
  n[3] = vsub_u32(zero, b[3]);
  n[4] = vsub_u32(zero, b[4]);
  n[5] = vsub_u32(zero, b[5]);
  n[6] = vsub_u32(zero, b[6]);
  n[7] = vsub_u32(zero, b[7]);
  n[8] = vsub_u32(zero, b[8]);
  n[9] = vsub_u32(zero, b[9]);
  n[10] = vsub_u32(zero, b[10]);
  n[11] = vsub_u32(zero, b[11]);
  n[12] = vsub_u32(zero, b[12]);
  n[13] = vsub_u32(zero, b[13]);
  n[14] = vsub_u32(zero, b[14]);
  n[15] = vsub_u32(zero, b[15]);
  carry = prop(n);
  mpy = vdup_n_u32(2u);
  n[0] = vsub_u32(n[0], vand_u32(mpy, carry));
  mpy = vdup_n_u32(2u);
  n[8] = vsub_u32(n[8], vand_u32(mpy, carry));
  mpy = vdup_n_u32(2u * q);
  n[15] = vadd_u32(n[15], vand_u32(mpy, carry));
  (void)prop(n);
}

// Overflow limit   = 18446744073709551616
// maximum possible = 1152921565004824591
// Modular multiplication, c=a*b mod 2p
static void modmul(const spint *a, const spint *b, spint *c) {
  dpint t = vdupq_n_u64(0);
  spint one = vdup_n_u32(1);
  spint q = vdup_n_u32(1 << 28u);
  spint mask = vdup_n_u32((1 << 28u) - 1);
  t = vmlal_u32(t, a[0], b[0]);
  spint v0 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[1]);
  t = vmlal_u32(t, a[1], b[0]);
  spint v1 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[2]);
  t = vmlal_u32(t, a[1], b[1]);
  t = vmlal_u32(t, a[2], b[0]);
  spint v2 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[3]);
  t = vmlal_u32(t, a[1], b[2]);
  t = vmlal_u32(t, a[2], b[1]);
  t = vmlal_u32(t, a[3], b[0]);
  spint v3 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[4]);
  t = vmlal_u32(t, a[1], b[3]);
  t = vmlal_u32(t, a[2], b[2]);
  t = vmlal_u32(t, a[3], b[1]);
  t = vmlal_u32(t, a[4], b[0]);
  spint v4 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[5]);
  t = vmlal_u32(t, a[1], b[4]);
  t = vmlal_u32(t, a[2], b[3]);
  t = vmlal_u32(t, a[3], b[2]);
  t = vmlal_u32(t, a[4], b[1]);
  t = vmlal_u32(t, a[5], b[0]);
  spint v5 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[6]);
  t = vmlal_u32(t, a[1], b[5]);
  t = vmlal_u32(t, a[2], b[4]);
  t = vmlal_u32(t, a[3], b[3]);
  t = vmlal_u32(t, a[4], b[2]);
  t = vmlal_u32(t, a[5], b[1]);
  t = vmlal_u32(t, a[6], b[0]);
  spint v6 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[7]);
  t = vmlal_u32(t, a[1], b[6]);
  t = vmlal_u32(t, a[2], b[5]);
  t = vmlal_u32(t, a[3], b[4]);
  t = vmlal_u32(t, a[4], b[3]);
  t = vmlal_u32(t, a[5], b[2]);
  t = vmlal_u32(t, a[6], b[1]);
  t = vmlal_u32(t, a[7], b[0]);
  spint v7 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[8]);
  t = vmlal_u32(t, a[1], b[7]);
  t = vmlal_u32(t, a[2], b[6]);
  t = vmlal_u32(t, a[3], b[5]);
  t = vmlal_u32(t, a[4], b[4]);
  t = vmlal_u32(t, a[5], b[3]);
  t = vmlal_u32(t, a[6], b[2]);
  t = vmlal_u32(t, a[7], b[1]);
  t = vmlal_u32(t, a[8], b[0]);
  t = vaddw_u32(t, vsub_u32(q, v0));
  spint v8 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[9]);
  t = vmlal_u32(t, a[1], b[8]);
  t = vmlal_u32(t, a[2], b[7]);
  t = vmlal_u32(t, a[3], b[6]);
  t = vmlal_u32(t, a[4], b[5]);
  t = vmlal_u32(t, a[5], b[4]);
  t = vmlal_u32(t, a[6], b[3]);
  t = vmlal_u32(t, a[7], b[2]);
  t = vmlal_u32(t, a[8], b[1]);
  t = vmlal_u32(t, a[9], b[0]);
  spint s = (spint)mask;
  s = vsub_u32(s, v1);
  t = vaddw_u32(t, s);
  spint v9 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[10]);
  t = vmlal_u32(t, a[1], b[9]);
  t = vmlal_u32(t, a[2], b[8]);
  t = vmlal_u32(t, a[3], b[7]);
  t = vmlal_u32(t, a[4], b[6]);
  t = vmlal_u32(t, a[5], b[5]);
  t = vmlal_u32(t, a[6], b[4]);
  t = vmlal_u32(t, a[7], b[3]);
  t = vmlal_u32(t, a[8], b[2]);
  t = vmlal_u32(t, a[9], b[1]);
  t = vmlal_u32(t, a[10], b[0]);
  s = (spint)mask;
  s = vsub_u32(s, v2);
  t = vaddw_u32(t, s);
  spint v10 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[11]);
  t = vmlal_u32(t, a[1], b[10]);
  t = vmlal_u32(t, a[2], b[9]);
  t = vmlal_u32(t, a[3], b[8]);
  t = vmlal_u32(t, a[4], b[7]);
  t = vmlal_u32(t, a[5], b[6]);
  t = vmlal_u32(t, a[6], b[5]);
  t = vmlal_u32(t, a[7], b[4]);
  t = vmlal_u32(t, a[8], b[3]);
  t = vmlal_u32(t, a[9], b[2]);
  t = vmlal_u32(t, a[10], b[1]);
  t = vmlal_u32(t, a[11], b[0]);
  s = (spint)mask;
  s = vsub_u32(s, v3);
  t = vaddw_u32(t, s);
  spint v11 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[12]);
  t = vmlal_u32(t, a[1], b[11]);
  t = vmlal_u32(t, a[2], b[10]);
  t = vmlal_u32(t, a[3], b[9]);
  t = vmlal_u32(t, a[4], b[8]);
  t = vmlal_u32(t, a[5], b[7]);
  t = vmlal_u32(t, a[6], b[6]);
  t = vmlal_u32(t, a[7], b[5]);
  t = vmlal_u32(t, a[8], b[4]);
  t = vmlal_u32(t, a[9], b[3]);
  t = vmlal_u32(t, a[10], b[2]);
  t = vmlal_u32(t, a[11], b[1]);
  t = vmlal_u32(t, a[12], b[0]);
  s = (spint)mask;
  s = vsub_u32(s, v4);
  t = vaddw_u32(t, s);
  spint v12 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[13]);
  t = vmlal_u32(t, a[1], b[12]);
  t = vmlal_u32(t, a[2], b[11]);
  t = vmlal_u32(t, a[3], b[10]);
  t = vmlal_u32(t, a[4], b[9]);
  t = vmlal_u32(t, a[5], b[8]);
  t = vmlal_u32(t, a[6], b[7]);
  t = vmlal_u32(t, a[7], b[6]);
  t = vmlal_u32(t, a[8], b[5]);
  t = vmlal_u32(t, a[9], b[4]);
  t = vmlal_u32(t, a[10], b[3]);
  t = vmlal_u32(t, a[11], b[2]);
  t = vmlal_u32(t, a[12], b[1]);
  t = vmlal_u32(t, a[13], b[0]);
  s = (spint)mask;
  s = vsub_u32(s, v5);
  t = vaddw_u32(t, s);
  spint v13 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[14]);
  t = vmlal_u32(t, a[1], b[13]);
  t = vmlal_u32(t, a[2], b[12]);
  t = vmlal_u32(t, a[3], b[11]);
  t = vmlal_u32(t, a[4], b[10]);
  t = vmlal_u32(t, a[5], b[9]);
  t = vmlal_u32(t, a[6], b[8]);
  t = vmlal_u32(t, a[7], b[7]);
  t = vmlal_u32(t, a[8], b[6]);
  t = vmlal_u32(t, a[9], b[5]);
  t = vmlal_u32(t, a[10], b[4]);
  t = vmlal_u32(t, a[11], b[3]);
  t = vmlal_u32(t, a[12], b[2]);
  t = vmlal_u32(t, a[13], b[1]);
  t = vmlal_u32(t, a[14], b[0]);
  s = (spint)mask;
  s = vsub_u32(s, v6);
  t = vaddw_u32(t, s);
  spint v14 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[0], b[15]);
  t = vmlal_u32(t, a[1], b[14]);
  t = vmlal_u32(t, a[2], b[13]);
  t = vmlal_u32(t, a[3], b[12]);
  t = vmlal_u32(t, a[4], b[11]);
  t = vmlal_u32(t, a[5], b[10]);
  t = vmlal_u32(t, a[6], b[9]);
  t = vmlal_u32(t, a[7], b[8]);
  t = vmlal_u32(t, a[8], b[7]);
  t = vmlal_u32(t, a[9], b[6]);
  t = vmlal_u32(t, a[10], b[5]);
  t = vmlal_u32(t, a[11], b[4]);
  t = vmlal_u32(t, a[12], b[3]);
  t = vmlal_u32(t, a[13], b[2]);
  t = vmlal_u32(t, a[14], b[1]);
  t = vmlal_u32(t, a[15], b[0]);
  s = (spint)mask;
  s = vsub_u32(s, v7);
  t = vaddw_u32(t, s);
  spint v15 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[1], b[15]);
  t = vmlal_u32(t, a[2], b[14]);
  t = vmlal_u32(t, a[3], b[13]);
  t = vmlal_u32(t, a[4], b[12]);
  t = vmlal_u32(t, a[5], b[11]);
  t = vmlal_u32(t, a[6], b[10]);
  t = vmlal_u32(t, a[7], b[9]);
  t = vmlal_u32(t, a[8], b[8]);
  t = vmlal_u32(t, a[9], b[7]);
  t = vmlal_u32(t, a[10], b[6]);
  t = vmlal_u32(t, a[11], b[5]);
  t = vmlal_u32(t, a[12], b[4]);
  t = vmlal_u32(t, a[13], b[3]);
  t = vmlal_u32(t, a[14], b[2]);
  t = vmlal_u32(t, a[15], b[1]);
  s = (spint)mask;
  s = vadd_u32(s, v0);
  s = vsub_u32(s, v8);
  t = vaddw_u32(t, s);
  spint v16 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[2], b[15]);
  t = vmlal_u32(t, a[3], b[14]);
  t = vmlal_u32(t, a[4], b[13]);
  t = vmlal_u32(t, a[5], b[12]);
  t = vmlal_u32(t, a[6], b[11]);
  t = vmlal_u32(t, a[7], b[10]);
  t = vmlal_u32(t, a[8], b[9]);
  t = vmlal_u32(t, a[9], b[8]);
  t = vmlal_u32(t, a[10], b[7]);
  t = vmlal_u32(t, a[11], b[6]);
  t = vmlal_u32(t, a[12], b[5]);
  t = vmlal_u32(t, a[13], b[4]);
  t = vmlal_u32(t, a[14], b[3]);
  t = vmlal_u32(t, a[15], b[2]);
  s = (spint)mask;
  s = vadd_u32(s, v1);
  s = vsub_u32(s, v9);
  t = vaddw_u32(t, s);
  c[0] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[3], b[15]);
  t = vmlal_u32(t, a[4], b[14]);
  t = vmlal_u32(t, a[5], b[13]);
  t = vmlal_u32(t, a[6], b[12]);
  t = vmlal_u32(t, a[7], b[11]);
  t = vmlal_u32(t, a[8], b[10]);
  t = vmlal_u32(t, a[9], b[9]);
  t = vmlal_u32(t, a[10], b[8]);
  t = vmlal_u32(t, a[11], b[7]);
  t = vmlal_u32(t, a[12], b[6]);
  t = vmlal_u32(t, a[13], b[5]);
  t = vmlal_u32(t, a[14], b[4]);
  t = vmlal_u32(t, a[15], b[3]);
  s = (spint)mask;
  s = vadd_u32(s, v2);
  s = vsub_u32(s, v10);
  t = vaddw_u32(t, s);
  c[1] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[4], b[15]);
  t = vmlal_u32(t, a[5], b[14]);
  t = vmlal_u32(t, a[6], b[13]);
  t = vmlal_u32(t, a[7], b[12]);
  t = vmlal_u32(t, a[8], b[11]);
  t = vmlal_u32(t, a[9], b[10]);
  t = vmlal_u32(t, a[10], b[9]);
  t = vmlal_u32(t, a[11], b[8]);
  t = vmlal_u32(t, a[12], b[7]);
  t = vmlal_u32(t, a[13], b[6]);
  t = vmlal_u32(t, a[14], b[5]);
  t = vmlal_u32(t, a[15], b[4]);
  s = (spint)mask;
  s = vadd_u32(s, v3);
  s = vsub_u32(s, v11);
  t = vaddw_u32(t, s);
  c[2] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[5], b[15]);
  t = vmlal_u32(t, a[6], b[14]);
  t = vmlal_u32(t, a[7], b[13]);
  t = vmlal_u32(t, a[8], b[12]);
  t = vmlal_u32(t, a[9], b[11]);
  t = vmlal_u32(t, a[10], b[10]);
  t = vmlal_u32(t, a[11], b[9]);
  t = vmlal_u32(t, a[12], b[8]);
  t = vmlal_u32(t, a[13], b[7]);
  t = vmlal_u32(t, a[14], b[6]);
  t = vmlal_u32(t, a[15], b[5]);
  s = (spint)mask;
  s = vadd_u32(s, v4);
  s = vsub_u32(s, v12);
  t = vaddw_u32(t, s);
  c[3] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[6], b[15]);
  t = vmlal_u32(t, a[7], b[14]);
  t = vmlal_u32(t, a[8], b[13]);
  t = vmlal_u32(t, a[9], b[12]);
  t = vmlal_u32(t, a[10], b[11]);
  t = vmlal_u32(t, a[11], b[10]);
  t = vmlal_u32(t, a[12], b[9]);
  t = vmlal_u32(t, a[13], b[8]);
  t = vmlal_u32(t, a[14], b[7]);
  t = vmlal_u32(t, a[15], b[6]);
  s = (spint)mask;
  s = vadd_u32(s, v5);
  s = vsub_u32(s, v13);
  t = vaddw_u32(t, s);
  c[4] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[7], b[15]);
  t = vmlal_u32(t, a[8], b[14]);
  t = vmlal_u32(t, a[9], b[13]);
  t = vmlal_u32(t, a[10], b[12]);
  t = vmlal_u32(t, a[11], b[11]);
  t = vmlal_u32(t, a[12], b[10]);
  t = vmlal_u32(t, a[13], b[9]);
  t = vmlal_u32(t, a[14], b[8]);
  t = vmlal_u32(t, a[15], b[7]);
  s = (spint)mask;
  s = vadd_u32(s, v6);
  s = vsub_u32(s, v14);
  t = vaddw_u32(t, s);
  c[5] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[8], b[15]);
  t = vmlal_u32(t, a[9], b[14]);
  t = vmlal_u32(t, a[10], b[13]);
  t = vmlal_u32(t, a[11], b[12]);
  t = vmlal_u32(t, a[12], b[11]);
  t = vmlal_u32(t, a[13], b[10]);
  t = vmlal_u32(t, a[14], b[9]);
  t = vmlal_u32(t, a[15], b[8]);
  s = (spint)mask;
  s = vadd_u32(s, v7);
  s = vsub_u32(s, v15);
  t = vaddw_u32(t, s);
  c[6] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[9], b[15]);
  t = vmlal_u32(t, a[10], b[14]);
  t = vmlal_u32(t, a[11], b[13]);
  t = vmlal_u32(t, a[12], b[12]);
  t = vmlal_u32(t, a[13], b[11]);
  t = vmlal_u32(t, a[14], b[10]);
  t = vmlal_u32(t, a[15], b[9]);
  s = (spint)mask;
  s = vadd_u32(s, v8);
  s = vsub_u32(s, v16);
  t = vaddw_u32(t, s);
  c[7] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[10], b[15]);
  t = vmlal_u32(t, a[11], b[14]);
  t = vmlal_u32(t, a[12], b[13]);
  t = vmlal_u32(t, a[13], b[12]);
  t = vmlal_u32(t, a[14], b[11]);
  t = vmlal_u32(t, a[15], b[10]);
  s = (spint)mask;
  s = vadd_u32(s, v9);
  t = vaddw_u32(t, s);
  c[8] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[11], b[15]);
  t = vmlal_u32(t, a[12], b[14]);
  t = vmlal_u32(t, a[13], b[13]);
  t = vmlal_u32(t, a[14], b[12]);
  t = vmlal_u32(t, a[15], b[11]);
  s = (spint)mask;
  s = vadd_u32(s, v10);
  t = vaddw_u32(t, s);
  c[9] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[12], b[15]);
  t = vmlal_u32(t, a[13], b[14]);
  t = vmlal_u32(t, a[14], b[13]);
  t = vmlal_u32(t, a[15], b[12]);
  s = (spint)mask;
  s = vadd_u32(s, v11);
  t = vaddw_u32(t, s);
  c[10] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[13], b[15]);
  t = vmlal_u32(t, a[14], b[14]);
  t = vmlal_u32(t, a[15], b[13]);
  s = (spint)mask;
  s = vadd_u32(s, v12);
  t = vaddw_u32(t, s);
  c[11] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[14], b[15]);
  t = vmlal_u32(t, a[15], b[14]);
  s = (spint)mask;
  s = vadd_u32(s, v13);
  t = vaddw_u32(t, s);
  c[12] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vmlal_u32(t, a[15], b[15]);
  s = (spint)mask;
  s = vadd_u32(s, v14);
  t = vaddw_u32(t, s);
  c[13] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  s = (spint)mask;
  s = vadd_u32(s, v15);
  t = vaddw_u32(t, s);
  c[14] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vaddw_u32(t, vsub_u32(v16, one));
  c[15] = vmovn_u64(t);
}

// Modular squaring, c=a*a  mod 2p
static void modsqr(const spint *a, spint *c) {
  udpint tot;
  udpint t = vdupq_n_u64(0);
  spint one = vdup_n_u32(1);
  spint q = vdup_n_u32(1 << 28u);
  spint mask = vdup_n_u32((1 << 28u) - 1);
  tot = vmull_u32(a[0], a[0]);
  t = tot;
  spint v0 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[1]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  spint v1 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[2]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[1], a[1]);
  t = vaddq_u64(t, tot);
  spint v2 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[3]);
  tot = vmlal_u32(tot, a[1], a[2]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  spint v3 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[4]);
  tot = vmlal_u32(tot, a[1], a[3]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[2], a[2]);
  t = vaddq_u64(t, tot);
  spint v4 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[5]);
  tot = vmlal_u32(tot, a[1], a[4]);
  tot = vmlal_u32(tot, a[2], a[3]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  spint v5 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[6]);
  tot = vmlal_u32(tot, a[1], a[5]);
  tot = vmlal_u32(tot, a[2], a[4]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[3], a[3]);
  t = vaddq_u64(t, tot);
  spint v6 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[7]);
  tot = vmlal_u32(tot, a[1], a[6]);
  tot = vmlal_u32(tot, a[2], a[5]);
  tot = vmlal_u32(tot, a[3], a[4]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  spint v7 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[8]);
  tot = vmlal_u32(tot, a[1], a[7]);
  tot = vmlal_u32(tot, a[2], a[6]);
  tot = vmlal_u32(tot, a[3], a[5]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[4], a[4]);
  t = vaddq_u64(t, tot);
  t = vaddw_u32(t, vsub_u32(q, v0));
  spint v8 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[9]);
  tot = vmlal_u32(tot, a[1], a[8]);
  tot = vmlal_u32(tot, a[2], a[7]);
  tot = vmlal_u32(tot, a[3], a[6]);
  tot = vmlal_u32(tot, a[4], a[5]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  spint s = (spint)mask;
  s = vsub_u32(s, v1);
  t = vaddw_u32(t, s);
  spint v9 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[10]);
  tot = vmlal_u32(tot, a[1], a[9]);
  tot = vmlal_u32(tot, a[2], a[8]);
  tot = vmlal_u32(tot, a[3], a[7]);
  tot = vmlal_u32(tot, a[4], a[6]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[5], a[5]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vsub_u32(s, v2);
  t = vaddw_u32(t, s);
  spint v10 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[11]);
  tot = vmlal_u32(tot, a[1], a[10]);
  tot = vmlal_u32(tot, a[2], a[9]);
  tot = vmlal_u32(tot, a[3], a[8]);
  tot = vmlal_u32(tot, a[4], a[7]);
  tot = vmlal_u32(tot, a[5], a[6]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vsub_u32(s, v3);
  t = vaddw_u32(t, s);
  spint v11 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[12]);
  tot = vmlal_u32(tot, a[1], a[11]);
  tot = vmlal_u32(tot, a[2], a[10]);
  tot = vmlal_u32(tot, a[3], a[9]);
  tot = vmlal_u32(tot, a[4], a[8]);
  tot = vmlal_u32(tot, a[5], a[7]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[6], a[6]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vsub_u32(s, v4);
  t = vaddw_u32(t, s);
  spint v12 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[13]);
  tot = vmlal_u32(tot, a[1], a[12]);
  tot = vmlal_u32(tot, a[2], a[11]);
  tot = vmlal_u32(tot, a[3], a[10]);
  tot = vmlal_u32(tot, a[4], a[9]);
  tot = vmlal_u32(tot, a[5], a[8]);
  tot = vmlal_u32(tot, a[6], a[7]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vsub_u32(s, v5);
  t = vaddw_u32(t, s);
  spint v13 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[14]);
  tot = vmlal_u32(tot, a[1], a[13]);
  tot = vmlal_u32(tot, a[2], a[12]);
  tot = vmlal_u32(tot, a[3], a[11]);
  tot = vmlal_u32(tot, a[4], a[10]);
  tot = vmlal_u32(tot, a[5], a[9]);
  tot = vmlal_u32(tot, a[6], a[8]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[7], a[7]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vsub_u32(s, v6);
  t = vaddw_u32(t, s);
  spint v14 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[0], a[15]);
  tot = vmlal_u32(tot, a[1], a[14]);
  tot = vmlal_u32(tot, a[2], a[13]);
  tot = vmlal_u32(tot, a[3], a[12]);
  tot = vmlal_u32(tot, a[4], a[11]);
  tot = vmlal_u32(tot, a[5], a[10]);
  tot = vmlal_u32(tot, a[6], a[9]);
  tot = vmlal_u32(tot, a[7], a[8]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vsub_u32(s, v7);
  t = vaddw_u32(t, s);
  spint v15 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[1], a[15]);
  tot = vmlal_u32(tot, a[2], a[14]);
  tot = vmlal_u32(tot, a[3], a[13]);
  tot = vmlal_u32(tot, a[4], a[12]);
  tot = vmlal_u32(tot, a[5], a[11]);
  tot = vmlal_u32(tot, a[6], a[10]);
  tot = vmlal_u32(tot, a[7], a[9]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[8], a[8]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v0);
  s = vsub_u32(s, v8);
  t = vaddw_u32(t, s);
  spint v16 = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[2], a[15]);
  tot = vmlal_u32(tot, a[3], a[14]);
  tot = vmlal_u32(tot, a[4], a[13]);
  tot = vmlal_u32(tot, a[5], a[12]);
  tot = vmlal_u32(tot, a[6], a[11]);
  tot = vmlal_u32(tot, a[7], a[10]);
  tot = vmlal_u32(tot, a[8], a[9]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v1);
  s = vsub_u32(s, v9);
  t = vaddw_u32(t, s);
  c[0] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[3], a[15]);
  tot = vmlal_u32(tot, a[4], a[14]);
  tot = vmlal_u32(tot, a[5], a[13]);
  tot = vmlal_u32(tot, a[6], a[12]);
  tot = vmlal_u32(tot, a[7], a[11]);
  tot = vmlal_u32(tot, a[8], a[10]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[9], a[9]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v2);
  s = vsub_u32(s, v10);
  t = vaddw_u32(t, s);
  c[1] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[4], a[15]);
  tot = vmlal_u32(tot, a[5], a[14]);
  tot = vmlal_u32(tot, a[6], a[13]);
  tot = vmlal_u32(tot, a[7], a[12]);
  tot = vmlal_u32(tot, a[8], a[11]);
  tot = vmlal_u32(tot, a[9], a[10]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v3);
  s = vsub_u32(s, v11);
  t = vaddw_u32(t, s);
  c[2] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[5], a[15]);
  tot = vmlal_u32(tot, a[6], a[14]);
  tot = vmlal_u32(tot, a[7], a[13]);
  tot = vmlal_u32(tot, a[8], a[12]);
  tot = vmlal_u32(tot, a[9], a[11]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[10], a[10]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v4);
  s = vsub_u32(s, v12);
  t = vaddw_u32(t, s);
  c[3] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[6], a[15]);
  tot = vmlal_u32(tot, a[7], a[14]);
  tot = vmlal_u32(tot, a[8], a[13]);
  tot = vmlal_u32(tot, a[9], a[12]);
  tot = vmlal_u32(tot, a[10], a[11]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v5);
  s = vsub_u32(s, v13);
  t = vaddw_u32(t, s);
  c[4] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[7], a[15]);
  tot = vmlal_u32(tot, a[8], a[14]);
  tot = vmlal_u32(tot, a[9], a[13]);
  tot = vmlal_u32(tot, a[10], a[12]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[11], a[11]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v6);
  s = vsub_u32(s, v14);
  t = vaddw_u32(t, s);
  c[5] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[8], a[15]);
  tot = vmlal_u32(tot, a[9], a[14]);
  tot = vmlal_u32(tot, a[10], a[13]);
  tot = vmlal_u32(tot, a[11], a[12]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v7);
  s = vsub_u32(s, v15);
  t = vaddw_u32(t, s);
  c[6] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[9], a[15]);
  tot = vmlal_u32(tot, a[10], a[14]);
  tot = vmlal_u32(tot, a[11], a[13]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[12], a[12]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v8);
  s = vsub_u32(s, v16);
  t = vaddw_u32(t, s);
  c[7] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[10], a[15]);
  tot = vmlal_u32(tot, a[11], a[14]);
  tot = vmlal_u32(tot, a[12], a[13]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v9);
  t = vaddw_u32(t, s);
  c[8] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[11], a[15]);
  tot = vmlal_u32(tot, a[12], a[14]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[13], a[13]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v10);
  t = vaddw_u32(t, s);
  c[9] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[12], a[15]);
  tot = vmlal_u32(tot, a[13], a[14]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v11);
  t = vaddw_u32(t, s);
  c[10] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[13], a[15]);
  tot = vaddq_u64(tot, tot);
  tot = vmlal_u32(tot, a[14], a[14]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v12);
  t = vaddw_u32(t, s);
  c[11] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[14], a[15]);
  tot = vaddq_u64(tot, tot);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v13);
  t = vaddw_u32(t, s);
  c[12] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  tot = vmull_u32(a[15], a[15]);
  t = vaddq_u64(t, tot);
  s = (spint)mask;
  s = vadd_u32(s, v14);
  t = vaddw_u32(t, s);
  c[13] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  s = (spint)mask;
  s = vadd_u32(s, v15);
  t = vaddw_u32(t, s);
  c[14] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28);
  t = vaddw_u32(t, vsub_u32(v16, one));
  c[15] = vmovn_u64(t);
}

// copy
static void modcpy(const spint *a, spint *c) {
  int i;
  for (i = 0; i < 16; i++) {
    c[i] = a[i];
  }
}

// square n times
static void modnsqr(spint *a, int n) {
  int i;
  for (i = 0; i < n; i++) {
    modsqr(a, a);
  }
}

// Calculate progenitor
static void modpro(const spint *w, spint *z) {
  spint x[16];
  spint t0[16];
  spint t1[16];
  modcpy(w, x);
  modsqr(x, z);
  modmul(x, z, z);
  modsqr(z, z);
  modmul(x, z, z);
  modcpy(z, t0);
  modnsqr(t0, 3);
  modmul(z, t0, z);
  modcpy(z, t0);
  modnsqr(t0, 6);
  modmul(z, t0, t0);
  modcpy(t0, t1);
  modnsqr(t1, 12);
  modmul(t0, t1, t0);
  modcpy(t0, t1);
  modnsqr(t1, 6);
  modmul(z, t1, z);
  modnsqr(t1, 18);
  modmul(t0, t1, t0);
  modcpy(t0, t1);
  modnsqr(t1, 48);
  modmul(t0, t1, t0);
  modcpy(t0, t1);
  modnsqr(t1, 96);
  modmul(t0, t1, t0);
  modnsqr(t0, 30);
  modmul(z, t0, z);
  modsqr(z, t0);
  modmul(x, t0, t0);
  modnsqr(t0, 223);
  modmul(z, t0, z);
}

// calculate inverse, provide progenitor h if available
static void modinv(const spint *x, const spint *h, spint *z) {
  spint s[16];
  spint t[16];
  if (h == NULL) {
    modpro(x, t);
  } else {
    modcpy(h, t);
  }
  modcpy(x, s);
  modnsqr(t, 2);
  modmul(s, t, z);
}

// Convert m to n-residue form, n=nres(m)
static void nres(const spint *m, spint *n) {
  spint c[16];
  c[0] = vdup_n_u32(0x0u);
  c[1] = vdup_n_u32(0x0u);
  c[2] = vdup_n_u32(0x2u);
  c[3] = vdup_n_u32(0x0u);
  c[4] = vdup_n_u32(0x0u);
  c[5] = vdup_n_u32(0x0u);
  c[6] = vdup_n_u32(0x0u);
  c[7] = vdup_n_u32(0x0u);
  c[8] = vdup_n_u32(0x0u);
  c[9] = vdup_n_u32(0x0u);
  c[10] = vdup_n_u32(0x3u);
  c[11] = vdup_n_u32(0x0u);
  c[12] = vdup_n_u32(0x0u);
  c[13] = vdup_n_u32(0x0u);
  c[14] = vdup_n_u32(0x0u);
  c[15] = vdup_n_u32(0x0u);
  modmul(m, c, n);
}

// Convert n back to normal form, m=redc(n)
static void redc(const spint *n, spint *m) {
  int i;
  spint c[16];
  c[0] = vdup_n_u32(1);
  for (i = 1; i < 16; i++) {
    c[i] = vdup_n_u32(0);
  }
  modmul(n, c, m);
  (void)modfsb(m);
}

// is unity?
static spint modis1(const spint *a) {
  int i;
  spint c[16];
  spint c0;
  spint d = vdup_n_u32(0);
  spint one = vdup_n_u32(1);
  redc(a, c);
  for (i = 1; i < 16; i++) {
    d = vorr_u32(d, c[i]);
  }
  return vand_u32(vand_u32(one, vshr_n_u32(vsub_u32(d, one), 28u)),
                  vshr_n_u32(vsub_u32(veor_u32(c[0], one), one), 28u));
}

// is zero?
static spint modis0(const spint *a) {
  int i;
  spint c[16];
  spint d = vdup_n_u32(0);
  spint one = vdup_n_u32(1);
  redc(a, c);
  for (i = 0; i < 16; i++) {
    d = vorr_u32(d, c[i]);
  }
  return vand_u32(vshr_n_u32(vsub_u32(d, one), 28u), one);
}

// set to zero
static void modzer(spint *a) {
  int i;
  for (i = 0; i < 16; i++) {
    a[i] = vdup_n_u32(0);
  }
}

// set to one
static void modone(spint *a) {
  int i;
  a[0] = vdup_n_u32(1);
  for (i = 1; i < 16; i++) {
    a[i] = vdup_n_u32(0);
  }
  nres(a, a);
}

// set to integer
static void modint(int x, spint *a) {
  int i;
  a[0] = vdup_n_u32(x);
  for (i = 1; i < 16; i++) {
    a[i] = vdup_n_u32(0);
  }
  nres(a, a);
}

// Modular multiplication by an integer, c=a*b mod 2p
// uses special method for trinomials, otherwise Barrett-Dhem reduction
static void modmli(const spint *a, int b, spint *c) {
  spint mask = vdup_n_u32((1 << 28u) - 1);
  udpint t = vdupq_n_u64(0);
  spint s;
  t = vmlal_n_u32(t, a[0], b);
  c[0] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[1], b);
  c[1] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[2], b);
  c[2] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[3], b);
  c[3] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[4], b);
  c[4] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[5], b);
  c[5] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[6], b);
  c[6] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[7], b);
  c[7] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[8], b);
  c[8] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[9], b);
  c[9] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[10], b);
  c[10] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[11], b);
  c[11] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[12], b);
  c[12] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[13], b);
  c[13] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[14], b);
  c[14] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  t = vmlal_n_u32(t, a[15], b);
  c[15] = vand_u32(vmovn_u64(t), mask);
  t = vshrq_n_u64(t, 28u);
  // reduction pass

  s = vmovn_u64(t);
  c[0] = vadd_u32(c[0], s);
  c[8] = vadd_u32(c[8], s);
}

// Test for quadratic residue
static spint modqr(const spint *h, const spint *x) {
  spint r[16];
  if (h == NULL) {
    modpro(x, r);
    modsqr(r, r);
  } else {
    modsqr(h, r);
  }
  modmul(r, x, r);
  return modis1(r) | modis0(x);
}

// conditional move g to f if d=1
// strongly recommend inlining be disabled using compiler specific syntax
static void __attribute__((noinline))
modcmv(spint b, const spint *g, volatile spint *f) {
  int i;
  spint c0, c1, s, t, w, aux;
  static uint32_t R[2] = {0, 0};
  spint one = vdup_n_u32(1);
  R[0] += 0x5aa5a55au;
  R[1] += 0x7447e88eu;
  w = vld1_u32(R);
  c0 = vand_u32(vmvn_u32(b), vadd_u32(w, one));
  c1 = vadd_u32(b, w);
  for (i = 0; i < 16; i++) {
    s = g[i];
    t = f[i];
    f[i] = aux = vadd_u32(vmul_u32(c0, t), vmul_u32(c1, s));
    f[i] = vsub_u32(aux, vmul_u32(w, vadd_u32(t, s)));
  }
}

// conditional swap g and f if d=1
// strongly recommend inlining be disabled using compiler specific syntax
static void __attribute__((noinline))
modcsw(spint b, volatile spint *g, volatile spint *f) {
  int i;
  spint c0, c1, s, t, w, v, aux;
  static uint32_t R[2] = {0, 0};
  spint one = vdup_n_u32(1);
  R[0] += 0x5aa5a55au;
  R[1] += 0x7447e88eu;
  w = vld1_u32(R);
  c0 = vand_u32(vmvn_u32(b), vadd_u32(w, one));
  c1 = vadd_u32(b, w);
  for (i = 0; i < 16; i++) {
    s = g[i];
    t = f[i];
    v = vmul_u32(w, vadd_u32(t, s));
    f[i] = aux = vadd_u32(vmul_u32(c0, t), vmul_u32(c1, s));
    f[i] = vsub_u32(aux, v);
    g[i] = aux = vadd_u32(vmul_u32(c0, s), vmul_u32(c1, t));
    g[i] = vsub_u32(aux, v);
  }
}

// Modular square root, provide progenitor h if available, NULL if not
static void modsqrt(const spint *x, const spint *h, spint *r) {
  spint s[16];
  spint y[16];
  if (h == NULL) {
    modpro(x, y);
  } else {
    modcpy(h, y);
  }
  modmul(y, x, s);
  modcpy(s, r);
}

// shift left by less than a word
static void modshl(unsigned int n, spint *a) {
  int i;
  spint mask = vdup_n_u32(0xfffffff);
  a[15] = vorr_u32(vshl_n_u32(a[15], n), vshr_n_u32(a[14], 28u - n));
  for (i = 14; i > 0; i--) {
    a[i] = vorr_u32(vand_u32(vshl_n_u32(a[i], n), mask),
                    vshr_n_u32(a[i - 1], 28u - n));
  }
  a[0] = vand_u32(vshl_n_u32(a[0], n), mask);
}

// shift right by less than a word. Return shifted out part
static spint modshr(unsigned int n, spint *a) {
  int i;
  spint mask = vdup_n_u32(0xfffffff);
  spint mskn = vdup_n_u32((1 << n) - 1);
  spint r = vand_u32(a[0], mskn);
  for (i = 0; i < 15; i++) {
    a[i] = vorr_u32(vshr_n_u32(a[i], n),
                    vand_u32(vshl_n_u32(a[i + 1], 28u - n), mask));
  }
  a[15] = vshr_n_u32(a[15], n);
  return r;
}

// set a= 2^r
static void mod2r(unsigned int r, spint *a) {
  unsigned int n = r / 28u;
  unsigned int m = r % 28u;
  spint one = vdup_n_u32(1);
  modzer(a);
  if (r >= 56 * 8)
    return;
  a[n] = one;
  a[n] = vshl_n_u32(a[n], m);
  nres(a, a);
}

// export to byte arrays
static void modexp(const spint *a, char *b, char *e) {
  int i;
  uint32_t s[2];
  spint mask = vdup_n_u32(0xff);
  spint c[16];
  redc(a, c);
  for (i = 55; i >= 0; i--) {
    vst1_u32(s, vand_u32(c[0], mask));
    b[i] = (char)s[0];
    if (e != NULL)
      e[i] = (char)s[1];
    (void)modshr(8, c);
  }
}

// import from byte arrays
// returns 1 if in range, else 0
static spint modimp(const char *b, const char *e, spint *a) {
  int i;
  spint res;
  uint32_t s[2];
  for (i = 0; i < 16; i++) {
    a[i] = vdup_n_u32(0);
  }
  for (i = 0; i < 56; i++) {
    modshl(8, a);
    s[0] = (uint32_t)b[i];
    if (e != NULL)
      s[1] = (uint32_t)e[i];
    else
      s[1] = 0;
    a[0] = vadd_u32(a[0], vld1_u32(s));
  }
  res = modfsb(a);
  nres(a, a);
  return res;
}

// determine sign
static spint modsign(const spint *a) {
  spint c[16];
  spint one = vdup_n_u32(1);
  redc(a, c);
  return vand_u32(c[0], one);
}

// return true if equal
static spint modcmp(const spint *a, const spint *b) {
  spint c[16], d[16];
  int i;
  spint one = vdup_n_u32(1);
  spint eq = one;
  redc(a, c);
  redc(b, d);
  for (i = 0; i < 16; i++) {
    eq = vand_u32(
        eq, vand_u32(vshr_n_u32(vsub_u32(veor_u32(c[i], d[i]), one), 28), one));
  }
  return eq;
}



/*** End of automatically generated code ***/

#define COUNT_CLOCKS
//#define USE_RDTSC

#ifdef COUNT_CLOCKS

#ifndef USE_RDTSC
#include <cpucycles.h>
#endif

#endif

#include <time.h>
#include <string.h>

static int char2int(char input)
{
    if ((input >= '0') && (input <= '9'))
        return input - '0';
    if ((input >= 'A') && (input <= 'F'))
        return input - 'A' + 10;
    if ((input >= 'a') && (input <= 'f'))
        return input - 'a' + 10;
    return 0;
}

static void byte2hex(char *ptr,unsigned char ch)
{
    int t=ch/16;
    int b=ch%16;
    if (t<10)
    	ptr[0]='0'+t;
    else
    	ptr[0]='a'+(t-10);
    if (b<10)
    	ptr[1]='0'+b;
    else
    	ptr[1]='a'+(b-10);    	
}

// Convert a byte array to a hex string 
static void toHex(const char *src, char *dst)
{
    int i;
    for (i = 0; i < Nbytes; i++)
    {
        unsigned char ch = src[i];
        byte2hex(&dst[i * 2],ch);
    }
    dst[2*Nbytes]='\0';
}

// Convert from a hex string to byte array 
static void fromHex(const char *src, char *dst)
{
    int i,lz,len=0;
    char pad[2*Nbytes];
    while (src[len]!=0) len++;
    lz=2*Nbytes-len;
    if (lz<0) lz=0;
    for (i=0;i<lz;i++) pad[i]='0';  // pad with leading zeros
    for (i=lz;i<2*Nbytes;i++) pad[i]=src[i-lz];

    for (i=0;i<Nbytes;i++)
    {
        dst[i] = (char2int(pad[2*i]) * 16) + char2int(pad[2*i + 1]);
    }
}

// reverse bytes. Useful when dealing with little-endian formats
static void reverse(char *w)
{
    int i;
    for (i = 0; i < (Nbytes/2); i++) {
        unsigned char ch = w[i];
        w[i] = w[Nbytes - i - 1]; 
        w[Nbytes - i - 1] = ch; 
    } 
}

// output a modulo number in hex
/*
static void output(spint *x) {
    char b[Nbytes+1];
    char buff[(2*Nbytes)+1];
    modexp(x,b);
    toHex(b,buff);
    puts(buff);
}
*/
// Describe Montgomery Curve parameters

#ifdef X25519
#define A24 121665  // Montgomery curve constant (A-2)/4
#define COF 3       // Montgomery curve cofactor = 2^cof (2 or 3)
#define GENERATOR 9
#endif

#ifdef X448
#define A24 39081   // Montgomery curve constant (A-2)/4
#define COF 2       // Montgomery curve cofactor = 2^cof (2 or 3)
#define GENERATOR 5
#endif

// clamp input - see RFC7748
static void clamp(char *bk) {
    int s=(8-(Nbits%8))%8;
    bk[0]&=-(1<<COF);
    char mask=(unsigned char)(0xffu>>s);
    bk[Nbytes-1]&=mask;
    bk[Nbytes-1]|=(unsigned char)(0x80u>>s);
}

// return nth bit of byte array
static int bit(int n,const char *a) {
    return (int)((a[n/8u]&((unsigned char)1u<<(n%8u)))>>(n%8u));
}

static char mask() {
    int r=Nbits%8;
    if (r==0) r=8;
    return (char)((1<<r)-1);
}

// RFC7748 - Montgomery curve
// bv=bk*bu, bu,bv are x coordinates on elliptic curve
void rfc7748(const char *bk1, const char *bu1,char *bv1,const char *bk2, const char *bu2,char *bv2) {
    int i;
    uint32_t kt[2];
    spint swap,sb;
    char ck1[Nbytes],ck2[Nbytes];
    char cu1[Nbytes],cu2[Nbytes];
    spint u[Nlimbs]; spint x1[Nlimbs]; spint x2[Nlimbs]; spint x3[Nlimbs]; spint z2[Nlimbs]; spint z3[Nlimbs];
    spint A[Nlimbs]; spint B[Nlimbs]; spint AA[Nlimbs]; spint BB[Nlimbs]; spint C[Nlimbs]; spint D[Nlimbs]; spint E[Nlimbs];
    char msk=mask();

    for (i=0;i<Nbytes;i++) {
        ck1[i]=bk1[i];
        cu1[i]=bu1[i];
        ck2[i]=bk2[i];
        cu2[i]=bu2[i];	
    }

    reverse(cu1);  // convert from little to big endian
    cu1[0]&=msk;  // Mask most significant bitS in the final byte
    reverse(cu2);  // convert from little to big endian
    cu2[0]&=msk;  // Mask most significant bitS in the final byte    
    
// clamp input
    clamp(ck1);
    clamp(ck2); 

// import into internal representation
    modimp(cu1,cu2,u);

    modcpy(u,x1);  // x_1=u
    modone(x2);    // x_2=1
    modzer(z2);    // z_2=0
    modcpy(u,x3);  // x_3=u
    modone(z3);    // z_3=1

    swap=vdup_n_u32(0);
    for (i=Nbits-1;i>=0;i--)
    {
        kt[0]=(uint32_t)bit(i,ck1);
	kt[1]=(uint32_t)bit(i,ck2);
	sb=vld1_u32(kt);
	swap=veor_u32(swap,sb);

        modcsw(swap,x2,x3);
        modcsw(swap,z2,z3);
        swap=sb;
            
        modadd(x2,z2,A);        // A = x_2 + z_2
        modsqr(A,AA);           // AA = A^2
        modsub(x2,z2,B);        // B = x_2 - z_2
        modsqr(B,BB);           // BB = B^2

        modsub(AA,BB,E);        // E = AA - BB
        modadd(x3,z3,C);        // C = x_3 + z_3
        
        modsub(x3,z3,D);        // D = x_3 - z_3
        modmul(D,A,D);          // DA = D * A
        modmul(C,B,C);          // CB = C * B
 
        modadd(D,C,x3); modsqr(x3,x3);    // x_3 = (DA + CB)^2
        
        modsub(D,C,z3); modsqr(z3,z3); modmul(z3,x1,z3);  // z_3 = x_1 * (DA - CB)^2
        modmul(AA,BB,x2);       // x_2 = AA * BB
        modmli(E,A24,z2);        
        modadd(z2,AA,z2); modmul(z2,E,z2);   // z_2 = E * (AA + a24 * E)
    }
    modcsw(swap,x2,x3);
    modcsw(swap,z2,z3);

    modpro(z2,A);       
    modinv(z2,A,z2);    // sufficient for twist secure curves like X25519 and X448 

    modmul(x2,z2,x2);   

    modexp(x2,bv1,bv2);
    reverse(bv1); // convert to little endian
    reverse(bv2);
}

// a test vector for x25519 or x448 from RFC7748
int main()
{
    char sv1[(Nbytes*2)+1];
    char sv2[(Nbytes*2)+1];
    uint64_t start,fin;
    uint16_t rnd=1;  // crude random number generator
    clock_t begin;
    int i,elapsed;
    char bk1[Nbytes],bv1[Nbytes];
    char bu1[Nbytes]={};
    char bk2[Nbytes],bv2[Nbytes];
    char bu2[Nbytes]={};    

#if defined(X25519) || defined(X448)
#ifdef X25519
    const char *sk1=(const char *)"77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a";
    const char *sk2=(const char *)"5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb";    
#endif
#ifdef X448
    const char *sk1=(const char *)"9a8f4925d1519f5775cf46b04b5800d4ee9ee8bae8bc5565d498c28dd9c9baf574a9419744897391006382a6f127ab1d9ac2d8c0a598726b";
    const char *sk2=(const char *)"1c306a7ac2a0e2e0990b294470cba339e6453772b075811d8fad0d1d6927c120bb5ee8972b0d3e21374c9c921b09d1b0366f10b65173992d";    
#endif

    bu1[0]=bu2[0]=GENERATOR;
// convert to byte array
    fromHex(sk1,bk1);
    fromHex(sk2,bk2);
    rfc7748(bk1,bu1,bv1,bk2,bu2,bv2);

// convert to Hex
    toHex(bv1,sv1);
    printf("RFC7748 Test Vector 1\n");
    puts(sk1); 
    puts(sv1); 
    
    toHex(bv2,sv2);
    printf("RFC7748 Test Vector 2\n");
    puts(sk2); 
    puts(sv2);    
#endif

#ifdef COUNT_CLOCKS
#ifdef USE_RDTSC
    start=__rdtsc();
#else
    start=cpucycles();
#endif    
#endif


    for (i=0;i<Nbytes;i++) {
	bu1[i]=bu2[i]=0;
        rnd=5*rnd+1; bk1[i]=rnd%256;
    }
    bu1[0]=bu2[0]=GENERATOR;    
    for (i-0;i<Nbytes;i++) {
	rnd=5*rnd+1; bk2[i]=rnd%256;
    }
    begin=clock();
    for (i=0;i<5000;i++) {
        rfc7748(bk1,bu1,bv1,bk2,bu2,bv2);
        rfc7748(bk1,bv1,bu1,bk2,bv2,bu2);
    }
    elapsed=100*(clock() - begin) / CLOCKS_PER_SEC;
#ifdef COUNT_CLOCKS
#ifdef USE_RDTSC
    fin=__rdtsc();
#else
    fin=cpucycles();
#endif
    printf("Clock cycles per point multiplication= %d\n",(int)((fin-start)/20000ULL));
#endif

    printf("Microseconds per point multiplication= %d\n",elapsed/2);
    toHex(bu1,sv1);
    puts(sv1);
    toHex(bu2,sv2);
    puts(sv2);    

// do a double D-H key exchange

    char alice1[Nbytes],bob1[Nbytes],apk1[Nbytes],bpk1[Nbytes],ssa1[Nbytes],ssb1[Nbytes];
    char alice2[Nbytes],bob2[Nbytes],apk2[Nbytes],bpk2[Nbytes],ssa2[Nbytes],ssb2[Nbytes];
    for (i=0;i<Nbytes;i++)
    {
        rnd=5*rnd+1; alice1[i]=rnd%256;
        rnd=5*rnd+1; bob1[i]=rnd%256;
        apk1[i]=0; bpk1[i]=0;
    }
    for (i=0;i<Nbytes;i++)
    {
        rnd=5*rnd+1; alice2[i]=rnd%256;
        rnd=5*rnd+1; bob2[i]=rnd%256;
        apk2[i]=0; bpk2[i]=0;
    }    
    apk1[0]=bpk1[0]=GENERATOR;
    apk2[0]=bpk2[0]=GENERATOR;    
    rfc7748(alice1,apk1,apk1,alice2,apk2,apk2);
    rfc7748(bob1,bpk1,bpk1,bob2,bpk2,bpk2);

    rfc7748(alice1,bpk1,ssa1,alice2,bpk2,ssa2);
    rfc7748(bob1,apk1,ssb1,bob2,apk2,ssb2);
    toHex(ssa1,sv1);
    printf("Alice first shared secret\n");
    puts(sv1);
    toHex(ssa2,sv2);
    printf("Alice second shared secret\n");
    puts(sv2);    
    
    toHex(ssb1,sv1);
    printf("Bob's first shared secret\n");
    puts(sv1);
    toHex(ssb2,sv2);
    printf("Bob's second shared secret\n");
    puts(sv2);

}
