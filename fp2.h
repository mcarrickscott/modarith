//
// C FP2 functions
// M.Scott 2026
//

#ifndef FP2_H
#define FP2_H

typedef struct {
	spint a[Nlimbs];
	spint b[Nlimbs];
} fp2;

// hopefully most of the fp2 operations that will be required

extern void fp2_zero(fp2 *x);                // x=0
extern void fp2_seti(const unsigned int r,const unsigned int i,fp2 *x);  // x={real, imaginary}
extern void fp2_set(spint *r,spint *i,fp2 *x); // x={real, imaginary}
extern void fp2_add(fp2 *x,fp2 *y,fp2 *z);  // z=x+y
extern void fp2_sub(fp2 *x,fp2 *y,fp2 *z);  // z=x-y
extern void fp2_mul(fp2 *x,fp2 *y,fp2 *z);  // z=x*y
extern void fp2_sqr(fp2 *x,fp2 *z);         // z=x*x
extern void fp2_conj(fp2 *x);               // z is conjugate of x  (Frobenius)
extern void fp2_neg(fp2 *x,fp2 *z);         // z=-x
extern void fp2_half(fp2 *x);                // half x
extern void fp2_imul(fp2 *x,int n,fp2 *z);  // z=n*x
extern void fp2_pmul(fp2 *x,spint *n,fp2 *z); // z=n*x
extern void fp2_times_i(fp2 *x,fp2 *z);     // x = z times i
extern int fp2_is0(fp2 *x);                 // test x=0? 
extern int fp2_is1(fp2 *x);                 // test x=1? 
extern int fp2_split(fp2 *x,spint *r, spint *i); // split z into "real" and "imaginary" parts
extern int fp2_cmp(fp2 *x,fp2 *y);          // compare for equality
extern void fp2_inv(fp2 *x,const spint *h,fp2 *z) ; // z=1/x (with possible hint, or NULL)
extern int fp2_qr(const spint *h,fp2 *x);          // test if x is quadratic residue (with possible hint, or NULL)
extern void fp2_cpy(fp2 *x,fp2 *z);          // z=x
extern int fp2_sign(fp2 *x);                 // extract "sign" of x 
extern void fp2_size(fp2 *x,spint *z);       // z=|x|
extern void fp2_cmv(int d,fp2 *x,fp2 *z);    // move x to z if d=1
extern void fp2_csw(int d,fp2 *x,fp2 *z);    // swap x and z if d=1
extern void fp2_sqrt(fp2 *x,const spint *h,fp2 *z); // z=square root of x (with possible hint, or NULL)
extern void fp2_export(const fp2 *x,char *outr,char *outi); // export to byte array
extern void fp2_import(const char *inr,const char *ini,fp2 *x);  // import from byte array

#endif