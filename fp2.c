//
// C Fp2 functions
// M.Scott 2026
//
// To compile:-
// gcc -O2 fp2.c -o fp2

// run monty.py or pseudo.py script for chosen modulus
// This code assumes modulus is prime = 3 mod 4


/*** Insert automatically generated code for modulus field.c here ***/



/*** End of automatically generated code ***/

#include "fp2.h"

void fp2_zero(fp2 *x)
{
	modzer(x->a); modzer(x->b);
}

void fp2_seti(const unsigned int r,const unsigned int i,fp2 *x)
{
	modint(r,x->a); modint(i,x->b);
}

void fp2_set(fp r,fp i,fp2 *x)
{
	modcpy(r,x->a); modcpy(i,x->b);
}

void fp2_add(fp2 *x,fp2 *y,fp2 *z)
{
	modadd(x->a,y->a,z->a);
	modadd(x->b,y->b,z->b);
}

void fp2_sub(fp2 *x,fp2 *y,fp2 *z)
{
	modsub(x->a,y->a,z->a);
	modsub(x->b,y->b,z->b);
}

// NOTE: faster "lazy" method may apply if prime modulus is not "shaped"
void fp2_mul(fp2 *x,fp2 *y,fp2 *z)
{ 
	fp c,d,f;
	modmul(x->a,y->a,c);  // c=xa*ya
	modmul(x->b,y->b,d);  // d=xb*yb
	modaddf(y->a,y->b,f); // f=ya+yb
	modaddf(x->b,x->a,z->b); // e=xa+xb
	modsub(c,d,z->a);     // za=xa*ya-xb*yb
	modadd(c,d,c);        // c=xa*ya+xb*yb
	modmul(z->b,f,z->b);        // e=(xa+xb)*(ya+yb)
	modsub(z->b,c,z->b);     // zb=(xa+xb)*(ya+yb) - (xa*ya+xb*yb) = xa*yb+xb*ya
}

void fp2_sqr(fp2 *x,fp2 *z)
{
	fp c,d;
	modaddf(x->a,x->a,c);
	modaddf(x->a,x->b,d);
	modsubf(x->a,x->b,z->a);
	modmul(z->a,d,z->a);
	modmul(x->b,c,z->b);
}

void fp2_conj(fp2 *x)
{
	modneg(x->b,x->b);
}

void fp2_neg(fp2 *x,fp2 *z)
{
	modneg(x->a,z->a);
	modneg(x->b,z->b);
}

void fp2_half(fp2 *x)
{
	modhaf(x->a); modhaf(x->b);
}

// multiply by integer
void fp2_imul(fp2 *x,int n,fp2 *z)
{
	modmli(x->a,n,z->a);
	modmli(x->b,n,z->b);
}

// multiply by fp
void fp2_pmul(fp2 *x,fp n,fp2 *z)
{
	modmul(x->a,n,z->a);
	modmul(x->b,n,z->b);
}

// multiply by i
void fp2_times_i(fp2 *x,fp2 *z)
{
	fp t;
	modcpy(x->a,t);
	modneg(x->b,z->a);
	modcpy(t,z->b);
}

int fp2_is0(fp2 *x) 
{
	return (modis0(x->a) && modis0(x->b));
}

int fp2_is1(fp2 *x) 
{
	return (modis1(x->a) && modis0(x->b));
}

// test x less than y
int fp2_ltn(fp2 *x,fp2 *y)
{
	int m=modltn(x->b,y->b);
	int n=modcmp(x->b,y->b);
	int w=modltn(x->a,y->a);
	return (m|(n&w));
}

int fp2_real(fp2 *x,fp r)
{
	modcpy(x->a,r);
}

int fp2_imag(fp2 *x,fp i)
{
	modcpy(x->b,i);
}

int fp2_cmp(fp2 *x,fp2 *y)
{
	return (modcmp(x->a,y->a) && modcmp(x->b,y->b));
}

void fp2_inv(fp2 *x,const fp h,fp2 *z) 
{
	fp w1,w2;
	modsqr(x->a,w1);
	modsqr(x->b,w2);
	modadd(w1,w2,w1);
	modinv(w1,h,w1);
	modmul(x->a,w1,z->a);
	modneg(w1,w1);
	modmul(x->b,w1,z->b);
}

int fp2_qr(const fp h,fp2 *x) 
{
	fp w1,w2;
	modsqr(x->a,w1);
	modsqr(x->b,w2);
	modadd(w1,w2,w1);
	return modqr(h,w1);
}

void fp2_cpy(fp2 *x,fp2 *z)
{
	modcpy(x->a,z->a); modcpy(x->b,z->b);
}

int fp2_sign(fp2 *x)
{
	int p1=modsign(x->a);
	int p2=modsign(x->b);
	int u=modis0(x->a);
	return p1^((p1^p2)&u);
}

void fp2_cmv(int d,fp2 *x,fp2 *z)
{
	modcmv(d,x->a,z->a);
	modcmv(d,x->b,z->b);
}

void fp2_csw(int d,fp2 *x,fp2 *z)
{
	modcsw(d,x->a,z->a);
	modcsw(d,x->b,z->b);
}

// z=|x|
void fp2_size(fp2 *x,fp z)
{
	fp w1,w2;
	modsqr(x->a,w1);
	modsqr(x->b,w2);
	modadd(w1,w2,z);
}

void fp2_sqrt(fp2 *x,const fp h,fp2 *z)
{
	fp w1,w2,w3,w4,hint;
	fp2 nz;
	int sgn,qr;
	if (fp2_is0(x) || fp2_is1(x)) return;
	fp2_cpy(x,z);
	
	modsqr(z->b,w1);
	modsqr(z->a,w2);
	modadd(w1,w2,w1);
	modsqrt(w1,h,w1);   // sqrt(a^2+b^2)  - could make use of an input hint to avoid exp!
	modadd(w1,z->a,w2);
	modhaf(w2);

	modcpy(z->b,w1);
	modhaf(w1);
	
	modpro(w2,hint);
	qr=modqr(hint,w2);  // only exponentiation??

	modneg(hint,w3);
	modneg(w2,w4);

	modcmv(1-qr,w4,w2);
	modcmv(1-qr,w3,hint);

	modsqrt(w2,hint,z->a);
	modinv(w2,hint,w3);
	modmul(w3,z->a,w3);
	modmul(w3,w1,z->b);
	modcpy(z->a,w4);

	modcmv(1-qr,z->b,z->a);
	modcmv(1-qr,w4,z->b);
 
	sgn=fp2_sign(z);
	fp2_neg(z,&nz);
	fp2_cmv(sgn,&nz,z);
}

void fp2_export(const fp2 *x,char *outr,char *outi)
{
	modexp(x->a,outr);
	modexp(x->b,outi);
}

void fp2_import(const char *inr,const char *ini,fp2 *x)
{
	modimp(inr,x->a);
	modimp(ini,x->b);
}

// ------------------------------------

// Give it a good work out

#include <stdlib.h>

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

int main() 
{
	char buffr[Nbytes],buffi[Nbytes],hex[(Nbytes*2)+1];
	int sgn1,sgn2;
	fp2 x,y,z,w1,w2;
	fp hint,t;

#if PMOD8==1 || PMOD8==5
	printf("The prime modulus must be 3 mod 4\n");
	exit(0);
#endif

	fp2_seti(0x18,0x34,&x);
	fp2_seti(0x87,0x22,&y);

	for (int i=0;i<10000;i++)
	{
		fp2_add(&x,&y,&w1);
		fp2_sub(&x,&y,&w2);
		fp2_mul(&w1,&w2,&z);   // z=(x+y)(x-y)
		sgn1=fp2_sign(&z);
		fp2_sqr(&z,&z);        // z=[(x+y)(x-y)]^2

		fp2_size(&z,t);        // t=|x|
		modpro(t,hint);        // create a hint for z
		fp2_inv(&z,hint,&z);   // z=1/[(x+y)(x-y)]^2
		fp2_half(&z);
		fp2_add(&z,&z,&z);

		fp2_sqrt(&z,NULL,&z);  // z=1/[(x+y)(x-y)]
		fp2_inv(&z,NULL,&z);   // z=(x+y)(x-y)
		sgn2=fp2_sign(&z);
		if (sgn1!=sgn2) fp2_neg(&z,&z); // wrong root
		
		fp2_sqr(&y,&w1);       // w1=y^2
		fp2_add(&z,&w1,&z);    // z=x^2
		fp2_sqr(&x,&w2);  
		if (fp2_cmp(&z,&w2))
			printf("passed\n");
		else
			{printf("failed\n"); exit(0);}
		fp2_cpy(&w1,&x);
		fp2_cpy(&z,&y);
	}
	fp2_export(&z,buffr,buffi);
	toHex(buffr,hex);
	printf("r= %s\n",hex);
	toHex(buffi,hex);
	printf("i= %s\n",hex);
	return 0;
}
