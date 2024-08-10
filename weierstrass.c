// Weierstrass curve support 
// Use python scripts to generate code for NIST curves, or your own curves
//
// Mike Scott 16th July 2024
// TII
//
// code for 16/32/64-bit processor for NIST curve can be generated  by 
//
// python monty.py 16/32/64 NIST256
//

// make sure decoration and generic are both set to False in monty.py or pseudo.py

/*** Insert automatically generated code for modulus code.c here ***/



/*** End of automatically generated code ***/

#include "curve.h"

#define BYTES Nbytes
#define LIMBS Nlimbs
#define TOPBIT (8*sizeof(int)-1)

// define weierstrass curve here y^2=x^3-3x+B, that is B and prime order generator (x,y)
// define ZEROA if curve is y^2=x^3+B, and constant_B = 3B
#ifdef NUMS256W  // the way it should have been done... see https://csrc.nist.gov/csrc/media/events/workshop-on-elliptic-curve-cryptography-standards/documents/papers/session4-costello-craig.pdf
#ifdef MULBYINT
#define CONSTANT_B 152961
#define CONSTANT_X 2
#else
#error "Only allowed if MULBYINT defined in Fp. Otherwise use make.py to precompute constant"
#endif
#endif

// utility make.py can be used to generate these constants

#ifdef NIST256  // the way it was done....
// get constant d, generator point (x,y)
#if Radix == 52
const spint constant_b[5]={0xdf6229c4bddfd,0xca8843090d89c,0x212ed6acf005c,0x83415a220abf7,0xc30061dd4874};
const spint constant_x[5]={0x30d418a9143c1,0xc4fedb60179e7,0x62251075ba95f,0x5c669fb732b77,0x8905f76b5375};
const spint constant_y[5]={0x5357ce95560a8,0x43a19e45cddf2,0x21f3258b4ab8e,0xd8552e88688dd,0x571ff18a5885};
#endif
#if Radix ==  29
const spint constant_b[9]={0x1897bbfb,0x1cdf6229,0x18486c4,0x1732821,0x1dad59e0,0xabf7212,0x1a06d110,0x17721d20,0x8600c3};
const spint constant_x[9]={0x15228783,0x730d418,0xdb00bcf,0x57f11fb,0xa20eb75,0x12b77622,0x330fdb9,0x1af4dd57,0x120bee};
const spint constant_y[9]={0x12aac150,0x125357ce,0xf22e6ef,0xe390e86,0x64b1695,0x88dd21f,0x2a97443,0x2962176,0xae3fe3};
#endif
#endif

#ifdef SECP256K1
#define ZEROA
#ifdef MULBYINT
#define CONSTANT_B 7
#define CONSTANT_B3 21 // needed if ZEROA
#else
#if Radix == 29 // monty
const spint constant_b[9]={0x356e0,0x700,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
const spint constant_b3[9]={0xa04a0,0x1500,0x0,0x0,0x0,0x0,0x0,0x0,0x0}; // needed if ZEROA
#endif

#endif

#if Radix == 52 // pseudo
const spint constant_x[5]={0x2815b16f81798,0xdb2dce28d959f,0xe870b07029bfc,0xbbac55a06295c,0x79be667ef9dc};
const spint constant_y[5]={0x7d08ffb10d4b8,0x48a68554199c4,0xe1108a8fd17b4,0xc4655da4fbfc0,0x483ada7726a3};
#endif
#if Radix ==  29 // monty
const spint constant_x[9]={0xfc45b63,0x162e5ae0,0x336deb9,0xa54ca6f,0x538463c,0xc033fd1,0x44bcfa4,0xfa4227d,0x303cc8};
const spint constant_y[9]={0x1b75dba9,0x1ea6d39b,0xe326d8a,0x175747c7,0x26d1bf8,0x19aac19c,0xb385b5a,0x1f52960b,0xe7f0a3};
#endif

#endif

// return 1 if b==c, no branching 
static int teq(int b, int c)
{
    int x = b ^ c;
    x -= 1; // if x=0, x now -1
    return ((x >> TOPBIT) & 1);
}

// copy point
void ecncpy(point *Q,point *P)
{
    modcpy(Q->x,P->x);
    modcpy(Q->y,P->y);
    modcpy(Q->z,P->z);
}

// negate P
void ecnneg(point *P)
{
    modneg(P->y,P->y);
}

// add Q to P
// complete formulae from https://eprint.iacr.org/2015/1060
void ecnadd(point *Q,point *P)
{
    spint B[Nlimbs],T0[Nlimbs],T1[Nlimbs],T2[Nlimbs],T3[Nlimbs],T4[Nlimbs];

    modmul(P->x,Q->x,T0);   // 1
    modmul(P->y,Q->y,T1);   // 2
    modmul(P->z,Q->z,T2);   // 3

    modadd(P->x,P->y,T3);   // 4
    modadd(Q->x,Q->y,T4);   // 5
    modmul(T3,T4,T3);       // 6

    modadd(T0,T1,T4);       // 7
    modsub(T3,T4,T3);       // 8
    modadd(P->y,P->z,T4);   // 9

    modadd(Q->y,Q->z,B); // 10
    modmul(T4,B,T4);     // 11
    modadd(T1,T2,B);        // 12  use B

    modsub(T4,B,T4);        // 13
    modadd(P->x,P->z,P->x); // 14
    modadd(Q->z,Q->x,P->y); // 15

    modmul(P->x,P->y,P->x); // 16
    modadd(T0,T2,P->y);     // 17
    modsub(P->x,P->y,P->y); // 18

#ifdef ZEROA
    modadd(T0,T0,P->x);   // 19
    modadd(T0,P->x,T0);   // 20

#ifdef CONSTANT_B3
#if CONSTANT_B3>0
    modmli(T2,CONSTANT_B3,T2);      // 21
    modmli(P->y,CONSTANT_B3,P->y);  // 24
#else
    modmli(T2,-CONSTANT_B3,T2);  modneg(T2,T2);    // 21
    modmli(P->y,-CONSTANT_B3,P->y);  modneg(P->y,P->y);// 24
#endif
#else
    modcpy(constant_b3,B);
    modmul(T2,B,T2);      // 21
    modmul(P->y,B,P->y);  // 24
#endif
    modadd(T1,T2,P->z);   // 22
    modsub(T1,T2,T1);     // 23

    modmul(P->y,T4,P->x); // 25
    modmul(T3,T1,T2);     // 26
    modsub(T2,P->x,P->x); // 27
    modmul(P->y,T0,P->y); // 28
    modmul(T1,P->z,T1);   // 29
    modadd(P->y,T1,P->y); // 30
    modmul(T0,T3,T0);     // 31
    modmul(P->z,T4,P->z); // 32
    modadd(P->z,T0,P->z); // 33
#else

#ifdef CONSTANT_B
#if CONSTANT_B>0
    modmli(T2,CONSTANT_B,P->z);      //19
    modsub(P->y,P->z,P->x); // 20
    modmli(P->y,CONSTANT_B,P->y);    // 25
#else
    modmli(T2,-CONSTANT_B,P->z);      //19
    modadd(P->y,P->z,P->x); // 20
    modmli(P->y,-CONSTANT_B,P->y); modneg(P->y,P->y);    // 25
#endif
#else
    modcpy(constant_b,B);
    modmul(B,T2,P->z);      //19
    modsub(P->y,P->z,P->x); // 20
    modmul(P->y,B,P->y);    // 25
#endif
    modadd(P->x,P->x,P->z); // 21

    modadd(P->x,P->z,P->x); // 22
    modsub(T1,P->x,P->z);   // 23
    modadd(P->x,T1,P->x);   // 24

    modadd(T2,T2,T1);       // 26
    modadd(T2,T1,T2);       // 27

    modsub(P->y,T2,P->y);   // 28
    modsub(P->y,T0,P->y);   // 29
    modadd(P->y,P->y,T1);   // 30

    modadd(P->y,T1,P->y);   // 31
    modadd(T0,T0,T1);       // 32
    modadd(T0,T1,T0);       // 33

    modsub(T0,T2,T0);       // 34
    modmul(T4,P->y,T1);     // 35
    modmul(T0,P->y,T2);     // 36

    modmul(P->x,P->z,P->y); // 37
    modadd(P->y,T2,P->y);   // 38
    modmul(P->x,T3,P->x);   // 39

    modsub(P->x,T1,P->x);   // 40
    modmul(P->z,T4,P->z);   // 41
    modmul(T3,T0,T1);       // 42

    modadd(P->z,T1,P->z);   // 43
#endif
}

// subtract Q from P
void ecnsub(point *Q,point *P)
{
    point W;
    ecncpy(Q,&W); ecnneg(&W);
    ecnadd(&W,P);
}

// double P
// complete formuale from https://eprint.iacr.org/2015/1060
void ecndbl(point *P)
{
    spint B[Nlimbs],T0[Nlimbs],T1[Nlimbs],T2[Nlimbs],T3[Nlimbs],T4[Nlimbs];

#ifdef ZEROA
    modsqr(P->y,T0); // 1
    modadd(T0,T0,T3);  // 2 T3=Z3
    modadd(T3,T3,T3);  // 3
    modadd(T3,T3,T3);  // 4
    modmul(P->x,P->y,T4); // 16 T4=X*Y
    modmul(P->y,P->z,T1); // 5
    modsqr(P->z,T2);      // 6
#ifdef CONSTANT_B3
#if CONSTANT_B3>0
    modmli(T2,CONSTANT_B3,T2);      // 7
#else
    modmli(T2,-CONSTANT_B3,T2); modneg(T2,T2);
#endif
#else
    modcpy(constant_b3,B);
    modmul(T2,B,T2);      // 7
#endif
    modmul(T2,T3,P->x);   // 8
    modadd(T0,T2,P->y);   // 9
    modmul(T3,T1,P->z); // 10
    modadd(T2,T2,T1);     // 11
    modadd(T2,T1,T2);     // 12
    modsub(T0,T2,T0);     // 13
    modmul(P->y,T0,P->y); // 14
    modadd(P->y,P->x,P->y); // 15
    modmul(T0,T4,P->x);    // 17
    modadd(P->x,P->x,P->x); // 18
#else
    modsqr(P->x,T0); // 1
    modsqr(P->y,T1); // 2
    modsqr(P->z,T2); // 3

    modmul(P->x,P->y,T3); // 4
    modmul(P->y,P->z,T4); // 28
    modadd(T3,T3,T3);     // 5
    modmul(P->z,P->x,P->z); // 6

    modadd(P->z,P->z,P->z); // 7

#ifdef CONSTANT_B
#if CONSTANT_B>0 
    modmli(T2,CONSTANT_B,P->y); // 8
    modsub(P->y,P->z,P->y); // 9
    modmli(P->z,CONSTANT_B,P->z); // 18
#else
    modmli(T2,-CONSTANT_B,P->y); modneg(P->y,P->y); // 8
    modsub(P->y,P->z,P->y); // 9
    modmli(P->z,-CONSTANT_B,P->z);  modneg(P->z,P->z); // 18
#endif
#else
    modcpy(constant_b,B);
    modmul(T2,B,P->y); // 8
    modsub(P->y,P->z,P->y); // 9
    modmul(P->z,B,P->z); // 18
#endif

    modadd(P->y,P->y,P->x); // 10
    modadd(P->y,P->x,P->y); // 11
    modsub(T1,P->y,P->x); // 12

    modadd(P->y,T1,P->y); // 13
    modmul(P->y,P->x,P->y); // 14
    modmul(P->x,T3,P->x); // 15

    modadd(T2,T2,T3); // 16
    modadd(T2,T3,T2); // 17

    modsub(P->z,T2,P->z); // 19
    modsub(P->z,T0,P->z); // 20
    modadd(P->z,P->z,T3); // 21

    modadd(P->z,T3,P->z); // 22
    modadd(T0,T0,T3); // 23
    modadd(T0,T3,T0); // 24

    modsub(T0,T2,T0); // 25
    modmul(T0,P->z,T0); // 26
    modadd(P->y,T0,P->y); // 27

    modadd(T4,T4,T4); // 29
    modmul(P->z,T4,P->z); // 30

    modsub(P->x,P->z,P->x); // 31
    modmul(T4,T1,P->z); // 32
    modadd(P->z,P->z,P->z); // 33

    modadd(P->z,P->z,P->z); // 34
#endif
}

// set to infinity
void ecninf(point *P)
{
    modzer(P->x);
    modone(P->y);
    modzer(P->z);
}

// test for infinity
int ecnisinf(point *P)
{
    return (modis0(P->x) && modis0(P->z));
}

// set to affine
void ecnaffine(point *P)
{
    spint I[Nlimbs];
    modinv(P->z,NULL,I);
    modone(P->z);
    modmul(P->x,I,P->x);
    modmul(P->y,I,P->y);
}

// move Q to P if d=1
void ecncmv(int d,point *Q,point *P)
{
    modcmv(d,Q->x,P->x);
    modcmv(d,Q->y,P->y);
    modcmv(d,Q->z,P->z);
}

// return 1 if equal, else 0
int ecncmp(point *P,point *Q)
{
    spint a[Nlimbs],b[Nlimbs];
    modmul(P->x,Q->z,a);
    modmul(Q->x,P->z,b);
    if (!modcmp(a,b)) return 0;
    modmul(P->y,Q->z,a);
    modmul(Q->y,P->z,b);
    if (!modcmp(a,b)) return 0;
    return 1;
}

// extract (x,y) from point, if y is NULL compress and just return x and sign of y
int ecnget(point *P,char *x,char *y)
{
    spint X[Nlimbs],Y[Nlimbs];
    ecnaffine(P);
    if (x!=NULL)
    {
        modcpy(P->x,X);
        modexp(X,x);
    }
    if (y!=NULL)
    {
        modcpy(P->y,Y);
        modexp(Y,y);
    }
    if (y==NULL) return modsign(P->y);
    if (x==NULL) return modsign(P->x);
    return 0;
}
/*
int ecngetxyz(point *P,char *x,char *y,char *z)
{
    spint X[Nlimbs],Y[Nlimbs];
    modexp(P->x,x);
    modexp(P->y,y);
    modexp(P->z,z);
    return 0;
}
*/
// weierstrass set point function
// sets P=O if point not on curve
// if y!=NULL tries to set (x,y)
// if y==NULL calculates y (decompresses x) and selects sign from s=0/1
static void setxy(int s,const spint *x,const spint *y,point *P)
{
    spint T[Nlimbs],V[Nlimbs],H[Nlimbs];
    if (x==NULL)
    {
        ecninf(P);
        return;
    }
    modcpy(x,P->x);
    modsqr(x,V);
    modmul(V,x,V); // x^3
#ifndef ZEROA
    modsub(V,x,V);
    modsub(V,x,V);
    modsub(V,x,V); // x^3-3x
#endif
#ifdef CONSTANT_B
#if CONSTANT_B>0
    modint(CONSTANT_B,T);
    modadd(V,T,V); // V=1+dx^2
#else
    modint(-CONSTANT_B,T);
    modsub(V,T,V); // V=1-dx^2
#endif      
#else
    modadd(V,constant_b,V);
#endif
    if (y!=NULL)
    {
        modsqr(y,T);
        if (modcmp(T,V)) {
            modcpy(y,P->y);
            modone(P->z);
        } else {
            ecninf(P);
        }
        return;
    }
    modpro(V,H);
    if (!modqr(H,V))
    { // point not on curve
        ecninf(P);
        return;
    }
    modsqrt(V,H,P->y);
    int d=(modsign(P->y)-s)&1;
    modneg(P->y,T);
    modcmv(d,T,P->y);
    modone(P->z);
}    

// multiply point by small curve cofactor (here assumed to be 1)
void ecncof(point *P)
{}

// api visible version, x and y are big endian byte arrays
void ecnset(int s,const char *x,const char *y,point *P)
{
    spint X[Nlimbs],Y[Nlimbs];
    if (x!=NULL && y!=NULL)
    {
        modimp(x,X);
        modimp(y,Y);
        setxy(s,X,Y,P);
        return;
    }
    if (x!=NULL)
    {
        modimp(x,X);
        setxy(s,X,NULL,P);
    }
    if (y!=NULL)
    {
        modimp(y,Y);
        setxy(s,NULL,X,P);
    }
}

// set generator
void ecngen(point *P)
{
#ifdef CONSTANT_X
    spint X[Nlimbs];
    modint(CONSTANT_X,X);
    setxy(0,X,NULL,P);
#else
    setxy(0,constant_x,constant_y,P);
#endif
}

// select point from precomputed array in constant time
static void select(int b,point W[],point *P)
{
    point MP;  
    int m = b >> TOPBIT;
    int babs = (b ^ m) - m;

    ecncmv(teq(babs, 0),&W[0],P); // conditional move
    ecncmv(teq(babs, 1),&W[1],P); // conditional move
    ecncmv(teq(babs, 2),&W[2],P); // conditional move
    ecncmv(teq(babs, 3),&W[3],P); // conditional move
    ecncmv(teq(babs, 4),&W[4],P); // conditional move
    ecncmv(teq(babs, 5),&W[5],P); // conditional move
    ecncmv(teq(babs, 6),&W[6],P); // conditional move
    ecncmv(teq(babs, 7),&W[7],P); // conditional move
    ecncmv(teq(babs, 8),&W[8],P); // conditional move
    
    ecncpy(P,&MP);
    ecnneg(&MP);  // minus P
    ecncmv((int)(m & 1),&MP,P);
}

// convert to double naf form
static void dnaf(const char *e,const char *f, char *w)
{
    int i,j,t;
    unsigned char ce=0;
    unsigned char cf=0;
    unsigned char m,n,p,q;
    for (i=0;i<Nbytes;i++)
    {
        m=n=e[Nbytes-i-1];
        t=3*(int)n+ce;
        ce=(unsigned char)(t>>8);
        n=(unsigned char)(t&0xff);
        p=q=f[Nbytes-i-1];
        t=3*(int)q+cf;
        cf=(unsigned char)(t>>8);
        q=(unsigned char)(t&0xff);
        for (j=0;j<8;j++)
        {
            w[8*i+j]=(n&1)-(m&1)+3*((q&1)-(p&1));
            n>>=1; m>>=1; p>>=1; q>>=1;
        }
    }
    for (j=0;j<8;j++)
    {
        w[8*Nbytes+j]=(ce&1)+3*(cf&1);
        ce>>=1; cf>>=1;
    }
}

// multiply point by scalar
// constant time
void ecnmul(const char *e,point *P) 
{
    int i,j;
    point Q,W[9];
    signed char w[2*Nbytes+1];

    ecninf(&W[0]);                         // O
    ecncpy(P,&W[1]);                       // P
    ecncpy(P,&W[2]); ecndbl(&W[2]);        // 2P
    ecncpy(&W[2],&W[3]); ecnadd(P,&W[3]);  // 3P
    ecncpy(&W[2],&W[4]); ecndbl(&W[4]);    // 4P
    ecncpy(&W[4],&W[5]); ecnadd(P,&W[5]);  // 5P
    ecncpy(&W[3],&W[6]); ecndbl(&W[6]);    // 6P
    ecncpy(&W[6],&W[7]); ecnadd(P,&W[7]);  // 7P
    ecncpy(&W[4],&W[8]); ecndbl(&W[8]);    // 8P

// convert exponent to signed digit
    for (i=j=0;i<Nbytes;i++,j+=2)
    {
        char c=e[Nbytes-i-1];
        w[j]=c&0xf;
        w[j+1]=(c>>4)&0xf;
    }
    w[2*Nbytes]=0;
    for (j=0;j<2*Nbytes;j++)
    {
        int t=7-w[j];
        int m=(t>>4)&1;
        w[j]-=(m<<4);
        w[j+1]+=m;
    }

//    printf("w= ");
//    for (i=0;i<2*Nbytes+1;i++) printf(" %d",(int)w[i]);
//    printf("\n");

    select(w[2*Nbytes],W,P);
    for (i = 2*Nbytes - 1; i >= 0; i--)
    {
        select(w[i],W,&Q);
        ecndbl(P);
        ecndbl(P);
        ecndbl(P);
        ecndbl(P);
        ecnadd(&Q,P);
    }
}

// double point multiplication R=eP+fQ
// not constant time
void ecnmul2(const char *e,point *P,const char *f,point *Q,point *R)
{
    int i;
    point T,W[5];
    signed char w[8*Nbytes+8];
    ecninf(&W[0]);     // O
    ecncpy(P,&W[1]);   // P
    ecncpy(Q,&W[3]);   // Q
    ecncpy(Q,&W[2]); ecnsub(P,&W[2]);  // Q-P
    ecncpy(Q,&W[4]); ecnadd(P,&W[4]);  // Q+P

    dnaf(e,f,w);

    i=8*Nbytes+7;
    while (w[i]==0) i--; // ignore leading zeros
    ecninf(R);
    while (i>=1)
    {
        ecndbl(R);
        if (w[i]!=0) {
            select(w[i],W,&T);
            ecnadd(&T,R);
        }
        i--;
    }
}
