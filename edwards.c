// Edwards curve support 
// Use python scripts to generate code for C25519 or C448, or your own curve
//
// Mike Scott 16th July 2024
// TII
//
// code for 16/32/64-bit processor for C25519 curve can be generated  by 
//
// python pseudo.py 16/32/64 C25519
// or
// python monty.py 16/32/64 C25519
//
// code for 16/32/64-bit processor for X448 curve can be generated  by
//
// python monty.py 16/32/64 C448

// make sure decoration and generic are both set to False in monty.py or pseudo.py

/*** Insert automatically generated code for modulus code.c here ***/



/*** End of automatically generated code ***/

#include "curve.h"

// define edwards curve here ax^2+y^2 = 1+dx^2y^2, that is d, cofactor and prime order generator (x,y)

#ifdef NUMS256E  // the way it should have been done...
#ifdef MULBYINT
#define CONSTANT_D -15342   
#define CONSTANT_X 34
#else
#error "Only allowed if MULBYINT defined in Fp. Otherwise use make.py to precompute constant"
#endif

#define COF 2
#endif

// utility make.py can be used to generate these constants

#ifdef Ed25519  // the way it was done....
// get constant d, generator point (x,y)
#if Radix == 51
const spint constant_d[5]={0x34dca135978a3,0x1a8283b156ebd,0x5e7a26001c029,0x739c663a03cbb,0x52036cee2b6ff};
const spint constant_x[5]={0x62d608f25d51a,0x412a4b4f6592a,0x75b7171a4b31d,0x1ff60527118fe,0x216936d3cd6e5};
const spint constant_y[5]={0x6666666666658,0x4cccccccccccc,0x1999999999999,0x3333333333333,0x6666666666666};
#endif
#if Radix ==  29
const spint constant_d[9]={0x135978a3,0xf5a6e50,0x10762add,0x149a82,0x1e898007,0x3cbbbc,0x19ce331d,0x1dc56dff,0x52036c};
const spint constant_x[9]={0xf25d51a,0xab16b04,0x969ecb2,0x198ec12a,0xdc5c692,0x1118feeb,0xffb0293,0x1a79adca,0x216936};
const spint constant_y[9]={0x6666658,0x13333333,0x19999999,0xccccccc,0x6666666,0x13333333,0x19999999,0xccccccc,0x666666};
#endif
#define COF 3
#define NEGA   // a is -1
#endif

#ifdef Ed448
// get constant d, generator point (x,y), converted to n-residue form
#ifdef MULBYINT
#define CONSTANT_D -39081
#else
#error "Only allowed if MULBYINT defined in Fp. Otherwise use make.py to precompute constant"
#endif
#if Radix == 56
const spint constant_x[8]={0x420685f0ea8836,0x35bf93b17aa383,0xb7bc2914f8fe6d,0xe44cd37ab765fa,0x34f39b1b69235e,0x44d6fb9be886a8,0xee96c7295e6eb4,0xd16ef0905d88b9};
const spint constant_y[8]={0xd81f4fba184177,0xac119c79a99632,0xda8e9ac23c2104,0x416ef259fc5486,0x46ff5902c1cc32,0x4fa9dd01223251,0xa1f0e6acaf9471,0x65f7687a33ab50};
//const spint constant_x[8]={0xaaaaaaaaaaaaaa,0xffffffffffffff,0xffffffffffffff,0xffffffffffffff,0xaaaaaaaaaaaaa9,0xaaaaaaaaaaaaa9,0xaaaaaaaaaaaaaa,0xaaaaaaaaaaaaaa};
//const spint constant_y[8]={0xe6e4977ea3136e,0x7349a5a276f460,0x756ed53c780c10,0x10ed8406b90c43,0x1fc34599fb4fd8,0xfb977d5a5a61d4,0x381ed8ab4fac45,0x9ea0dddcd95cca};
#endif
#if Radix == 28
const spint constant_x[16]={0x420685f,0x17aa383,0x35bf93b,0x4f8fe6d,0xb7bc291,0xab765fa,0xe44cd37,0xa7e9b28,0x34f39b1,0xbe886a8,0x44d6fb9,0x95e6eb4,0xee96c72,0x5d88b9,0xd16ef09,0xea8836};
const spint constant_y[16]={0xd81f4fb,0x9a99632,0xac119c7,0x23c2104,0xda8e9ac,0x9fc5486,0x416ef25,0x8a98abb,0x46ff58f,0x1223251,0x4fa9dd0,0xcaf9471,0xa1f0e6a,0xa33ab50,0x65f7687,0xa184177};
//const spint constant_x[16]={0xaaaaaaa,0xfffffff,0xfffffff,0xfffffff,0xfffffff,0xfffffff,0xfffffff,0xfffffff,0xaaaaaa9,0xaaaaaa9,0xaaaaaaa,0xaaaaaaa,0xaaaaaaa,0xaaaaaaa,0xaaaaaaa,0xaaaaaaa};
//const spint constant_y[16]={0xe6e4977,0x276f460,0x7349a5a,0xc780c10,0x756ed53,0x6b90c43,0x10ed840,0xb583c6a,0x1fc3458,0xa5a61d4,0xfb977d5,0xb4fac45,0x381ed8a,0xcd95cca,0x9ea0ddd,0xea3136e};
#endif
#define COF 2
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
    modneg(P->x,P->x);
}

// add Q to P
void ecnadd(point *Q,point *P)
{
    spint A[Nlimbs],B[Nlimbs],C[Nlimbs],D[Nlimbs],E[Nlimbs],F[Nlimbs],G[Nlimbs];
    modmul(Q->z,P->z,A);
    modsqr(A,B);
    modmul(Q->x,P->x,C);
    modmul(Q->y,P->y,D);
    modmul(C,D,E);
#ifdef CONSTANT_D
#if CONSTANT_D>0
    modmli(E,CONSTANT_D,E);
    modsub(B,E,F);
    modadd(B,E,G);
#else
    modmli(E,-CONSTANT_D,E);
    modadd(B,E,F);
    modsub(B,E,G);
#endif
#else
    modmul(E,constant_d,E);
    modsub(B,E,F);
    modadd(B,E,G);
#endif
    modadd(P->x,P->y,B);
    modadd(Q->x,Q->y,E);
    modmul(B,E,P->x);
    modsub(P->x,C,P->x);
    modsub(P->x,D,P->x);
    modmul(P->x,F,P->x);
    modmul(P->x,A,P->x);
#ifdef NEGA
    modadd(D,C,P->y);
#else
    modsub(D,C,P->y);
#endif
    modmul(P->y,A,P->y);
    modmul(P->y,G,P->y);
    modmul(F,G,P->z);
}

// subtract Q from P
void ecnsub(point *Q,point *P)
{
    point W;
    ecncpy(Q,&W); ecnneg(&W);
    ecnadd(&W,P);
}

// double P
void ecndbl(point *P)
{
    spint B[Nlimbs],C[Nlimbs],D[Nlimbs],E[Nlimbs],F[Nlimbs],H[Nlimbs],J[Nlimbs];
    modadd(P->x,P->y,B);
    modsqr(B,B);
    modsqr(P->x,C);
    modsqr(P->y,D);
    modsqr(P->z,H);
    modadd(H,H,H);
#ifdef NEGA
    modneg(C,E);
#else
    modcpy(C,E);
#endif
    modadd(E,D,F);
    modsub(F,H,J);
    modsub(B,C,P->x);
    modsub(P->x,D,P->x);
    modmul(P->x,J,P->x);
    modsub(E,D,P->y);
    modmul(P->y,F,P->y);
    modmul(F,J,P->z);
}
/*
      B = (X1+Y1)^2
      C = X1^2
      D = Y1^2
      H = 2Z1^2

      E = C                 
      F = E+D               F=C+D
      J = F-H               J=F-H
      X3 = (B-C-D)*J        X3=(B-F)*J
      Y3 = F*(E-D)          Y3=(C-D)*F
      Z3 = F*J              Z3=F*J

      E = -C
      F = E+D               F=-C+D          F=D-C           F=D-C
      J = F-H               J=F-H           J=-F-H          J=F+H
      X3 = (B-C-D)*J        X3=(B-C-D)*J    X3=(B-C-D)*J    X3=(C+D-B)*J
      Y3 = F*(E-D)          Y3=F*(-C-D)     Y3=F*(C+D)      Y3=F*(C+D)
      Z3 = F*J              Z3=F*J          Z3=-F*J         Z3=F*J


*/

// set to infinity
void ecninf(point *P)
{
    modzer(P->x);
    modone(P->y);
    modone(P->z);
}

// test for infinity
int ecnisinf(point *P)
{
    return (modis0(P->x) && modcmp(P->y,P->z));
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

// extract (x,y) from point, if y is NULL compress and just return x and sign of y, if x is NULL compress and just return y and sign of x
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

// general purpose set point function
// sets P=O if point not on curve
// if x and y are not NULL tries to set (x,y)
// if y==NULL tries to set from x and sign s of y (decompression)
// if x==NULL tries to set from y and sign s of x
static void setxy(int s,const spint *x,const spint *y,point *P)
{
    spint X[Nlimbs],Y[Nlimbs],O[Nlimbs],U[Nlimbs],V[Nlimbs],H[Nlimbs];
    modone(O);

    if (x!=NULL && y!=NULL)
    {
        modsqr(x,X);
        modsqr(y,Y);
#ifdef NEGA
        modsub(Y,X,U);
#else
        modadd(Y,X,U);  //lhs
#endif
        modmul(X,Y,V);  //rhs
#ifdef CONSTANT_D
#if CONSTANT_D>0
        modmli(V,CONSTANT_D,V);
        modadd(O,V,V); // V=1+dx^2
#else
        modmli(V,-CONSTANT_D,V);
        modsub(O,V,V); // V=1-dx^2
#endif      
#else
        modmul(V,constant_d,V);
        modadd(O,V,V);
#endif
        if (modcmp(U,V)) {
            modcpy(x,P->x);
            modcpy(y,P->y);
            modone(P->z);
        } else {
            ecninf(P);
        }
        return;
    }
    if (y==NULL)
    {
        modsqr(x,X);
#ifdef NEGA             // U=1-ax^2
        modadd(O,X,U);
#else
        modsub(O,X,U); 
#endif
        modcpy(X,V);   // V=x^2
    } else {
        modsqr(y,Y);
        modsub(O,Y,U); // U=1-y^2
#ifdef NEGA
        modneg(O,O);   // O=-1
#endif
        modcpy(Y,V);   // V=y^2
    }
#ifdef CONSTANT_D
#if CONSTANT_D>0
    modmli(V,CONSTANT_D,V);
    modsub(O,V,V); // V=1-dV^2
#else
    modmli(V,-CONSTANT_D,V);
    modadd(O,V,V); // V=1+dV^2
#endif
#else
    modmul(V,constant_d,V);
    modsub(O,V,V); // V=1-dV^2      
#endif

    modsqr(U,O);  // O=U^2
    modmul(U,O,U); // U=U^3
    modmul(U,V,U); // U=U^3*V
    modpro(U,H);
    if (!modqr(H,U))
    { // point not on curve
        ecninf(P);
        return;
    }
    modsqrt(U,H,V); // V=sqrt
    modinv(U,H,U);  // U=inv
    modmul(U,V,U);
    modmul(U,O,U);
    int d=(modsign(U)-s)&1;
    modneg(U,V);
    modcmv(d,V,U);
    if (y==NULL)
    {
        modcpy(U,P->y);
        modcpy(x,P->x);
    } else {
        modcpy(U,P->x);
        modcpy(y,P->y);
    }
    modone(P->z);
}

// multiply point by small curve cofactor (here assumed to be 4 or 8)
void ecncof(point *P)
{
    int i;
    for (i=0;i<COF;i++)
        ecndbl(P);
}

// Is (x,y) of the right order? Must be checked by calling program!
// api visible version, x and y are big endian byte arrays
void ecnset(int s,const char *x,const char *y,point *P)
{
    spint X[Nlimbs],Y[Nlimbs];
    if (x!=NULL && y!=NULL)
    {
        modimp(x,X);
        modimp(y,Y);
        setxy(s,X,Y,P);
    }
    if (x!=NULL)
    {
        modimp(x,X);
        setxy(s,X,NULL,P);
    }
    if (y!=NULL)
    {
        modimp(y,Y);
        setxy(s,NULL,Y,P);
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

//    printf("w= ");
//    for (i=0;i<8*Nbytes+8;i++) printf(" %d",(int)w[i]);
//    printf("\n");


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
