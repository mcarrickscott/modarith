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



/*** End of automatically generated code ***/

#define COUNT_CLOCKS
//#define USE_RDTSC

#ifdef COUNT_CLOCKS

#ifdef USE_RDTSC
#include <x86intrin.h>
#else
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
static void output(spint *x) {
    char b[Nbytes+1];
    char buff[(2*Nbytes)+1];
    modexp(x,b);
    toHex(b,buff);
    puts(buff);
}

// Describe Montgomery Curve parameters

#ifdef X25519
#define A24 121665  // Montgomery curve constant (A-2)/4
#define COF 3       // Montgomery curve cofactor = 2^cof (2 or 3)
#define GENERATOR 9
#define TWIST_SECURE // If it is a twist secure curve
#endif

#ifdef X448
#define A24 39081   // Montgomery curve constant (A-2)/4
#define COF 2       // Montgomery curve cofactor = 2^cof (2 or 3)
#define GENERATOR 5
#define TWIST_SECURE // If it is a twist secure curve
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

// RFC7748 - Montgomery curve
// bv=bk*bu, bu,bv are x coordinates on elliptic curve
void rfc7748(const char *bk, const char *bu,char *bv) {
    int i;
    int kt;
    int swap = 0;
    char ck[Nbytes];
    char cu[Nbytes];
    spint u[Nlimbs]; spint x1[Nlimbs]; spint x2[Nlimbs]; spint x3[Nlimbs]; spint z2[Nlimbs]; spint z3[Nlimbs];
    spint A[Nlimbs]; spint B[Nlimbs]; spint AA[Nlimbs]; spint BB[Nlimbs]; spint C[Nlimbs]; spint D[Nlimbs]; spint E[Nlimbs];

    for (i=0;i<Nbytes;i++) {
        ck[i]=bk[i];
        cu[i]=bu[i];
    }

    reverse(cu);  // convert from little to big endian
#ifdef X25519
    cu[0]&=0x7f;  // implementations of X25519 (but not X448) MUST mask the most significant bit in the final byte
#endif     
// clamp input
    clamp(ck);

// import into internal representation
    modimp(cu,u);

    modcpy(u,x1);  // x_1=u
    modone(x2);    // x_2=1
    modzer(z2);    // z_2=0
    modcpy(u,x3);  // x_3=u
    modone(z3);    // z_3=1

    for (i=Nbits-1;i>=0;i--)
    {
        kt=bit(i,ck);
        swap^=kt;
        modcsw(swap,x2,x3);
        modcsw(swap,z2,z3);
        swap=kt;
            
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

#ifdef TWIST_SECURE
    modpro(z2,A);       
    modinv(z2,A,z2);    // sufficient for twist secure curves like X25519 and X448 
#else
    // Do cheap point validation here - see https://eprint.iacr.org/2020/1497
    modmul(u,z2,B);     // wZ
    modmul(B,z2,A);     // wZ^2
    modpro(A,E);        // y
    modcpy(A,C);
    modmul(E,z2,D);     // y.Z2
    modsqr(D,D);     
    modmul(D,u,D);      // w.(y.z2)^2 
    for (i=0;i<(COF-2);i++) {  // COF is 2 or 3
        modsqr(C,C);
        modmul(C,A,C);
    }
    for (i=0;i<COF;i++) {
        modsqr(E,E);
    }
    modmul(C,E,C);
    modmul(C,B,z2);
    for (i=0;i<(COF-2);i++) {
        modsqr(D,D);
    }
    modone(A); modadd(D,A,D); modfsb(D); modshr(1,D); // 1 for QR, else 0
    modmul(x2,D,x2); // set to zero for bad input point
#endif
    modmul(x2,z2,x2);   

    modexp(x2,bv);
    reverse(bv); // convert to little endian
}

// a test vector for x25519 or x448 from RFC7748
int main()
{
#ifdef X25519
    const char *sk=(const char *)"77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a";
#endif
#ifdef X448
    const char *sk=(const char *)"9a8f4925d1519f5775cf46b04b5800d4ee9ee8bae8bc5565d498c28dd9c9baf574a9419744897391006382a6f127ab1d9ac2d8c0a598726b";
#endif
    uint64_t start,fin;
    clock_t begin;
    int i,elapsed;
    char sv[(Nbytes*2)+1];
    char bk[Nbytes],bv[Nbytes];
    char bu[Nbytes]={};

    bu[0]=GENERATOR;
// convert to byte array
    fromHex(sk,bk);

    rfc7748(bk,bu,bv);

// convert to Hex
    toHex(bv,sv);

    puts(sk); 
    puts(sv); 

#ifdef COUNT_CLOCKS
#ifdef USE_RDTSC
    start=__rdtsc();
#else
    start=cpucycles();
#endif    
#endif
    begin=clock();
    for (i=0;i<5000;i++) {
        rfc7748(bk,bu,bv);
        rfc7748(bk,bv,bu);
    }
    elapsed=100*(clock() - begin) / CLOCKS_PER_SEC;
#ifdef COUNT_CLOCKS
#ifdef USE_RDTSC
    fin=__rdtsc();
#else
    fin=cpucycles();
#endif
    printf("Clock cycles= %d\n",(int)((fin-start)/10000ULL));
#endif

    printf("Microseconds= %d\n",elapsed);
    toHex(bu,sv);
    puts(sv);
}
