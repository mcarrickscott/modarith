
// ECDSA Implementation for curve P-256
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf
// python curve.py 64 NIST256
// This completes weierstrass.c for this curve. Then
// gcc -O2 EC256.c weierstrass.c hash.c -o EC256

#include <stdio.h>
#include <stdint.h>
#include "hash.h"  // Some useful hash functions

#include "curve.h"   // elliptic curve API

/*** Insert automatically generated code for P-256 prime group order group.c here ***/
/* Note that much of this code is not needed and can be deleted */



/*** End of automatically generated code ***/

// number of bytes in representation
#define BYTES Nbytes
typedef spint gel[Nlimbs];  // group element definition

// Some utility functions for I/O and debugging

// reverse bytes of buff - for little endian
static void reverse(char *buff) {
    int n=BYTES;
    for (int i = 0; i < n/2; i++) { 
        char ch = buff[i]; 
        buff[i] = buff[n - i - 1]; 
        buff[n - i - 1] = ch; 
    } 
}

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

// Convert from a hex string to byte array 
static void fromHex(int ilen, const char *src, char *dst)
{
    int i,lz,len=0;
    char pad[128];
    while (src[len]!=0) len++;
    lz=2*ilen-len;
    if (lz<0) lz=0;
    for (i=0;i<lz;i++) pad[i]='0';  // pad with leading zeros
    for (i=lz;i<2*ilen;i++) pad[i]=src[i-lz];

    for (i=0;i<ilen;i++)
    {
        dst[i] = (char2int(pad[2*i]) * 16) + char2int(pad[2*i + 1]);
    }
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
static void toHex(int len, const char *src, char *dst)
{
    int i;
    for (i = 0; i < len; i++)
    {
        unsigned char ch = src[i];
        byte2hex(&dst[i * 2],ch);
    }
    dst[2*len]='\0';
}

// I/O debug code
// output a modulo number in hex
/*
static void output(spint *x) {
    char b[Nbytes+1];
    char buff[(2*Nbytes)+1];
    modexp(x,b);
    toHex(Nbytes,b,buff);
    puts(buff);
}

// output a point (x,y)
void outputxy(point *P)
{
    if (ecnisinf(P)) {
        printf("P= O\n");
    } else {
        char x[BYTES],y[BYTES];
        char buff[(2*BYTES)+1];
        ecnget(P,x,y);
        toHex(BYTES,x,buff);
        printf("Px= "); puts(buff);
        toHex(BYTES,y,buff);
        printf("Py= "); puts(buff);
    }
}
*/

// reduce 40 byte array h to integer r modulo group order q, in constant time
// Consider h as 2^248.x + y, where x and y < q (x is top 9 bytes, y is bottom 31 bytes)
// Important that x and y < q
// precalculate c=nres(2^248 mod q) - see ec256_order.py
#if Wordlength==64
const spint constant_c[5]={0x57b7a0c28d8f5,0x31e8132cb7905,0x92b6bec16b3c6,0xd9571e2845b23,0xfe66e12c96f3};
#endif

#if Wordlength==32
const spint constant_c[9]={0xc28d8f5,0x2abdbd0,0x4cb2de4,0x78c63d0,0x16bec16b,0x2d91c95,0x55c78a1,0x592de7b,0xfe66e1};
#endif

static void reduce(char *h,spint *r)
{
    int i;
    char buff[BYTES];
    gel x,y,z;

    for (i=0;i<31;i++)
        buff[i]=h[i];  // little endian
    buff[31]=0;
    reverse(buff);
    modimp(buff,y);

    for (i=0;i<9;i++)
        buff[i]=h[31+i];
    for (i=9;i<32;i++)
        buff[i]=0;
    reverse(buff);
    modimp(buff,x);

    modmul(x,constant_c,x);  // 2^248.x 
    modadd(x,y,r);
}

//#define COMPRESS 
#define PREHASHED   // define for test vectors

// Input private key - 32 random bytes
// Output public key - 65 bytes (0x04<x>|<y>), or 33 if compressed (0x02<x>.. or 0x03<x>)
void EC256_KEY_GEN(char *prv,char *pub)
{
    point P;
    ecngen(&P);

    ecnmul(prv,&P); 

#ifdef COMPRESS
    pub[0]=0x02+ecnget(&P,&pub[1],NULL); // 0x02 or 0x03
#else
    pub[0]=0x04; // no compression
    ecnget(&P,&pub[1],&pub[BYTES+1]);  // get x and y
#endif
}

// input private key, per-message random number, message to be signed. Output signature.
// ran must be Nbytes+8 in length, in this case 40 bytes
void EC256_SIGN(char *prv,char *ran,char *m,char *sig)
{
    char h[BYTES];
    point R;
    gel e,r,s,k;

#ifdef PREHASHED
    modimp(m,e);
#else
    int i;
    char h[BYTES];
    hash256 sh256;
    HASH256_init(&sh256);
    i=0;
    while (m[i]!=0)
        HASH256_process(&sh256,m[i++]);
    HASH256_hash(&sh256,h); 

    modimp(h,e);
#endif

    ecngen(&R);
    modimp(prv,s);

    reduce(ran,k);
    modexp(k,h);
    ecnmul(h,&R);
    modinv(k,NULL,k);

    //ecnmul(ran,&R);
    //modimp(ran,k);
    //modinv(k,NULL,k);
    ecnget(&R,h,NULL);
    modimp(h,r);

    modmul(s,r,s);
    modadd(s,e,s);
    modmul(s,k,s);
    modzer(k);

    modexp(r,sig);
    modexp(s,&sig[BYTES]);
}

// input public key, message and signature
int EC256_VERIFY(char *pub,char *m,char *sig) 
{
    point G,Q;
    int i,res;
    char rb[BYTES],u[BYTES],v[BYTES];
    gel e,r,s;
#ifdef PREHASHED
    modimp(m,e);
#else
    char h[BYTES];
    hash256 sh256;
    HASH256_init(&sh256);
    i=0;
    while (m[i]!=0)
        HASH256_process(&sh256,m[i++]);
    HASH256_hash(&sh256,h); 

    modimp(h,e);
#endif

    ecngen(&G);

// import from signature
    if (!modimp(sig,r)) return 0; // if not in range
    if (!modimp(&sig[BYTES],s)) return 0;

    if (modis0(r) || modis0(s)) return 0;
    modinv(s,NULL,s);
    modmul(r,s,r); modexp(r,v);  // export to byte array
    modmul(s,e,s); modexp(s,u); 

    if (pub[0]==0x04) {
        ecnset(0,&pub[1],&pub[BYTES+1],&Q);
    } else {
        ecnset((int)pub[0]&1,&pub[1],NULL,&Q);
    }

    ecnmul2(u,&G,v,&Q,&Q);
    if (ecnisinf(&Q)) return 0;

    ecnget(&Q,rb,NULL);

    res=1;
    for (i=0;i<BYTES;i++) {
        if (sig[i]!=rb[i]) res=0;
    }
    
    return res;
}

// test for FIPS 186-3 ECDSA Signature Generation

int main()
{
    const char *sk= (const char *)"519b423d715f8b581f4fa8ee59f4771a5b44c8130b4e3eacca54a56dda72b464";
    const char *ran=(const char *)"94a1bbb14b906a61a280f245f9e93c7f3b4a6247824f5d33b9670787642a68deb9670787642a68de";
    const char *msg=(const char *)"44acf6b7e36c1342c2c5897204fe09504e1e2efb1a900377dbc4e7a6a133ec56";
    char prv[BYTES],pub[2*BYTES+1];
    char buff[256],m[BYTES],k[BYTES+8],sig[2*BYTES];
    int res;
    printf("Run test vector\n");
    printf("private key= "); puts(sk); 
    fromHex(BYTES,sk,prv);
    fromHex(BYTES+8,ran,k);
    fromHex(BYTES,msg,m);
    EC256_KEY_GEN(prv,pub);
#ifdef COMPRESS
    toHex(BYTES+1,pub,buff);
#else
    toHex(2*BYTES+1,pub,buff);
#endif
    printf("public key= "); puts(buff);
    EC256_SIGN(prv,k,m,sig);
    toHex(2*BYTES,sig,buff);
    printf("signature=  "); puts(buff);

    res=EC256_VERIFY(pub,m,sig);
    if (res)
        printf("Signature is valid\n");
    else
        printf("Signature is NOT valid\n");
}
