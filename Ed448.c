// Ed448 Implementation
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf
// and RFC8032
// python curve.py 64 ED448
// This completes edwards.rs for this curve. Then
// gcc -O2 Ed448.c edwards.c hash.c -o Ed448

#include <stdio.h>
#include <stdint.h>
#include "hash.h"  // Some useful hash functions

#include "curve.h"   // elliptic curve API

/*** Insert code automatically generated from group.c here ***/
/* Note that much of this code is not needed and can be deleted */



/*** End of automatically generated code ***/

// number of limbs and bytes in representation
#define BYTES Nbytes
#define LIMBS Nlimbs

// Some utility functions for I/O and debugging

// general purpose SHAKE256 hash function
// Input ilen bytes, output olen bytes
static void H(int ilen,int olen,char *s,char *digest)
{
    sha3 SHA3;
    SHA3_init(&SHA3,SHAKE256);
    for (int i=0;i<ilen;i++) 
        SHA3_process(&SHA3,s[i]);
    SHA3_shake(&SHA3,digest,olen); 
}

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
    toHex(56,b,buff);
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

// reduce 114 byte array h to integer r modulo group order q, in constant time
// Consider h as 2^472.x + 2^440.y + z, where x,y and z < q (z is first 55 bytes, y is next 4 bytes, x is last 55 bytes)
// precalculate c1=nres(2^472 mod q) and c2=nres(2^440 mod q)
// using utility ed448_order.py 
#if Wordlength==64
const spint constant_c1[8]={0xe3033c23525654,0x7624da8d5b86ce,0x3a503352aa569,0x3337c35f209580,0xae17cf72c9860f,0x9cc14ba3c47c44,0xbcb7e4d070af1a,0x39f823b7292052};
const spint constant_c2[8]={0xecfdd2acf1d7f0,0x1d5087ce8f0058,0xbe9c8459d676ba,0xff8fc444c3d266,0xa3c47c44ae17ce,0xd070af1a9cc14b,0xb7292052bcb7e4,0x383402a939f823};
#endif

#if Wordlength==32
const spint constant_c1[16]={0x3525654,0xe3033c2,0xd5b86ce,0x7624da8,0x52aa569,0x3a5033,0xf209580,0x3337c35,0x2c9860f,0xae17cf7,0x3c47c44,0x9cc14ba,0x70af1a,0xbcb7e4d,0x7292052,0x39f823b};
const spint constant_c2[16]={0xcf1d7f0,0xecfdd2a,0xe8f0058,0x1d5087c,0x9d676ba,0xbe9c845,0x4c3d266,0xff8fc44,0x4ae17ce,0xa3c47c4,0xa9cc14b,0xd070af1,0x2bcb7e4,0xb729205,0x939f823,0x383402a};
#endif

static void reduce(char *h,spint *r)
{
    int i;
    char buff[BYTES];
    spint x[LIMBS],y[LIMBS],z[LIMBS];

    for (int i=0;i<55;i++)
        buff[i]=h[i];
    buff[55]=0;
    reverse(buff);
    modimp(buff,z);

    for (i=0;i<4;i++ )
        buff[i]=h[55+i];
    for (i=4;i<56;i++)
        buff[i]=0;
    reverse(buff);
    modimp(buff,y);

    for (i=0;i<55;i++)
        buff[i]=h[59+i];
    buff[55]=0;
    reverse(buff);
    modimp(buff,x);

    modmul(x,constant_c1,x);
    modmul(y,constant_c2,y);
    modadd(x,y,r);
    modadd(r,z,r);
}

// Input private key - 57 random bytes
// Output public key - 57 bytes
void ED448_KEY_GEN(char *prv,char *pub)
{
    int sign;
    point P;
    ecngen(&P);
    char s[BYTES];

    H(BYTES+1,BYTES,prv,s);
// clamp s
    s[0]&=0xFC;
    s[55]|=0x80;

    reverse(s);  // little endian to big endian
    ecnmul(s,&P); 

    sign=ecnget(&P,NULL,pub);  // get y coordinate and sign
    reverse(pub);              // big endian to little endian
    pub[56]=(char)(sign<<7);
}

const char dom4[10]={'S','i','g','E','d','4','4','8',0,0};

// input private key, public key, message to be signed. Output signature
void ED448_SIGN(char *prv,char *pub,char *m,char *sig)
{
    int i,sign;
    char h[2*BYTES+2];
    spint r[LIMBS],s[LIMBS],d[LIMBS];
    sha3 SHA3;
    point R;
    ecngen(&R);

    SHA3_init(&SHA3,SHAKE256);
    H(BYTES+1,2*BYTES+2,prv,h);

// derive and clamp s
    h[0]&=0xFC;
    h[55]|=0x80;
    reverse(h); 
    modimp(h,s);

    for (i=0;i<10;i++)
        SHA3_process(&SHA3,dom4[i]);
    for (i=BYTES+1;i<2*BYTES+2;i++ )
        SHA3_process(&SHA3,h[i]);
    i=0;
    while (m[i]!=0)
        SHA3_process(&SHA3,m[i++]);
    SHA3_shake(&SHA3,h,2*BYTES+2); 

    reduce(h,r);
    modexp(r,h);  // convert to big endian array
    ecnmul(h,&R);

    sign=ecnget(&R,NULL,sig);  // get y coordinate and sign
    reverse(sig);              // big endian to little endian
    sig[BYTES]=(char)(sign<<7); // first part of signature

    SHA3_init(&SHA3,SHAKE256);
    for (i=0;i<10;i++)
        SHA3_process(&SHA3,dom4[i]);
    for (i=0;i<BYTES+1;i++ )
        SHA3_process(&SHA3,sig[i]);  // R
    for (int i=0;i<BYTES+1;i++)
        SHA3_process(&SHA3,pub[i]);  // Q
    i=0;
    while (m[i]!=0)
        SHA3_process(&SHA3,m[i++]);  // M
    SHA3_shake(&SHA3,h,2*BYTES+2); 

    reduce(h,d);
    modmul(d,s,d);
    modadd(d,r,d);

    modexp(d,&sig[BYTES+1]);
    reverse(&sig[BYTES+1]);
    sig[2*BYTES+1]=0;           // second part of signature
}

int ED448_VERIFY(char *pub,char *m,char *sig) 
{
    int i;
    char buff[BYTES];
    point G,R,Q;
    spint u[LIMBS];
    int sign;
    sha3 SHA3;
    char h[2*BYTES+2];

    ecngen(&G);

// reconstruct point R
    for (i=0;i<BYTES;i++)
        buff[i]=sig[i];
    reverse(buff);
    sign=sig[BYTES]>>7;
    ecnset(sign,NULL,buff,&R);

// reconstruct point Q 
    for (i=0;i<BYTES;i++)
        buff[i]=pub[i];
    reverse(buff);
    sign=(pub[BYTES]>>7)&1;
    ecnset(sign,NULL,buff,&Q);

    for (i=0;i<BYTES;i++)
        buff[i]=sig[i+BYTES+1];
    reverse(buff);

    SHA3_init(&SHA3,SHAKE256);
    for (i=0;i<10;i++)
        SHA3_process(&SHA3,dom4[i]);
    for (i=0;i<BYTES+1;i++ )
        SHA3_process(&SHA3,sig[i]);  // R
    for (int i=0;i<BYTES+1;i++)
        SHA3_process(&SHA3,pub[i]);  // Q
    i=0;
    while (m[i]!=0)
        SHA3_process(&SHA3,m[i++]);  // M
    SHA3_shake(&SHA3,h,2*BYTES+2); 

    reduce(h,u); modneg(u,u); modexp(u,h);

    if (!modimp(buff,u)) return 0;  // out of range

    ecncof(&G); ecncof(&R); ecncof(&Q);
    ecnmul2(buff,&G,h,&Q,&Q);

    if (ecncmp(&R,&Q))
        return 1;
    return 0;
}

// Test vector from RFC8032

int main()
{
    const char *sk=(const char *)"c4eab05d357007c632f3dbb48489924d552b08fe0c353a0d4a1f00acda2c463afbea67c5e8d2877c5e3bc397a659949ef8021e954e0a12274e";
    char prv[BYTES+1],pub[BYTES+1],sig[2*BYTES+2];
    char buff[256],m[2];
    int res;
    printf("Run RFC8032 test vector\n");
    printf("private key= "); puts(sk); 
    fromHex(BYTES+1,sk,prv);
    ED448_KEY_GEN(prv,pub);
    toHex(BYTES+1,pub,buff);
    printf("public key=  "); puts(buff);

    m[0]=0x03; // message to be signed
    m[1]=0;    // null terminated string
    ED448_SIGN(prv,pub,m,sig);
    toHex(2*BYTES+2,sig,buff);
    printf("signature=  "); puts(buff); 

    res=ED448_VERIFY(pub,m,sig);
    if (res)
        printf("Signature is valid\n");
    else
        printf("Signature is NOT valid\n");
}
