// Ed448 Implementation
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf
// and RFC8032

#include <stdio.h>
#include <stdint.h>

/* SHA3 */

/**
 * @brief SHA3 hash function instance */
typedef struct
{
    int length;   /**< 64-bit input length */
    uint64_t S[25];  /**< Internal state */
    int rate;          /**< TODO */
    int len;           /**< Hash length in bytes */
} sha3;

#define SHA3_HASH224 28 /**< SHA3 224 bit hash */
#define SHA3_HASH256 32 /**< SHA3 256 bit hash */
#define SHA3_HASH384 48 /**< SHA3 384 bit hash */
#define SHA3_HASH512 64 /**< SHA3 512 bit hash */

#define SHAKE128 16 /**< SHAKE128   hash */
#define SHAKE256 32 /**< SHAKE256 hash */

#define SHA3_ROUNDS 24
#define rotl(x,n) (((x)<<n) | ((x)>>(64-n)))

/* round constants */

static const uint64_t RC[24] =
{
    0x0000000000000001UL, 0x0000000000008082UL, 0x800000000000808AUL, 0x8000000080008000UL,
    0x000000000000808BUL, 0x0000000080000001UL, 0x8000000080008081UL, 0x8000000000008009UL,
    0x000000000000008AUL, 0x0000000000000088UL, 0x0000000080008009UL, 0x000000008000000AUL,
    0x000000008000808BUL, 0x800000000000008BUL, 0x8000000000008089UL, 0x8000000000008003UL,
    0x8000000000008002UL, 0x8000000000000080UL, 0x000000000000800AUL, 0x800000008000000AUL,
    0x8000000080008081UL, 0x8000000000008080UL, 0x0000000080000001UL, 0x8000000080008008UL
};

/* permutation */

static void SHA3_transform(sha3 *sh)
{
    int k;
    uint64_t B00,B01,B02,B03,B04,B10,B11,B12,B13,B14,B20,B21,B22,B23,B24,B30,B31,B32,B33,B34,B40,B41,B42,B43,B44;
    uint64_t C0,C1,C2,C3,C4,D0,D1,D2,D3,D4;

    for (k = 0; k < SHA3_ROUNDS; k++)
    {

        C0=sh->S[0] ^ sh->S[5] ^ sh->S[10] ^ sh->S[15] ^ sh->S[20];
        C1=sh->S[1] ^ sh->S[6] ^ sh->S[11] ^ sh->S[16] ^ sh->S[21];
        C2=sh->S[2] ^ sh->S[7] ^ sh->S[12] ^ sh->S[17] ^ sh->S[22];
        C3=sh->S[3] ^ sh->S[8] ^ sh->S[13] ^ sh->S[18] ^ sh->S[23];
        C4=sh->S[4] ^ sh->S[9] ^ sh->S[14] ^ sh->S[19] ^ sh->S[24];

        D0 = C4 ^ rotl(C1, 1);
        D1 = C0 ^ rotl(C2, 1);
        D2 = C1 ^ rotl(C3, 1);
        D3 = C2 ^ rotl(C4, 1);
        D4 = C3 ^ rotl(C0, 1);

        B00 =      sh->S[0]^D0;
        B10 = rotl(sh->S[6]^D1, 44);
        B20 = rotl(sh->S[12]^D2, 43);
        B30 = rotl(sh->S[18]^D3, 21);
        B40 = rotl(sh->S[24]^D4, 14);    

        B01 = rotl(sh->S[3]^D3, 28);
        B11 = rotl(sh->S[9]^D4, 20);
        B21 = rotl(sh->S[10]^D0, 3);
        B31 = rotl(sh->S[16]^D1, 45);
        B41 = rotl(sh->S[22]^D2, 61);

        B02 = rotl(sh->S[1]^D1, 1);
        B12 = rotl(sh->S[7]^D2, 6);
        B22 = rotl(sh->S[13]^D3, 25);
        B32 = rotl(sh->S[19]^D4, 8);
        B42 = rotl(sh->S[20]^D0, 18);

        B03 = rotl(sh->S[4]^D4, 27);
        B13 = rotl(sh->S[5]^D0, 36);
        B23 = rotl(sh->S[11]^D1, 10);
        B33 = rotl(sh->S[17]^D2, 15);
        B43 = rotl(sh->S[23]^D3, 56);

        B04 = rotl(sh->S[2]^D2, 62);
        B14 = rotl(sh->S[8]^D3, 55);
        B24 = rotl(sh->S[14]^D4, 39);
        B34 = rotl(sh->S[15]^D0, 41);
        B44 = rotl(sh->S[21]^D1, 2);

        sh->S[0]=B00^(~B10&B20);
        sh->S[1]=B10^(~B20&B30);
        sh->S[2]=B20^(~B30&B40);
        sh->S[3]=B30^(~B40&B00);
        sh->S[4]=B40^(~B00&B10);

        sh->S[5]=B01^(~B11&B21);
        sh->S[6]=B11^(~B21&B31);
        sh->S[7]=B21^(~B31&B41);
        sh->S[8]=B31^(~B41&B01);
        sh->S[9]=B41^(~B01&B11);

        sh->S[10]=B02^(~B12&B22);
        sh->S[11]=B12^(~B22&B32);
        sh->S[12]=B22^(~B32&B42);
        sh->S[13]=B32^(~B42&B02);
        sh->S[14]=B42^(~B02&B12);

        sh->S[15]=B03^(~B13&B23);
        sh->S[16]=B13^(~B23&B33);
        sh->S[17]=B23^(~B33&B43);
        sh->S[18]=B33^(~B43&B03);
        sh->S[19]=B43^(~B03&B13);

        sh->S[20]=B04^(~B14&B24);
        sh->S[21]=B14^(~B24&B34);
        sh->S[22]=B24^(~B34&B44);
        sh->S[23]=B34^(~B44&B04);
        sh->S[24]=B44^(~B04&B14);

        sh->S[0] ^= RC[k];
    }
}

/* Re-Initialize. olen is output length in bytes -
   should be 28, 32, 48 or 64 (224, 256, 384, 512 bits resp.) */

void SHA3_init(sha3 *sh, int olen)
{
    int i;
    for (i = 0; i < 25; i++)
        sh->S[i] = 0;  /* 5x5x8 bytes = 200 bytes of state */
    sh->length = 0;
    sh->len = olen;
    sh->rate = 200 - 2 * olen; /* number of bytes consumed in one gulp. Note that some bytes in the
                            state ("capacity") are not touched. Gulps are smaller for larger digests.
                            Important that olen<rate */
}

/* process a single byte */
void SHA3_process(sha3 *sh, int byt)
{
    int cnt = (int)(sh->length);
    int b = cnt % 8;
    cnt /= 8;
    sh->S[cnt] ^= ((uint64_t)(byt&0xff) << (8 * b));
    sh->length++;
    if (sh->length == sh->rate) {
        sh->length=0;
        SHA3_transform(sh);
    }
}

/* squeeze the sponge */
void SHA3_squeeze(sha3 *sh, char *buff, int len)
{
    int i, j, k, m = 0;
    uint64_t el;
    int nb=len/sh->rate;

    for (j=0;j<nb;j++ )
    {
        for (i=0;i<sh->rate/8;i++)
        {
            el=sh->S[i];
            for (k=0;k<8;k++)
            {
               buff[m++] = (el & 0xff);
               el >>= 8;
            }
        }
        SHA3_transform(sh);    
    }
   
    i=0;
    while (m<len)
    {
        el = sh->S[i++];
        for (k = 0; k < 8; k++)
        {
            buff[m++] = (el & 0xff);
            if (m >= len) break;
            el >>= 8;
        }    
    } 
}

void SHA3_hash(sha3 *sh, char *hash)
{
    /* generate a SHA3 hash of appropriate size */
    int q = sh->rate - sh->length;
    if (q == 1) SHA3_process(sh, 0x86);
    else
    {
        SHA3_process(sh, 0x06);  /* 0x06 for SHA-3 */
        while ((int)sh->length != sh->rate - 1) SHA3_process(sh, 0x00);
        SHA3_process(sh, 0x80); /* this will force a final transform */
    }
    SHA3_squeeze(sh, hash, sh->len);
}

/* return intermediate hash */
void SHA3_continuing_hash(sha3 *sh,char *digest)
{
    sha3 cp=*sh;
    SHA3_hash(&cp,digest);
}

void SHA3_shake(sha3 *sh, char *buff, int len)
{
    /* SHAKE out a buffer of variable length len */
    int q = sh->rate - sh->length;
    if (q == 1) SHA3_process(sh, 0x9f);
    else
    {
        SHA3_process(sh, 0x1f);  // 0x06 for SHA-3 !!!!
        while ((int) sh->length != sh->rate - 1) SHA3_process(sh, 0x00);
        SHA3_process(sh, 0x80); /* this will force a final transform */
    }
    SHA3_squeeze(sh, buff, len);
}

/* return intermediate hash */
void SHA3_continuing_shake(sha3 *sh,char *digest,int len)
{
    sha3 cp=*sh;
    SHA3_shake(&cp,digest,len);
}

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

#include "curve.h"   // elliptic curve API

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


/*** Insert automatically generated code for Ed448 prime group order code.c here ***/

// python monty.py 64 0x3fffffffffffffffffffffffffffffffffffffffffffffffffffffff7cca23e9c44edb49aed63690216cc2728dc58f552378c292ab5844f3

/* Note that much of this code is not needed and can be deleted */



/*** End of automatically generated code ***/

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
    s[55]!=0x80;

    reverse(s);  // little endian to big endian
    ecnmul(s,&P);

    sign=ecnget(&P,NULL,pub);  // get y coordinate and sign
    reverse(pub);              // big endian to little endian
    pub[56]=(char)(sign<<7);
}

const char dom4[10]={'S','i','g','E','d','4','4','8',0,0};

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
    h[55]!=0x80;
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
    spint t[LIMBS],u[LIMBS];
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
