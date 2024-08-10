
// ECDSA Implementation for curve P-256
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf

#include <stdio.h>
#include <stdint.h>

// SHA256

/**
 * @brief SHA256 hash function instance */
typedef struct
{
    uint32_t length[2]; /**< 64-bit input length */
    uint32_t h[8];      /**< Internal state */
    uint32_t w[64];	/**< Internal state */
    int hlen;		/**< Hash length in bytes */
} hash256;

#define H0_256 0x6A09E667L
#define H1_256 0xBB67AE85L
#define H2_256 0x3C6EF372L
#define H3_256 0xA54FF53AL
#define H4_256 0x510E527FL
#define H5_256 0x9B05688CL
#define H6_256 0x1F83D9ABL
#define H7_256 0x5BE0CD19L

static const uint32_t K_256[64] =
{
    0x428a2f98L, 0x71374491L, 0xb5c0fbcfL, 0xe9b5dba5L, 0x3956c25bL, 0x59f111f1L, 0x923f82a4L, 0xab1c5ed5L,
    0xd807aa98L, 0x12835b01L, 0x243185beL, 0x550c7dc3L, 0x72be5d74L, 0x80deb1feL, 0x9bdc06a7L, 0xc19bf174L,
    0xe49b69c1L, 0xefbe4786L, 0x0fc19dc6L, 0x240ca1ccL, 0x2de92c6fL, 0x4a7484aaL, 0x5cb0a9dcL, 0x76f988daL,
    0x983e5152L, 0xa831c66dL, 0xb00327c8L, 0xbf597fc7L, 0xc6e00bf3L, 0xd5a79147L, 0x06ca6351L, 0x14292967L,
    0x27b70a85L, 0x2e1b2138L, 0x4d2c6dfcL, 0x53380d13L, 0x650a7354L, 0x766a0abbL, 0x81c2c92eL, 0x92722c85L,
    0xa2bfe8a1L, 0xa81a664bL, 0xc24b8b70L, 0xc76c51a3L, 0xd192e819L, 0xd6990624L, 0xf40e3585L, 0x106aa070L,
    0x19a4c116L, 0x1e376c08L, 0x2748774cL, 0x34b0bcb5L, 0x391c0cb3L, 0x4ed8aa4aL, 0x5b9cca4fL, 0x682e6ff3L,
    0x748f82eeL, 0x78a5636fL, 0x84c87814L, 0x8cc70208L, 0x90befffaL, 0xa4506cebL, 0xbef9a3f7L, 0xc67178f2L
};

#define PAD  0x80
#define ZERO 0

/* functions */

#define S(m,n,x) (((x)>>n) | ((x)<<(m-n)))
#define R(n,x) ((x)>>n)

#define Ch(x,y,z)  ((x&y)^(~(x)&z))
#define Maj(x,y,z) ((x&y)^(x&z)^(y&z))
#define Sig0_256(x)    (S(32,2,x)^S(32,13,x)^S(32,22,x))
#define Sig1_256(x)    (S(32,6,x)^S(32,11,x)^S(32,25,x))
#define theta0_256(x)  (S(32,7,x)^S(32,18,x)^R(3,x))
#define theta1_256(x)  (S(32,17,x)^S(32,19,x)^R(10,x))

#define Sig0_512(x)    (S(64,28,x)^S(64,34,x)^S(64,39,x))
#define Sig1_512(x)    (S(64,14,x)^S(64,18,x)^S(64,41,x))
#define theta0_512(x)  (S(64,1,x)^S(64,8,x)^R(7,x))
#define theta1_512(x)  (S(64,19,x)^S(64,61,x)^R(6,x))


/* SU= 72 */
static void HASH256_transform(hash256 *sh)
{
    /* basic transformation step */
    uint32_t a, b, c, d, e, f, g, h, t1, t2;
    int j;
    for (j = 16; j < 64; j++)
        sh->w[j] = theta1_256(sh->w[j - 2]) + sh->w[j - 7] + theta0_256(sh->w[j - 15]) + sh->w[j - 16];

    a = sh->h[0];
    b = sh->h[1];
    c = sh->h[2];
    d = sh->h[3];
    e = sh->h[4];
    f = sh->h[5];
    g = sh->h[6];
    h = sh->h[7];

    for (j = 0; j < 64; j++)
    {
        /* 64 times - mush it up */
        t1 = h + Sig1_256(e) + Ch(e, f, g) + K_256[j] + sh->w[j];
        t2 = Sig0_256(a) + Maj(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }

    sh->h[0] += a;
    sh->h[1] += b;
    sh->h[2] += c;
    sh->h[3] += d;
    sh->h[4] += e;
    sh->h[5] += f;
    sh->h[6] += g;
    sh->h[7] += h;
}

/* Initialise Hash function */
void HASH256_init(hash256 *sh)
{
    /* re-initialise */
    int i;
    for (i = 0; i < 64; i++) sh->w[i] = 0L;
    sh->length[0] = sh->length[1] = 0L;
    sh->h[0] = H0_256;
    sh->h[1] = H1_256;
    sh->h[2] = H2_256;
    sh->h[3] = H3_256;
    sh->h[4] = H4_256;
    sh->h[5] = H5_256;
    sh->h[6] = H6_256;
    sh->h[7] = H7_256;

    sh->hlen = 32;
}

/* process a single byte */
void HASH256_process(hash256 *sh, int byt)
{
    /* process the next message byte */
    int cnt;
    cnt = (int)((sh->length[0] / 32) % 16);

    sh->w[cnt] <<= 8;
    sh->w[cnt] |= (uint32_t)(byt & 0xFF);

    sh->length[0] += 8;
    if (sh->length[0] == 0L)
    {
        sh->length[1]++;
        sh->length[0] = 0L;
    }
    if ((sh->length[0] % 512) == 0) HASH256_transform(sh);
}

/* SU= 24 */
/* Generate 32-byte Hash */
void HASH256_hash(hash256 *sh, char *digest)
{
    /* pad message and finish - supply digest */
    int i;
    uint32_t len0, len1;
    len0 = sh->length[0];
    len1 = sh->length[1];
    HASH256_process(sh, PAD);
    while ((sh->length[0] % 512) != 448) HASH256_process(sh, ZERO);
    sh->w[14] = len1;
    sh->w[15] = len0;
    HASH256_transform(sh);
    for (i = 0; i < sh->hlen; i++)
    {
        /* convert to bytes */
        digest[i] = (char)((sh->h[i / 4] >> (8 * (3 - i % 4))) & 0xffL);
    }
    HASH256_init(sh);
}

#include "curve.h"   // elliptic curve API

/*** Insert automatically generated code for P-256 prime group order code.c here ***/

// python monty.py 64 nist256

/* Note that much of this code is not needed and can be deleted */



/*** End of automatically generated code ***/

// number of limbs and bytes in representation
#define BYTES Nbytes
#define LIMBS Nlimbs

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

//#define COMPRESS 
#define PREHASHED   // define for test vectors

// Input private key - 32 random bytes
// Output public key - 65 bytes, or 33 is compressed
void EC256_KEY_GEN(char *prv,char *pub)
{
    point P;
    ecngen(&P);

    ecnmul(prv,&P); 

#ifdef COMPRESS
    pub[0]=0x02+ecnget(&P,&pub[1],NULL);
#else
    pub[0]=0x04; // no compression
    ecnget(&P,&pub[1],&pub[BYTES+1]);  // get x and y
#endif
}

// input private key, per-message random number, message to be signed. Output signature.
void EC256_SIGN(char *prv,char *ran,char *m,char *sig)
{
    char rb[BYTES];
    point R;
    spint e[LIMBS],r[LIMBS],s[LIMBS],k[LIMBS];

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
    ecnmul(ran,&R);
    modimp(ran,k);
    modinv(k,NULL,k);
    ecnget(&R,rb,NULL);
    modimp(rb,r);

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
    spint e[LIMBS],r[LIMBS],s[LIMBS];
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

    if (!modimp(sig,r)) return 0; // if not in range
    if (!modimp(&sig[BYTES],s)) return 0;

    if (modis0(r) || modis0(s)) return 0;
    modinv(s,NULL,s);
    modmul(r,s,r); modexp(r,v); 
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

int main()
{
    const char *sk= (const char *)"519b423d715f8b581f4fa8ee59f4771a5b44c8130b4e3eacca54a56dda72b464";
    const char *ran=(const char *)"94a1bbb14b906a61a280f245f9e93c7f3b4a6247824f5d33b9670787642a68de";
    const char *msg=(const char *)"44acf6b7e36c1342c2c5897204fe09504e1e2efb1a900377dbc4e7a6a133ec56";
    char prv[BYTES],pub[2*BYTES+1];
    char buff[256],m[BYTES],k[BYTES],sig[2*BYTES];
    int res;
    printf("Run test vector\n");
    printf("private key= "); puts(sk); 
    fromHex(BYTES,sk,prv);
    fromHex(BYTES,ran,k);
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
