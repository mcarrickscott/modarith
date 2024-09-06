// Some useful hash functions, SHA2 and SHA3
#include "hash.h"

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
