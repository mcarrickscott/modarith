
// ECDSA Implementation for curve P-256
// see https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-5.pdf
// Uses random nonce generation process as described in A.3.1
// python curve.py 64 NIST256
// This completes weierstrass.c and nist256.c for this curve. Then
// gcc -O2 nist256.c weierstrass.c hash.c -o nist256

#include <stdio.h>
#include <stdint.h>
#include "hash.h"  // Some useful hash functions

#include "curve.h"   // elliptic curve API

/*** Insert automatically generated code for P-256 prime group order group.c here ***/
/* Note that much of this code is not needed and can be deleted */

@group@

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
    if (ecnXXXisinf(P)) {
        printf("P= O\n");
    } else {
        char x[BYTES],y[BYTES];
        char buff[(2*BYTES)+1];
        ecnXXXget(P,x,y);
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
static void reduce(char *h,spint *r)
{
    int i;
    char buff[BYTES];
    gel x,y,c;
    
    mod2r(8*(BYTES-1),c); // 2^248

    for (i=0;i<BYTES-1;i++)
        buff[i]=h[i];  // little endian
    buff[BYTES-1]=0;
    reverse(buff);
    modimp(buff,y);

    for (i=0;i<9;i++)
        buff[i]=h[BYTES-1+i];
    for (i=9;i<BYTES;i++)
        buff[i]=0;
    reverse(buff);
    modimp(buff,x);

    modmul(x,c,x);  
    modadd(x,y,r); // 2^248.x + y
}

// Input private key - 32 random bytes
// Output public key - 65 bytes (0x04<x>|<y>), or 33 if compressed (0x02<x>.. or 0x03<x>)
void NIST256_KEY_PAIR(int compress,char *prv,char *pub)
{
    point P;
    ecnXXXgen(&P);

    ecnXXXmul(prv,&P); 

    if (compress) {
        pub[0]=0x02+ecnXXXget(&P,&pub[1],NULL); // 0x02 or 0x03
    } else {
        pub[0]=0x04; // no compression
        ecnXXXget(&P,&pub[1],&pub[BYTES+1]);  // get x and y
    }
}

// choice of hash functions
int NIST256_PREHASH(int sha,int mlen,char *m,char * th)
{
    int i;
    char h[64];
    if (sha==32)
    {
        hash256 sh256;
        HASH256_init(&sh256);
        for (i=0;i<mlen;i++)
            HASH256_process(&sh256,m[i]);
        HASH256_hash(&sh256,h);       
        for (i=0;i<32;i++) th[i]=h[i]; 
        return 1;
    }
    if (sha==48)
    {
        hash384 sh384;
        HASH384_init(&sh384);
        for (i=0;i<mlen;i++)
            HASH384_process(&sh384,m[i]);
        HASH384_hash(&sh384,h);       
        for (i=0;i<32;i++) th[i]=h[i];  
        return 1;
    }
    return 0;
}

// input private key, per-message random number, message to be signed. Output signature.
// ran must be Nbytes+8 in length, in this case 40 bytes
void NIST256_SIGN(char *prv,char *ran,char *thm,char *sig)
{
    char h[BYTES];
    point R;
    gel e,r,s,k;

    modimp(thm,e);

    ecnXXXgen(&R);
    modimp(prv,s);

    reduce(ran,k);
    modexp(k,h);
    ecnXXXmul(h,&R);
    modinv(k,NULL,k);

    ecnXXXget(&R,h,NULL);
    modimp(h,r);

    modmul(s,r,s);
    modadd(s,e,s);
    modmul(s,k,s);
    modzer(k);

    modexp(r,sig);
    modexp(s,&sig[BYTES]);
}

// input public key, message and signature
// NOTE signatures that are of the wrong length should be rejected prior to calling this function
int NIST256_VERIFY(char *pub,char *thm,char *sig) 
{
    point G,Q;
    int i;
    char rb[BYTES],u[BYTES],v[BYTES];
    gel e,r,s,rds;

    modimp(thm,e);

    ecnXXXgen(&G);

// import from signature
    if (!modimp(sig,r)) return 0; // if not in range
    if (!modimp(&sig[BYTES],s)) return 0;

    if (modis0(r) || modis0(s)) return 0;
    modinv(s,NULL,s);
    modmul(r,s,rds); modexp(rds,v);  // export to byte array
    modmul(s,e,s); modexp(s,u); 

    if (pub[0]==0x04) {
        ecnXXXset(0,&pub[1],&pub[BYTES+1],&Q);
    } else {
        ecnXXXset((int)pub[0]&1,&pub[1],NULL,&Q);
    }

    ecnXXXmul2(u,&G,v,&Q,&Q);
    if (ecnXXXisinf(&Q)) return 0;

    ecnXXXget(&Q,rb,NULL);

    modimp(rb,e);
    if (modcmp(r,e)) return 1;
    return 0;
}

// test for FIPS 186-5 ECDSA Signature Generation. msg is SHA256 hash

int main()
{
    const char *sk= (const char *)"519b423d715f8b581f4fa8ee59f4771a5b44c8130b4e3eacca54a56dda72b464";
    const char *ran=(const char *)"94a1bbb14b906a61a280f245f9e93c7f3b4a6247824f5d33b9670787642a68deb9670787642a68de";
    const char *msg=(const char *)"44acf6b7e36c1342c2c5897204fe09504e1e2efb1a900377dbc4e7a6a133ec56";
    char prv[BYTES],pub[2*BYTES+1];
    char buff[256],m[BYTES],k[BYTES+8],sig[2*BYTES];
    int res,compress=1;
    printf("Run test vector\n");
    printf("private key= "); puts(sk); 
    fromHex(BYTES,sk,prv);
    fromHex(BYTES+8,ran,k);
    fromHex(BYTES,msg,m);
    NIST256_KEY_PAIR(compress,prv,pub);
    if (compress)
        toHex(BYTES+1,pub,buff);
    else
        toHex(2*BYTES+1,pub,buff);

    printf("public key= "); puts(buff);
    NIST256_SIGN(prv,k,m,sig);
    toHex(2*BYTES,sig,buff);
    printf("signature=  "); puts(buff);

    res=NIST256_VERIFY(pub,m,sig);
    if (res)
        printf("Signature is valid\n");
    else
        printf("Signature is NOT valid\n");
}
