// test and timing program for Edwards and Weierstrass curves
//
// run curve.py script to complete edwards.rs or weierstrass.rs 
// For example python curve_rust.py 64 ED25519
// make sure curve.h is in the path
//
// clang -O3 -march=native -mtune=native testcurve.c edwards.c -lcpucycles -o testcurve
// or
// clang -O3 -march=native -mtune=native testcurve.c weierstrass.c -lcpucycles -o testcurve

#include <stdio.h>
#include "curve.h"

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

#ifdef NIST256
    const char order[]="FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551"; // order
    const char r1[]=   "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B"; // random
    const char r2[]=   "E99789339379389A9F9998765C890986B39059D7021039135CE26D61EE5687D6"; // order-r1
    const char n1[]=   "C20347457078878f77b707c070707077a07707b7b07070707223252357134272"; // random
    const char n2[]=   "D35279279432f249b298a876788d86294e02842092769136c086038b1812383a"; // random
#endif
#ifdef NIST384
    const char order[]="ffffffffffffffffffffffffffffffffffffffffffffffffc7634d81f4372ddf581a0db248b0a77aecec196accc52973";
    const char r1[]=   "bd9c66b3ad3c2d6d1a3d1fa7bc8960a923b8c1e9392456de3eb13b9046685257bdd640fb06671ad11c80317fa3b1799d";
    const char r2[]=   "4263994c52c3d292e5c2e05843769f56dc473e16c6dba92188b211f1adcedb879a43ccb742498ca9d06be7eb2913afd6";
    const char n1[]=   "9a1de644815ef6d13b8faa1837f8a88b17fc695a07a0ca6e0822e8f36c031199972a846916419f828b9d2434e465e150";
    const char n2[]=   "4737819096da1dac72ff5d2a386ecbe06b65a6a48b8148f6b38a088ca65ed389b74d0fb132e706298fadc1a606cb0fb3";
#endif
#ifdef NIST521
    const char order[]="1fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffa51868783bf2f966b7fcc0148f709a5d03bb5c9b8899c47aebb6fb71e91386409";
    const char r1[]=   "d8972a846916419f828b9d2434e465e150bd9c66b3ad3c2d6d1a3d1fa7bc8960a923b8c1e9392456de3eb13b9046685257bdd640fb06671ad11c80317fa3b1799d";
    const char r2[]=   "12768d57b96e9be607d7462dbcb1b9a1eaf4263994c52c3d292e5c2e05843769f512dcdc59a860b3f8d411ac5b8b0a153787ddf88bd83352cdd9eef859eed86ea6c";
    const char n1[]=   "e5386ecbe06b65a6a48b8148f6b38a088ca65ed389b74d0fb132e706298fadc1a606cb0fb39a1de644815ef6d13b8faa1837f8a88b17fc695a07a0ca6e0822e8f3";
    const char n2[]=   "acc37459eef50bea63371ecd7b27cd813047229389571aa8766c307511b2b9437a28df6ec4ce4a2bbdc241330b01a9e71fde8a774bcf36d58b4737819096da1dac";
#endif
#ifdef NUMS256E
    const char order[]="4000000000000000000000000000000041955AA52F59439B1A47B190EEDD4AF5";
    const char r1[]=   "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
    const char r2[]=   "29978934937938999F9998765C890987383EB9CE8A51DE298370542FE0D0AD7A";
    const char n1[]=   "21347457078878f77b707c070707077a07707b7b070707072232523571342729";
    const char n2[]=   "35279279432f249b298a876788d86294e02842092769136c086038b1812383a5";
#endif
#ifdef ED25519
    const char order[]="1000000000000000000000000000000014DEF9DEA2F79CD65812631A5CF5D3ED";
    const char r1[]=   "66876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
    const char r2[]=   "9978934937938999F9998765C8909870B885907FDF03764C13B05B94EE93672";
    const char n1[]=   "20347457078878f77b707c070707077a07707b7b07070707223252357134272";
    const char n2[]=   "35279279432f249b298a876788d86294e02842092769136c086038b1812383a";
#endif
#ifdef ED448
    const char order[]="3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF7CCA23E9C44EDB49AED63690216CC2728DC58F552378C292AB5844F3";
    const char r1[]=   "1F6868465567867838578786787978654735836726356262066876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
    const char r2[]=   "209797B9AA987987C7A878798786879AB8CA7C98D9CA9D9DF997893410435C8363E873C00B5F40171816219BE8BE29E38CA165319D4BA778";
    const char n1[]=   "21347457078878f77b707c070707077a07707b7b0707070722325235713427294582948924358f9898529852985298085690868659860866";
    const char n2[]=   "35279279432f249b298a876788d86294e02842092769136c086038b1812383a5875793576874589355835357878467862937894039158588";
#endif
#ifdef NUMS256W
    const char order[]="FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE43C8275EA265C6020AB20294751A825";
    const char r1[]=   "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
    const char r2[]=   "E9978934937938999F9998765C890986DAE5E19F451EF6EE89D3C2C839450AAA";
    const char n1[]=   "120347457078878f77b707c070707077a07707b7b07070707223252357134272";
    const char n2[]=   "235279279432f249b298a876788d86294e02842092769136c086038b1812383a";    
#endif
#ifdef SECP256K1
    const char order[]="FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";
    const char r1[]=   "166876CB6C86C76660666789A376F6790956A0D6A507657196D75D610E0C9D7B";
    const char r2[]=   "E9978934937938999F9998765C890985B1583C100A413ACA28FB012BC229A3C6";
    const char n1[]=   "120347457078878f77b707c070707077a07707b7b07070707223252357134272";
    const char n2[]=   "235279279432f249b298a876788d86294e02842092769136c086038b1812383a";   
#endif

#define BYTES (sizeof(order)/2)

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
    for (i = 0; i < BYTES; i++)
    {
        unsigned char ch = src[i];
        byte2hex(&dst[i * 2],ch);
    }
    dst[2*BYTES]='\0';
}

// Convert from a hex string to byte array 
static void fromHex(const char *src, char *dst)
{
    int i,lz,len=0;
    char pad[2*BYTES];
    while (src[len]!=0) len++;
    lz=2*BYTES-len;
    if (lz<0) lz=0;
    for (i=0;i<lz;i++) pad[i]='0';  // pad with leading zeros
    for (i=lz;i<2*BYTES;i++) pad[i]=src[i-lz];

    for (i=0;i<BYTES;i++)
    {
        dst[i] = (char2int(pad[2*i]) * 16) + char2int(pad[2*i + 1]);
    }
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
        toHex(x,buff);
        printf("Px= "); puts(buff);
        toHex(y,buff);
        printf("Py= "); puts(buff);
    }
}

int main() {
    uint64_t start,fin;
    clock_t begin;
    int i,elapsed;
    char r[BYTES],a[BYTES],b[BYTES];
    point P,Q;

#ifdef NIST256
    printf("NIST256 Weierstrass curve, %d bit code\n",(int)WORDLENGTH);
#endif
#ifdef NUMS256E
    printf("NUMS256E Edwards curve, %d bit code\n",(int)WORDLENGTH);
#endif
#ifdef ED25519
    printf("ED25519 Edwards curve, %d bit code\n",(int)WORDLENGTH);
#endif
#ifdef ED448
    printf("ED448 Edwards curve, %d bit code\n",(int)WORDLENGTH);
#endif
#ifdef NUMS256W
    printf("NUMS256W Weierstrass curve, %d bit code\n",(int)WORDLENGTH);   
#endif
#ifdef SECP256K1
    printf("SECP256K1 Weierstrass curve, %d bit code\n",(int)WORDLENGTH);    
#endif

    //order=r1+r2
    printf("Generator= \n");
    ecnXXXgen(&P);
    outputxy(&P);
    
    ecnXXXcpy(&P,&Q);

    // these should be the same
    //ecnXXXadd(&Q,&P);
    //printf("P+P=\n"); outputxy(&P);
    //ecnXXXdbl(&Q);
    //printf("2P=\n"); outputxy(&Q);

    fromHex(order,r); // exponent is big endian byte array  
    ecnXXXmul(r,&P);
    if (ecnXXXisinf(&P))
        printf("MUL test passed OK\n");
    else
        printf("MUL test FAILED\n");

    fromHex(r1,a);
    fromHex(r2,b);
    ecnXXXmul2(a,&Q,b,&Q,&P);
    if (ecnXXXisinf(&P))
        printf("MUL2 test passed OK\n");
    else
        printf("MUL2 test FAILED\n");

    fromHex(n1,a); // random
    fromHex(n2,b);

    ecnXXXcpy(&Q,&P);

    printf("Timing point multiplication\n");
#ifdef COUNT_CLOCKS
#ifdef USE_RDTSC
    start=__rdtsc();
#else
    start=cpucycles();
#endif    
#endif
    begin=clock();
    for (i=0;i<10000;i++) {
        ecnXXXmul(a,&P);
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

    outputxy(&P);

    printf("Timing double point multiplication\n");
#ifdef COUNT_CLOCKS
#ifdef USE_RDTSC
    start=__rdtsc();
#else
    start=cpucycles();
#endif    
#endif
    begin=clock();
    for (i=0;i<10000;i++) {
        ecnXXXmul2(a,&P,b,&Q,&P);
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

    outputxy(&P);

    return 0;
}
