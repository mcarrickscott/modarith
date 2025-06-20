// Program to check testvectors for modular arithmetic
#include <stdio.h>
#include <stdlib.h>

/*** automatically generated code field.c for modulus inserted here ***/

@field@

/*** End of automatically generated code ***/

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

// output a modulo number in hex
static void output(spint *x) {
    char b[Nbytes+1];
    char buff[(2*Nbytes)+1];
    modexp(x,b);
    toHex(b,buff);
    puts(buff);
}

int main()
{
    char a[Nbytes],b[Nbytes],t[Nbytes];
    spint x[Nlimbs],y[Nlimbs],z[Nlimbs],w[Nlimbs];
    FILE *fp=fopen("edge.txt","r");
    size_t len;
    char line[256];
    int success=0;
    int failures=0;
    while(fgets(line, 256, fp)) {
        printf("a= %s",line);
        fromHex(line,a);
        modimp(a,x);
        if (!fgets(line, 256, fp)) break;
        printf("b= %s",line);
        fromHex(line,b);
        modimp(b,y);

        printf("Testing 1/(1/a)=a \n");
        modinv(x,NULL,z); modinv(z,NULL,z);
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("a mod p passed\n");
            success++;
        } else {
            printf("a mod p failed\n");
            failures++;
        }

        printf("Testing modadd(a,b)\n");
        modadd(x,y,z);  
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modadd(a,b) passed\n");
            success++;
        } else {
            printf("modadd(a,b) failed\n");
            failures++;
        }

        printf("Testing modsub(a,b)\n");
        modsub(x,y,z);
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modsub(a,b) passed\n");
            success++;
        } else {
            printf("modsub(a,b) failed\n");
            failures++;
        }

        printf("Testing modsub(b,a)\n");
        modsub(y,x,z); 
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modsub(b,a) passed\n");
            success++;
        } else {
            printf("modsub(b,a) failed\n");
            failures++;
        }

        printf("Testing modmul(a,b)\n");
        modmul(x,y,z);  
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modmul(a,b) passed\n");
            success++;
        } else {
            printf("modmul(a,b) failed\n");
            failures++;
        }

        printf("Testing modsqr(modsqrt(modsqr(a)))\n");
        modsqr(x,z); modsqrt(z,NULL,z); modsqr(z,z);
        if (!fgets(line, 256, fp)) break; 
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modsqr(a)   passed\n");
            success++;
        } else {
            printf("modsqr(a)   failed\n");
            failures++;
        }


        modadd(x,x,x); modadd(y,y,y); // double them and try again. Doubling will not trigger reduction.
        
        printf("Testing 1/(1/2a)=2a \n");
        modinv(x,NULL,z); modinv(z,NULL,z); 
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("2a mod p passed\n");
            success++;
        } else {
            printf("2a mod p failed\n");
            failures++;
        }

        printf("Testing modadd(2a,2b)\n");
        modadd(x,y,z); 
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modadd(2a,2b) passed\n");
            success++;
        } else {
            printf("modadd(2a,2b) failed\n");
            failures++;
        }

        printf("Testing modsub(2a,2b)\n");
        modsub(x,y,z);
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modsub(2a,2b) passed\n");
            success++;
        } else {
            printf("modsub(2a,2b) failed\n");
            failures++;
        }

        printf("Testing modsub(2b,2a)\n");
        modsub(y,x,z);
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modsub(2b,2a) passed\n");
            success++;
        } else {
            printf("modsub(2b,2a) failed\n");
            failures++;
        }

        printf("Testing modmul(2a,2b)\n");
        modmul(x,y,z); 
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modmul(2a,2b) passed\n");
            success++;
        } else {
            printf("modmul(2a,2b) failed\n");
            failures++;
        }

        printf("Testing modsqr(modsqrt(modsqr(2a)))\n");
        modsqr(x,z); modsqrt(z,NULL,z); modsqr(z,z); 
        if (!fgets(line, 256, fp)) break;
        fromHex(line,t); modimp(t,w); 
        if (modcmp(w,z)) {
            printf("modsqr(2a)   passed\n");
            success++;
        } else {
            printf("modsqr(2a)   failed\n");
            failures++;
        }
        printf("\n");
    }

    printf("Succeeded %d times\n",success);
    printf("Failed %d times\n",failures);
}
