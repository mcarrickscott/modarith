# Python program to generate reasonably efficient C/C++ modular arithmetic code for any prime on a 64-bit processor
# uses unsaturated radix of 2^52
# Special version for AVX-IFMA Intel SIMD extensions
# 
# In particular this script generates code for the NIST field primes
#
# NIST256 = 2^256-2^224+2^192+2^96-1
# NIST384 = 2^384-2^128-2^96+2^32-1
#
# and the "Goldilocks" field prime
#
# X448=2^448-2^224-1
#
# requires addchain utility in the path - see https://github.com/mmcloughlin/addchain 
#
# How to use. 
# (1) First execute this program: python3 monty_ifma.py X448. Output code is written to file field.c or group.c
# (2) All constants and inputs must be converted to Montgomery nresidue form by calling nres()
# (3) All final outputs must be converted back to integer form by calling redc()
#
# The modulus can either be entered directly, or derived from the hard-wired name of an elliptic curve
# Normally this will be the field prime modulus
# However if the modulus is entered in decimal with two leading 0s, it is instead assumed to be a group 
# order (a large prime factor of the number of points on a curve)
#
# Note that even though a modulus is represented using an unsaturated base, it may still retain some shape
# For example on a 32-bit processor using a radix of 2^29 the NIST384 prime is
# {0x1FFFFFFF,0x7,0x0,0x1FFFFE00,0x1FFFEFFF,0x1FFFFFFF,0x1FFFFFFF,0x1FFFFFFF,0x1FFFFFFF,0x1FFFFFFF,0x1FFFFFFF,0x1FFFFFFF,0x1FFFFFFF,0x7F}
# which can also be represented as
# {0x1FFFFFFF,0x7,0x0,0x1FFFFE00,0x1FFFEFFF,      -0x1,       0x0,       0x0,       0x0,       0x0,       0x0,       0x0,       0x0,0x80}
#
# and now the majority of the multiplications in the Montgomery reduction are by 0, and can be eliminated.
#
# Mike Scott 22nd April 2024
# TII
#

# Some default settings

cyclesorsecs=True     # count cycles otherwise seconds
compiler="gcc" # gcc, clang or icx (inlining can sometimes cause icx to hang)
cyclescounter=True # use Bernstein's cpucycle counter, otherwise just provide timings
use_rdtsc=False # use rdtsc directly, x86 only, for better comparison with other implementations
if use_rdtsc :
    cyclescounter=False
if cyclescounter :
    use_rdtsc=False
decoration=False # decorate function names to avoid name clashes
formatted=True # pretty up the final output
inline=True # consider encouraging inlining
generic=True # set to False if algorithm is known in advance, in which case modadd and modsub can be faster - see https://eprint.iacr.org/2017/437. Set False for RFC7748 implementation.
check=False # run cppcheck on the output
scale=1 # set to 10 or 100 for faster timing loops. Default to 1
PSCR=False # Power Side Channel Resistant conditional moves and swaps

import sys
import subprocess

def ispowerof2(n) :
    if (n & (n-1) == 0) and n>0 :
        e=0
        while n>0 :
            if (n&1)!=0 :
                return e
            e=e+1
            n>>=1
    return -1

# get number of limbs
def getN(n) :
    N=int(n/base)
    if n%base!=0 :
        N=N+1
    return N

def bits(n) :
    b=0
    m=n
    while m!=0 :
        b+=1
        m>>=1
    return b

# gcd
def gcd(x, y):
    a = x
    b = y
    while b != 0:
        a, b = b, a % b
    return a

# inverse mod prime p
def inverse(a, p) :
    n = p
    x = a % n
    if x == 0:
        return x
    kn = n
    if x < 0:
        x += n
    if gcd(x, n) != 1:
        return 0
    a = 1
    la = 0
    while x > 1:
        q, r = divmod(n, x)
        t = la - a * q
        la = a
        a = t
        n = x
        x = r
    if a < 0:
        a += kn
    return a

# get Montgomery R
def getR(n) :
    m=int(n/base)
    if n%base!=0 :
        m+=1
    R=2**(m*base)
    return R

# is it a lucky trinomial?
def trinomial(p,base) :
    n=p.bit_length()
    m=2**n - p
    m = m-1
    k=20
    while k < n :
        if 2**k > m :
            return 0
        if 2**k == m: break
        k=k+1
    if k%base == 0 :
        return k//base
    return 0

# convert to radix array
def makebig(p,base,N) :
    pw=[]
    i=0
    b=2**base
    tp=p
    while i<N :
        pw.append(tp%b)
        tp=tp>>base
        i=i+1
    return pw

# process the modulus given radix and number of limbs
def process_prime(p,base,N) :
    b=2**base
    pw=makebig(p,base,N)

# Here pw is further processed to convert 2^base-1 values to -1, and carry +1 to the next digit
# generated code must be modified to process a -1 multiplier in the reduction function
# observe that for NIST primes this will increase the number of 0 entries

# ppw contains processed modulus
    ppw=[]
    for i in range(0,N) :
        ppw.append(0)

    carry=0
    i=0
    while i<N :
        if carry==1 :
            ppw[i]=pw[i]+carry
            if ppw[i] == b-1 :
                ppw[i]=-1
            else :
                carry=0
                if ppw[i]==b :
                    ppw[i]=0
                    carry=1
        else : 
            ppw[i]=pw[i]
            if PM and i==0 :
                if pw[i] == b-m :
                    ppw[i]=-m
                    carry=1
            else :
                if pw[i] == b-1 :
                    ppw[i]=-1
                    carry=1
        i=i+1
    E=False
    if carry==1 :
        E=True
        ppw.append(carry)
    return ppw,E

def intrinsics() :
    str="// set both lanes to same constant\n"
    str+="static inline spint _mm_set2_epi64(int64_t c) {\n"
    str+="\treturn _mm_set_epi64x((int64_t)c,(int64_t)c);\n"
    str+="}\n" 
    
    str+="// set each lane to a constant\n"
    str+="static inline spint _mm_setc_epi64(int64_t c0,int64_t c1) {\n"
    str+="\treturn _mm_set_epi64x(c1,c0);\n"
    str+="}\n" 

    str+="// t+=a*b\n"
    str+="static inline void accum(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\t*tl=_mm_madd52lo_epu64(*tl,a,b);\n"
    str+="\t*th=_mm_madd52hi_epu64(*th,a,b);\n"
    str+="}\n\n"

    str+="// t+=s\n"
    str+="static inline void add(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\t*tl=_mm_add_epi64(*tl,a);\n"
    str+="\t*th=_mm_add_epi64(*th,b);\n"
    str+="}\n\n"

    str+="// t-=s\n"
    str+="static inline void sub(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\t*tl=_mm_sub_epi64(*tl,a);\n"
    str+="\t*th=_mm_sub_epi64(*th,b);\n"
    str+="}\n\n"

    str+="// t+=(lo+hi)*spm. Note (lo+hi) could be > 2^52\n"
    str+="static inline void accumx(spint *tl,spint *th,spint lo,spint hi,spint spm) {\n"
    str+="\tspint s=_mm_add_epi64(lo,hi);\n"
    str+="\t*tl=_mm_madd52lo_epu64(*tl,s,spm);\n"
    str+="\t*th=_mm_madd52hi_epu64(*th,s,spm);\n"
    str+="\t*th=_mm_madd52lo_epu64(*th,_mm_srli_epi64(s,52),spm);\n"
    str+="}\n\n"

    str+="// t+= (v<<s)\n"
    str+="static inline void accumsl(spint *tl,spint *th,spint v,char s) {\n"
    str+="\tspint wl, wh;\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"
    str+="\twl=_mm_and_si128(_mm_slli_epi64(v,s),mask);\n"
    str+="\twh=_mm_srli_epi64(v,52-s);\n"
    str+="\tadd(tl,th,wl,wh);\n"
    str+="}\n\n"

    str+="// t=a*b\n"
    str+="static inline void mul(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\tspint r=_mm_setzero_si128();\n"
    str+="\t*tl=_mm_madd52lo_epu64(r,a,b);\n"
    str+="\t*th=_mm_madd52hi_epu64(r,a,b);\n"
    str+="}\n\n"

    str+="// t>>52\n"
    str+="static inline void shiftr(spint *tl,spint *th){\n"
    str+="\t*tl=_mm_srai_epi64(*tl,52);\n"  # arithmetic as could be negative
    str+="\t*tl=_mm_add_epi64(*tl,*th);\n"
    str+="\t*th=_mm_setzero_si128();\n"
    str+="}\n\n"

    str+="// r=t>>52\n"
    str+="static inline spint shiftout(spint tl,spint th,char s) {\n"
    str+="\tspint r=_mm_srai_epi64(tl,s);\n"
    str+="\tr=_mm_add_epi64(r,_mm_slli_epi64(th,52-s));\n"
    str+="\treturn r;\n"
    str+="}\n\n"

# bit of a problem here! How to efficiently multiply 64-bit value by small number

    str+="// t*=m\n"
    str+="static inline spint muls(spint tb,int64_t m) {\n"
    str+="\tspint spm=_mm_set2_epi64(m);\n"
    str+="\tspint tp=_mm_srli_epi64(tb,32);\n"
    str+="\tspint pp1=_mm_mul_epu32(tb,spm);\n"
    str+="\tspint pp2=_mm_mul_epu32(tp,spm);\n"
    str+="\treturn _mm_add_epi64(pp1,_mm_slli_epi64(pp2,32));\n"
    str+="}\n"

    str+="// return high part of a*b\n"
    str+="static inline spint mulhi(spint a,spint b) {\n"
    str+="\tspint zero=_mm_setzero_si128();\n"
    str+="\treturn _mm_madd52hi_epu64(zero,a,b);\n"
    str+="}\n"

    str+="// t*=m\n"
    str+="static inline void muli(spint *tl,spint *th,int64_t m) {\n"
    str+="\t*tl=muls(*tl,m);\n"
    str+="\t*th=muls(*th,m);\n"
    str+="}\n"

    str+="// t<<s\n"
    str+="static inline void shiftl(spint *tl,spint *th,char s) {\n"
    str+="\t*th = _mm_slli_epi64(*th,s);\n"
    str+="\t*tl = _mm_slli_epi64(*tl,s);\n"
    str+="}\n\n"

    str+="// set each lane to a constant\n"
    str+="static inline spint tospint(int c0,int c1) {\n"
    str+="\treturn _mm_set_epi64x(c1,c0);\n"
    str+="}\n" 
    return str

#conditional add of x*p
def caddp(x) :
    str=""
    for i in range(0,N) :
        if ppw[i]==0 :
            continue
        if ppw[i]==-1:
            str+="\tmpy=_mm_set2_epi64({}u);\n".format(x)
            str+="\tn[{}]=_mm_sub_epi64(n[{}],_mm_and_si128(mpy,carry));\n".format(i,i)
            #str+="\tn[{}]-=(spint){}u&carry;\n".format(i,x)
            continue
        if ppw[i]<0:
            str+="\tmpy=_mm_set2_epi64(0x{:x}u);\n".format(-ppw[i]*x)
            str+="\tn[{}]=_mm_sub_epi64(n[{}],_mm_and_si128(mpy,carry));\n".format(i,i)
            #str+="\tn[{}]-=((spint)0x{:x}u)&carry;\n".format(i,-ppw[i]*x)
            continue
        str+="\tmpy=_mm_set2_epi64(0x{:x}u);\n".format(ppw[i]*x)
        str+="\tn[{}]=_mm_add_epi64(n[{}],_mm_and_si128(mpy,carry));\n".format(i,i)
        #str+="\tn[{}]+=((spint)0x{:x}u)&carry;\n".format(i,ppw[i]*x)
    if E:
        str+="\tmpy=_mm_set2_epi64({}u*q);\n".format(x)
        str+="\tn[{}]=_mm_add_epi64(n[{}],mpy);\n".format(N-1,N-1)
        #str+="\tn[{}]+=((spint){}u*q)&carry;\n".format(N-1,x)
    return str

#add x*p
def addp(x) :
    str=""
    for i in range(0,N) :
        if ppw[i]==0 :
            continue
        if ppw[i]==-1:
            str+="\tmpy=_mm_set2_epi64({}u);\n".format(x)
            str+="\tn[{}]=_mm_sub_epi64(n[{}],mpy);\n".format(i,i)
            #str+="\tn[{}]-=(spint){};\n".format(i,x)
            continue
        if ppw[i]<0:
            str+="\tmpy=_mm_set2_epi64(0x{:x}u);\n".format(-ppw[i]*x)
            str+="\tn[{}]=_mm_sub_epi64(n[{}],mpy);\n".format(i,i)
            #str+="\tn[{}]-=((spint)0x{:x}u);\n".format(i,-ppw[i]*x)
            continue
        str+="\tmpy=_mm_set2_epi64(0x{:x}u);\n".format(ppw[i]*x)
        str+="\tn[{}]=_mm_add_epi64(n[{}],mpy);\n".format(i,i)
        #str+="\tn[{}]+=((spint)0x{:x}u);\n".format(i,ppw[i]*x)
    if E:
        str+="\tmpy=_mm_set2_epi64({}u*q);\n".format(x)
        str+="\tn[{}]=_mm_add_epi64(n[{}],mpy);\n".format(N-1,N-1)
        #str+="\tn[{}]+=((spint){}u*q);\n".format(N-1,x)
    return str

#subtract x*p
def subp(x) :
    str=""
    for i in range(0,N) :
        if ppw[i]==0 :
            continue
        if ppw[i]==-1:
            str+="\tmpy=_mm_set2_epi64({}u);\n".format(x)
            str+="\tn[{}]=_mm_add_epi64(n[{}],mpy);\n".format(i,i)
            #str+="\tn[{}]+=(spint){}u;\n".format(i,x)
            continue
        if ppw[i]<0:
            str+="\tmpy=_mm_set2_epi64(0x{:x}u);\n".format(-ppw[i]*x)
            str+="\tn[{}]=_mm_add_epi64(n[{}],mpy);\n".format(i,i)
            #str+="\tn[{}]+=(spint)0x{:x}u;\n".format(i,-ppw[i]*x)
            continue
        str+="\tmpy=_mm_set2_epi64(0x{:x}u);\n".format(ppw[i]*x)
        str+="\tn[{}]=_mm_sub_epi64(n[{}],mpy);\n".format(i,i)
        #str+="\tn[{}]-=(spint)0x{:x}u;\n".format(i,ppw[i]*x)
    if E:
        str+="\tmpy=_mm_set2_epi64({}u*q);\n".format(x)
        str+="\tn[{}]=_mm_sub_epi64(n[{}],mpy);\n".format(N-1,N-1)
        #str+="\tn[{}]-=(spint){}u*q;\n".format(N-1,x)
    return str

# Propagate carries. Modified to please MISRA requirements and clang bug
def prop(n) :
    str="//propagate carries\n"
    str+="static spint inline prop(spint *n) {\n"
    str+="\tint i;\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<{}u)-1);\n".format(base)
    #str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

    str+="\tsspint carry=(sspint)n[0];\n"
    str+="\tcarry=_mm_srai_epi64(carry,{}u);\n".format(base)
    #str+="\tcarry>>={}u;\n".format(base)

    str+="\tn[0]=_mm_and_si128(n[0], mask);\n"
    #str+="\tn[0]&=mask;\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N-1)

    str+="\t\tcarry=_mm_add_epi64(carry,n[i]);\n"
    str+="\t\tn[i]=_mm_and_si128(carry,mask);\n"
    str+="\t\tcarry=_mm_srai_epi64(carry, {}u);\n".format(base)

    #str+="\t\tcarry+=(sspint)n[i];\n"
    #str+="\t\tn[i] = (spint)carry & mask;\n"
    #str+="\t\tcarry>>={}u;\n".format(base)

    str+="\t}\n"
    str+="\tn[{}]=_mm_add_epi64(n[{}],carry);\n".format(N-1,N-1)
    #str+="\tn[{}]+=(spint)carry;\n".format(N-1)
    str+="\treturn (_mm_srai_epi64(n[{}],{}u));\n}}\n".format(N-1,WL-1)
    #str+="\treturn -((n[{}]>>1)>>{}u);\n}}\n".format(N-1,WL-2)
    return str

#propagate carries and add p if negative, propagate carries again
def flat(n) :
    str="//propagate carries and add p if negative, propagate carries again\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="spint inline flatten(spint *n) {\n"
    else :
        str+="spint flatten(spint *n) {\n"
    if E :
        str+="\tint q=(1<<{}u);\n".format(base)
    str+="\tspint one=_mm_set2_epi64(1);\n"
    str+="\tspint carry=prop(n);\n"
    str+="\tspint mpy;\n"
    str+=caddp(1)
    str+="\t(void)prop(n);\n"
    str+="\treturn _mm_and_si128(carry,one);\n"
    #str+="\treturn (carry&1);\n"
    str+="}\n"
    return str

#final subtract
def modfsb(n) : 
    str="//Montgomery final subtract\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="spint inline modfsb{}(spint *n) {{\n".format(DECOR)
    else :
        str+="spint modfsb{}(spint *n) {{\n".format(DECOR)
    if E: 
        str+="\tint q=(1<<{}u);\n".format(base)
    str+="\tspint mpy;\n" 
    str+=subp(1)
    str+="\treturn flatten(n);\n}\n"
    #str+="\treturn;\n}\n"
    return str

#modular addition
def modadd(n) :
    str="//Modular addition - reduce less than 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic :
        str+="void inline modadd{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)
    else :
        str+="void modadd{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)
    str+="\tspint mpy;\n"
    if not algorithm :
        if E: 
            str+="\tint q=(1<<{}u);\n".format(base)
        str+="\tspint carry;\n"
    for i in range(0,N) :
        str+="\tn[{}]=_mm_add_epi64(a[{}],b[{}]);\n".format(i,i,i)
        #str+="\tn[{}]=a[{}]+b[{}];\n".format(i,i,i)
    if not algorithm :
        str+=subp(2)
        str+="\tcarry=prop(n);\n"
        str+=caddp(2)
    str+="\t(void)prop(n);\n"  
    str+="}\n"    
    return str

#modular subtraction
def modsub(n) :
    str="//Modular subtraction - reduce less than 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modsub{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)
    else :
        str+="void modsub{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)

    str+="\tspint mpy;\n" 
    if not algorithm :
        if E: 
            str+="\tint q=(1<<{}u);\n".format(base)
        str+="\tspint carry;\n"
    else :
        str+="\tint q=(1<<{}u);\n".format(base)
    for i in range(0,N) :
        str+="\tn[{}]=_mm_sub_epi64(a[{}],b[{}]);\n".format(i,i,i)
        #str+="\tn[{}]=a[{}]-b[{}];\n".format(i,i,i)
    if not algorithm :
        str+="\tcarry=prop(n);\n"
        str+=caddp(2)
    else :
        str+=addp(mp)
    str+="\t(void)prop(n);\n" 
    str+="}\n"
    return str

#modular negation
def modneg(n) :
    str="//Modular negation\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modneg{}(const spint *b,spint *n) {{\n".format(DECOR)
    else :
        str+="void modneg{}(const spint *b,spint *n) {{\n".format(DECOR)
    str+="\tspint mpy;\n"
    str+="\tspint zero=_mm_setzero_si128();\n"
    if not algorithm :
        if E: 
            str+="\tint q=(1<<{}u);\n".format(base)
        str+="\tspint carry;\n"
    else :
        str+="\tint q=(1<<{}u);\n".format(base)
    for i in range(0,N) :
        str+="\tn[{}]=_mm_sub_epi64(zero,b[{}]);\n".format(i,i)
        #str+="\tn[{}]=(spint)0-b[{}];\n".format(i,i)
    if not algorithm :
        str+="\tcarry=prop(n);\n"
        str+=caddp(2)
    else :
        str+=addp(mp)
    str+="\t(void)prop(n);\n" 
    str+="}\n"
    return str

# add column of partial products from multiplication on way up
def getZMU(str,i) :
    first=True
    global maxnum

    k=0
    while (k<=i) :
        if first :
            str+="\taccum(&tl,&th,a[{}],b[{}]);".format(k,i-k)    #"\tt+=(dpint)a[{}]*b[{}];".format(k,i-k)
            first=False
        else :
            str+=" accum(&tl,&th,a[{}],b[{}]);".format(k,i-k)     #" t+=(dpint)a[{}]*b[{}];".format(k,i-k)
        k+=1
        maxnum+=maxdigit*maxdigit
    return str

# add column of partial products from multiplication on way down
def getZMD(str,i) :

    first=True
    k=i-(N-1)
    while (k<=N-1) :
        if first :
            first=False
            str+="\taccum(&tl,&th,a[{}],b[{}]);".format(k,i-k)     #"\tt+=(dpint)a[{}]*b[{}];".format(k,i-k)
        else :
            str+=" accum(&tl,&th,a[{}],b[{}]);".format(k,i-k)    #" t+=(dpint)a[{}]*b[{}];".format(k,i-k)
        k+=1
    return str

# add column of partial products from squaring on way up
def getZSU(str,i) :
    first=True
    k=0
    j=i
    hap=False
    while k<j :
        hap=True
        if first :
            str+="\tmul(&totl,&toth,a[{}],a[{}]);".format(k,i-k)    #"\ttot=(udpint)a[{}]*a[{}];".format(k,i-k)
            first=False
        else :
            str+=" accum(&totl,&toth,a[{}],a[{}]);".format(k,i-k)   #" tot+=(udpint)a[{}]*a[{}];".format(k,i-k)
        k+=1
        j-=1
    if hap:
        str+=" add(&totl,&toth,totl,toth);"        #" tot*=2;"
    if i%2==0:
        if first :
            str+="\tmul(&totl,&toth,a[{}],a[{}]);".format(int(i/2),int(i/2))      #"\ttot=(udpint)a[{}]*a[{}];".format(int(i/2),int(i/2))
        else :
            str+=" accum(&totl,&toth,a[{}],a[{}]);".format(int(i/2),int(i/2))      #" tot+=(udpint)a[{}]*a[{}];".format(int(i/2),int(i/2))
    if i==0 :
        str+="tl=totl; th=toth;"   #" t=tot;"
    else :
        str+="add(&tl,&th,totl,toth); "   #" t+=tot; "
    return str

# add column of partial products from squaring on way down
def getZSD(str,i) :
    first=True
    k=i-(N-1)
    j=N
    hap=False
    while k<i-k :
        hap=True
        if first :
            str+="\tmul(&totl,&toth,a[{}],a[{}]);".format(k,i-k)     #"\ttot=(udpint)a[{}]*a[{}];".format(k,i-k)
            first=False
        else :
            str+=" accum(&totl,&toth,a[{}],a[{}]);".format(k,i-k)    #" tot+=(udpint)a[{}]*a[{}];".format(k,i-k)
        k+=1
        j-=1
    if hap :
        str+=" add(&totl,&toth,totl,toth);"    #" tot*=2;"
    if i%2==0:
        if first :
            str+="\tmul(&totl,&toth,a[{}],a[{}]);".format(int(i/2),int(i/2))      #"\ttot=(udpint)a[{}]*a[{}];".format(int(i/2),int(i/2))
        else :
            str+=" accum(&totl,&toth,a[{}],a[{}]);".format(int(i/2),int(i/2))     #" tot+=(udpint)a[{}]*a[{}];".format(int(i/2),int(i/2))
    str+="add(&tl,&th,totl,toth); "   #" t+=tot; "
    return str

# process ppw[i] element, where i is index of the prime p, j is index of v
# do a double-shift instead of a multiply if possible
# Add in q for appearance of negative v (maybe need to add in q again for such every appearance?)
# try to accumulate contributions in single-precision if possible

def mul_process(i,j,str,gone_neg,mask_set) :
    global maxnum
    n=ppw[i]   # do nothing if n=0
    if abs(n)>1 :
        e=ispowerof2(n)
        if e > 0  :
            str+=" accumsl(&tl,&th,v{},{}u); ".format(j,e)   #" t+=(dpint)(udpint)((udpint)v{}<<{}u); ".format(j,e)
            maxnum+=maxdigit*2**e
        else :
            str+=" accum(&tl,&th,v{},p{}); ".format(j,i)  #" t+=(dpint)v{}*(dpint)p{}; ".format(j,i)
            maxnum+=maxdigit*ppw[i]
    if n == 1 :
        if mask_set :
            str+=" s=_mm_add_epi64(s,v{}); ".format(j)
            #str+=" s+=v{}; ".format(j)
        else :
            str+=" tl=_mm_add_epi64(tl,v{}); ".format(j)
            #str+=" add(&tl,&th,v{},0); ".format(j)   #" t+=(dpint)v{}; ".format(j)
            maxnum+=maxdigit
    if n == -1 :
        if mask_set :
            if not gone_neg :
                str+=" s=_mm_add_epi64(s,_mm_sub_epi64(q,v{}));".format(j)
                #str+=" s+=q-v{};".format(j)
            else :
                str+=" s=_mm_sub_epi64(s,v{}); ".format(j)
                #str+=" s-=v{}; ".format(j)
        else: 
            if not gone_neg :
                str+=" tl=_mm_add_epi64(tl,_mm_sub_epi64(q,v{}));".format(j)
                #str+=" add(&tl,&th,q-v{},0);".format(j)       #" t+=(dpint)(spint)(q-v{});".format(j)
                maxnum+=maxdigit
            else :
                str+=" tl=_mm_sub_epi64(tl,v{});".format(j)
                #str+=" sub(&tl,&th,v{},0);".format(j)       #" t-=(dpint)v{}; ".format(j)
        gone_neg=True
    return str,gone_neg

def sqr_process(i,j,str,gone_neg,mask_set) :
    global maxnum
    n=ppw[i]   # do nothing if n=0
    if abs(n)>1 :
        e=ispowerof2(n)
        if e > 0  :
            str+=" accumsl(&tl,&th,v{},{}u); ".format(j,e)  #" t+=(udpint)v{}<<{}u; ".format(j,e)
            maxnum+=maxdigit*2**e
        else :
            str+=" accum(&tl,&th,v{},p{}); ".format(j,i)    #" t+=(udpint)v{}*p{}; ".format(j,i)
            maxnum+=maxdigit*ppw[i]
    if n == 1 :
        if mask_set :
            str+=" s=_mm_add_epi64(s,v{}); ".format(j)
            #str+=" s+=v{}; ".format(j)
        else :
            str+=" tl=_mm_add_epi64(tl,v{}); ".format(j)
            #str+=" add(&tl,&th,v{},0); ".format(j)      #" t+=(udpint)v{}; ".format(j)
            maxnum+=maxdigit
    if n == -1 :
        if mask_set :
            if not gone_neg :
                str+=" s=_mm_add_epi64(s,_mm_sub_epi64(q,v{}));".format(j)
                #str+=" s+=q-v{};".format(j)
            else :
                str+=" s=_mm_sub_epi64(s,v{}); ".format(j)
                #str+=" s-=v{}; ".format(j)
        else: 
            if not gone_neg :
                str+=" tl=_mm_add_epi64(tl,_mm_sub_epi64(q,v{}));".format(j)
                #str+=" add(&tl,&th,q-v{},0);".format(j)   #" t+=(udpint)(spint)(q-v{});".format(j)
                maxnum+=maxdigit
            else :
                str+=" tl=_mm_sub_epi64(tl,v{});".format(j)
                #str+=" sub(&tl,&th,v{},0); ".format(j)    #" t-=(udpint)v{}; ".format(j)
        gone_neg=True
    return str,gone_neg

# modular multiplication, modulo p. Exploits 0 digits in p.
# Note that allowing inlining gives significant speed-up
def modmul(n) :
    global maxnum
    mask=(1<<base)-1
    s_is_declared=False
    str="// Modular multiplication, c=a*b mod 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modmul{}(const spint *a,const spint *b,spint *c) {{\n".format(DECOR)
    else :
        str+="void modmul{}(const spint *a,const spint *b,spint *c) {{\n".format(DECOR)

    str+="\tspint tl=_mm_setzero_si128(), th=_mm_setzero_si128();\n"

    for i in range(0,N) :
        if i==0 and PM :
            continue
        n=ppw[i]
        if abs(n)>1 :
            if i==0 or ispowerof2(n)<0 :
                str+="\tspint p{}=_mm_set2_epi64({}u);\n".format(i,hex(ppw[i]))
                #str+="\tspint p{}={}u;\n".format(i,hex(ppw[i]))
    str+="\tspint q=_mm_set2_epi64((int64_t)1<<{}u);\n".format(base)
    #str+="\tspint q=((spint)1<<{}u); // q is unsaturated radix \n".format(base)
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<{}u)-1);\n".format(base)
    #str+="\tspint mask=(spint)(q-(spint)1);\n"
    str+="\tspint zero=_mm_setzero_si128();\n"
    str+="\tspint one=_mm_set2_epi64(1);\n"
    if fullmonty :
        str+="\tspint ndash=_mm_set2_epi64(0x{:x}u);\n".format(ndash)
        #str+="\tspint ndash=0x{:x}u;\n".format(ndash)
    maxnum=0
    str=getZMU(str,0)
    gone_neg=False 

    if PM :
        str+="spint wl,wh;\n"
        str+="spint sm=_mm_set2_epi64({});\n".format(M)
    
    if fullmonty :
        str+=" spint v0=_mm_and_si128(_mm_madd52lo_epu64(zero,tl,ndash),mask);"
        #str+=" spint v0=((tl*ndash)&mask);"    #" spint v0=(((spint)t*ndash)&mask);"
        if PM :
            gone_neg=True
            str+=" mul(&wl,&wh,sm,_mm_sub_epi64(q,v0)); "
            str+=" add(&tl,&th,wl,wh); "
            #str+=" add(&tl,&th,(spint){}*(q-v0),0); ".format(M)    #" t+=(dpint)(spint)((spint){}*(q-v0)); ".format(M)
            maxnum+=M*maxdigit
        else :
            if ppw[0]==1 :
                str+=" tl=_mm_add_epi64(tl,v0); "
                #str+=" add(&tl,&th,v0,0); "    #" t+=(dpint)v0;"
            else :
                str+=" accum(&tl,&th,v0,p0);"   #" t+=(dpint)v0 * (dpint)p0;"
            maxnum+=ppw[0]*maxdigit
    else :
        str+="spint v0=_mm_and_si128(tl,mask);"
        #str+="spint v0=(tl & mask);"   #" spint v0=((spint)t & mask);"
    str+=" shiftr(&tl,&th);\n"  #" t>>={};\n".format(base)
    maxnum=2**(2*WL-base)
    str=getZMU(str,1)
    mask_set=False

    for i in range(1,N) :
        if gone_neg :
            if PM :
                str+=" mul(&wl,&wh,sm,mask); "
                str+=" add(&tl,&th,wl,wh);"
                #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)    #" t+=(dpint)(spint)((spint){}*mask);".format(M)
                maxnum+=M*maxdigit
            else :
                if not s_is_declared :
                    str+=" spint s=(spint)mask;"
                    s_is_declared=True
                else :
                    str+=" s=(spint)mask;"
                mask_set=True
                #str+=" t+=mask;"
        k=1
        str,gone_neg=mul_process(i,0,str,gone_neg,mask_set)
        while k<i :
            str,gone_neg=mul_process(i-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" tl=_mm_add_epi64(tl,s);"
            #str+=" add(&tl,&th,s,0);"   #" t+=(dpint)s;"
            maxnum+=maxdigit
            mask_set=False

        if fullmonty :
            str+=" spint v{}=_mm_and_si128(_mm_madd52lo_epu64(zero,tl,ndash),mask); ".format(i)
            #str+=" spint v{}=((tl*ndash) & mask); ".format(i)    #" spint v{}=(((spint)t*ndash) & mask); ".format(i)
            if PM :
                str+=" mul(&wl,&wh,sm,v{}); ".format(i)
                str+=" sub(&tl,&th,wl,wh); "
                #str+=" sub(&tl,&th,(spint){}*v{},0); ".format(M,i)    #" t-=(dpint)(spint)((spint){}*v{}); ".format(M,i)
            else :
                if ppw[0]==1 :
                    str+=" tl=_mm_add_epi64(tl,v{});".format(i)
                    #str+=" add(&tl,&th,v{},0);".format(i)  #" t+=(dpint)v{};".format(i)
                else :
                    str+=" accum(&tl,&th,v{},p0);".format(i)  #" t+=(dpint)v{} * (dpint)p0;".format(i)
                maxnum+=ppw[0]*maxdigit
        else :
            str+=" spint v{}=_mm_and_si128(tl,mask); ".format(i)
            #str+=" spint v{}=(tl & mask); ".format(i)    #" spint v{}=((spint)t & mask); ".format(i)
        if i==N-1 :
            print("// Overflow limit   =",2**(2*WL))
            print("// maximum possible =",maxnum)
            if maxnum >= 2**(2*WL) :
                print("//Warning: Overflow possibility detected - change radix ")
        str+=" shiftr(&tl,&th);\n"   #" t>>={};\n".format(base)
        maxnum=2**(2*WL-base)
        if i<N-1 :
            str= getZMU(str,i+1)
     
    str=getZMD(str,N)
    if gone_neg :
        if PM :
            str+=" mul(&wl,&wh,sm,mask); "
            str+=" add(&tl,&th,wl,wh);"
            #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)  #" t+=(dpint)(spint)((spint){}*mask);".format(M)
        else :
            if not s_is_declared :
                str+=" spint s=(spint)mask;"
                s_is_declared=True
            else :
                str+=" s=(spint)mask;"
            mask_set=True
            #str+=" t+=mask;"

    if E :  
        k=0
        while k<N :
            str,gone_neg=mul_process(N-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" tl=_mm_add_epi64(tl,s);"
            #str+=" add(&tl,&th,s,0);"  #" t+=(dpint)s;"
            mask_set=False        

        if fullmonty :
            str+=" spint v{}=_mm_and_si128(_mm_madd52lo_epu64(zero,tl,ndash),mask); ".format(N)
            #str+=" spint v{}=((tl*ndash) & mask); ".format(N)    #" spint v{}=(((spint)t*ndash) & mask); ".format(N)
            if PM :
                str+=" mul(&wl,&wh,sm,v{}); ".format(N)
                str+=" sub(&tl,&th,wl,wh); "
                #str+=" sub(&tl,&th,(spint){}*v{},0); ".format(M,N)   #" t-=(dpint)(spint)((spint){}*v{}); ".format(M,N)
            else :
                if ppw[0]==1 :
                    str+=" tl=_mm_add_epi64(tl,v{});".format(N)
                    #str+=" add(&tl,&th,v{},0);".format(N)   #" t+=(dpint)v{};".format(N)
                else :
                    str+=" accum(&tl,&th,v{},p0);".format(N)   #" t+=(dpint)v{} * (dpint)p0;".format(N)
        else :
            str+=" spint v{}=_mm_and_si128(tl,mask); ".format(N)
            #str+=" spint v{}=(tl & mask); ".format(N)    #" spint v{}=((spint)t & mask); ".format(N)
        str+=" shiftr(&tl,&th);\n"       #" t>>={};\n".format(base)

        str=getZMD(str,N+1)

        for i in range(N+1,2*N) :
            if gone_neg :
                if PM :
                    str+=" mul(&wl,&wh,sm,mask); "
                    str+=" add(&tl,&th,wl,wh);"
                    #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)  #" t+=(dpint)(spint)((spint){}*mask);".format(M)
                else :
                    if not s_is_declared :
                        str+=" spint s=(spint)mask;"
                        s_is_declared=True
                    else :
                        str+=" s=(spint)mask;"
                    mask_set=True
            k=i-N 
            while k<=N : 
                str,gone_neg=mul_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" tl=_mm_add_epi64(tl,s);"
                #str+=" add(&tl,&th,s,0);"  #" t+=(dpint)s;"
                mask_set=False        

            str+=" c[{}]=_mm_and_si128(tl,mask); ".format(i-N-1)
            #str+=" c[{}]=(tl & mask); ".format(i-N-1)    #" c[{}]=((spint)t & mask); ".format(i-N-1)
            str+=" shiftr(&tl,&th);\n"   #" t>>={};\n".format(base)

            #if i<=2*N-1 :
            if i<2*N-2 :
                str=getZMD(str,i+1)
            else :
                str+="\t"
        
        if gone_neg :
            if PM :
                str+=" tl=_mm_add_epi64(tl,_mm_sub_epi64(v{},sm));".format(N)
                #str+=" add(&tl,&th,v{}-(spint){},0);".format(N,M)     #" t+=(dpint)(spint)(v{}-(spint){});".format(N,M)
            else :
                str+=" tl=_mm_add_epi64(tl,_mm_sub_epi64(v{},one));".format(N) 
                #str+=" add(&tl,&th,v{}-(spint)1,0);".format(N)    #" t+=(dpint)(spint)(v{}-(spint)1);".format(N)
        else :
            str+="  tl=_mm_add_epi64(tl,v{});".format(N)
            #str+="  add(&tl,&th,v{},0);".format(N)  #" t+=(dpint)v{};".format(N)

        str+=" c[{}] = tl;\n".format(N-1)    #" c[{}] = (spint)t;\n".format(N-1)
    else :
        for i in range(N,2*N-1) :
            k=i-(N-1) 
            while k<=N-1 :
                str,gone_neg=mul_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" tl=_mm_add_epi64(tl,s);"
                #str+=" add(&tl,&th,s,0);"  #" t+=(dpint)s;"
                mask_set=False   

            if i==2*N-1 :
                break
            str+=" c[{}]=_mm_and_si128(tl,mask); ".format(i-N)
            #str+=" c[{}]=(tl & mask); ".format(i-N)    #" c[{}]=((spint)t & mask); ".format(i-N)
            str+=" shiftr(&tl,&th);\n" #" t>>={};\n".format(base)
            if i<=2*N-3 :
                str=getZMD(str,i+1)
                if gone_neg :
                    if PM :
                        str+=" mul(&wl,&wh,sm,mask); "
                        str+=" add(&tl,&th,wl,wh);"
                        #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)   #" t+=(dpint)(spint)((spint){}*mask);".format(M)
                    else :
                        if not s_is_declared :
                            str+=" spint s=(spint)mask;"
                            s_is_declared=True
                        else :
                            str+=" s=(spint)mask;"
                        mask_set=True
        if gone_neg :
            if PM :
                str+="\ttl=_mm_sub_epi64(tl,sm);"
                #str+="\tsub(&tl,&th,{},0);".format(M)   #"\tt-=(dpint){};".format(M)
            else :
                str+="\ttl=_mm_sub_epi64(tl,one);"
                #str+="\tsub(&tl,&th,1u,0);"  #"\tt-=(dpint)1u;"
        str+="\tc[{}] = tl;\n".format(N-1)   #"\tc[{}] = (spint)t;\n".format(N-1)
    str+="}\n"
    return str


# multiply by an integer
def modmli(n) :
    N=getN(n)
    mask=(1<<base)-1
    xcess=N*base-n
    str="// Modular multiplication by an integer, c=a*b mod 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modmli{}(const spint *a,int b,spint *c) {{\n".format(DECOR)
    else :
        str+="void modmli{}(const spint *a,int b,spint *c) {{\n".format(DECOR)
    str+="\tspint bw=_mm_set2_epi64(b);\n"
    if xcess>0 :
        smask=(1<<(base-xcess))-1
        str+="\tspint smask=_mm_set2_epi64({});\n".format(smask)
        str+="\tspint zero=_mm_setzero_si128();\n"
    if trin>0 :
        str+="\tspint tl=_mm_set2_epi32(0),th=_mm_set2_epi32(0);\n"

        str+="\tspint s;\n"
        str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<{}u)-1);\n".format(base)
        #str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

        for i in range(0,N) :
            str+="\taccum(&tl,&th,a[{}],bw); ".format(i)     #"\tt+=(udpint)a[{}]*(udpint)b; ".format(i)
            str+="c[{}]=_mm_and_si128(tl,mask); shiftr(&tl,&th);\n".format(i)
            #str+="c[{}]=tl & mask; shiftr(&tl,&th,{}u);\n".format(i,base)    #"c[{}]=(spint)t & mask; t=t>>{}u;\n".format(i,base)

        str+="// reduction pass\n\n"  
        str+="\ts=tl;\n"  

        if xcess>0 :
            smask=(1<<(base-xcess))-1
            str+= "\ts=_mm_add_epi64(_mm_slli_epi64(s,{}),_mm_srli_epi64(c[{}],{}u)); c[{}]=_mm_and_si128(c[{}],smask);\n".format(xcess,N-1,base-xcess,N-1,N-1)
            #str+= "\ts=(s<<{})+(c[{}]>>{}u); c[{}]&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)
        str+="\tc[0]=_mm_add_epi64(c[0],s);\n"
        #str+="\tc[0]+=s;\n"
        str+="\tc[{}]=_mm_add_epi64(c[{}],s);\n".format(trin)
        #str+="\tc[{}]+=s;\n".format(trin)
    else :        
# Barrett/Dhem method
#        str+="/*\n"
        for i in range(0,N) :
            d=ppw[i]
            if abs(d)>1 :
                if d>1 :
                    if i==0 or ispowerof2(d)<0 :
                        str+="\tspint p{}=_mm_set2_epi64({}u);\n".format(i,hex(d))
                if d<1 :
                    if i==0 :
                        str+="\tspint p{}=_mm_set2_epi64({}u);\n".format(i,hex(-d))
        str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<{}u)-1);\n".format(base)
        #str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)
        str+="\tspint tl=_mm_setzero_si128(),th=_mm_setzero_si128();\n"
        str+="\tspint q,h;\n"
        str+="\tspint r=_mm_set2_epi64(0x{:x});\n".format((2**(n+base-2))//p)
        for i in range(0,N-1) :
            str+="\taccum(&tl,&th,a[{}],bw); ".format(i) #"\tt+=(udpint)a[{}]*(udpint)b; ".format(i)
            str+="c[{}]=_mm_and_si128(tl,mask); shiftr(&tl,&th);\n".format(i)
            #str+="c[{}]=tl & mask; shiftr(&tl,&th,{}u);\n".format(i,base) #"c[{}]=(spint)t & mask; t=t>>{}u;\n".format(i,base)
        str+="\taccum(&tl,&th,a[{}],bw); ".format(N-1)   #"\tt+=(udpint)a[{}]*(udpint)b; ".format(N-1)
        str+="c[{}]=tl;\n".format(N-1)  #"c[{}]=(spint)t & mask;\n".format(N-1)
        
        str+="\t\n//Barrett-Dhem reduction\n"
        str+="\th=shiftout(tl,th,{}u);\n".format((n-2)%base)
        #str+="\tshiftr(&tl,&th,{}u); h=tl;\n".format((n-WL)%base)   #"\th = (spint)(t>>{}u);\n".format((n-WL)%base)
        str+="\tq=mulhi(h,r);\n"     #"\tq=(spint)(((udpint)h*(udpint)r)>>{}u);\n".format(WL)  

# required to propagate carries?

        for i in range(0,N) :
            d=ppw[i]
            if d==0 :
                continue
            if d==-1 :
                str+="\tc[{}]=_mm_add_epi64(c[{}],q);\n".format(i,i)
                #str+="\tc[{}]+=q;\n".format(i)
                continue
            if d==1 :
                str+="\tc[{}]=_mm_sub_epi64(c[{},q);\n".format(i,i)
                #str+="\tc[{}]-=q;\n".format(i)
            e=ispowerof2(d)
            if e>0 :
                if i<N-1 :
                    str+="\ttl=q; th=zero; shiftl(&tl,&th,{}u); c[{}]=_mm_sub_epi64(c[{}],_mm_and_si128(tl,mask)); shiftr(&tl,&th); c[{}]=_mm_sub_epi64(c[{}],tl);\n".format(e,i,i,i+1,i+1) 
                    #str+="\ttl=q; shiftl(&tl,&th,{}u); c[{}]-=tl&mask; shiftr(&tl,&th,{}u); c[{}]-=tl;\n".format(e,i,base,i+1)   #"\tt=(udpint)q<<{}u; c[{}]-=(spint)t&mask; c[{}]-=(spint)(t>>{}u);\n".format(e,i,i+1,base)
                else :
                    str+="\tc[{}]=_mm_sub_epi64(c[{}],_mm_slli_epi64(q,{}u));\n".format(i,i,e)
                    #str+="\tc[{}]-=q<<{}u;\n".format(i,e)
            else :
                if d<0 :
                    str+="\tmul(&tl,&th,q,p{}); c[{}]=_mm_add_epi64(c[{}],_mm_and_si128(tl,mask)); shiftr(&tl,&th); c[{}]=_mm_add_epi64(c[{}],tl);\n".format(i,i,i,i+1,i+1)
                    #str+="\tmul(&tl,&th,q,p{}); c[{}]+=tl&mask; shiftr(&tl,&th,{}u); c[{}]+=tl;\n".format(i,i,base,i+1)  #"\tt=(udpint)q*p{}; c[{}]+=(spint)t&mask; c[{}]+=(spint)(t>>{}u);\n".format(i,i,i+1,base)
                else :
                    if i<N-1 :
                        str+="\tmul(&tl,&th,q,p{}); c[{}]=_mm_sub_epi64(c[{}],_mm_and_si128(tl,mask)); shiftr(&tl,&th); c[{}]=_mm_sub_epi64(c[{}],tl);\n".format(i,i,i,i+1,i+1)
                        #str+="\tmul(&tl,&th,q,p{}); c[{}]-=tl&mask; shiftr(&tl,&th,{}u); c[{}]-=tl;\n".format(i,i,base,i+1)   #"\tt=(udpint)q*p{}; c[{}]-=(spint)t&mask; c[{}]-=(spint)(t>>{}u);\n".format(i,i,i+1,base)
                    else :
                        str+="\ttl=_mm_madd52lo_epu64(zero,q,p{});\n".format(i)
                        str+="\tc[{}]=_mm_sub_epi64(c[{}],tl);\n".format(i,i)
                        #str+="\tc[{}]-=q*p{};\n".format(i,i)
                        #str+="\tc[{}]=(c[{}]-(q*p{}))&mask;\n".format(i,i,i)
#        str+="*/\n"
#        str+="\tspint t[{}];\n".format(N)
#        str+="\tmodint{}(b,t);\n".format(DECOR)
#        str+="\tmodmul{}(a,t,c);\n".format(DECOR)
    str+="\t(void)prop(c);\n"   # fully propagate carries
    str+="}\n"
    return str

# modular squaring
# Note that allowing inlining gives significant speed-up
def modsqr(n) :
    mask=(1<<base)-1
    s_is_declared=False
    str="// Modular squaring, c=a*a  mod 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modsqr{}(const spint *a,spint *c) {{\n".format(DECOR)
    else :
        str+="void modsqr{}(const spint *a,spint *c) {{\n".format(DECOR)

    str+="\tspint totl,toth;\n"   #"\tudpint tot;\n"
    str+="\tspint tl=_mm_setzero_si128(),th=_mm_setzero_si128();\n"   #"\tudpint t=0;\n"
    for i in range(0,N) :
        if i==0 and PM :
            continue
        n=ppw[i]
        if abs(n)>1 :
            if i==0 or ispowerof2(n)<0 :
                str+="\tspint p{}=_mm_set2_epi64({}u);\n".format(i,hex(ppw[i]))
    str+="\tspint q=_mm_set2_epi64((int64_t)1<<{}u);\n".format(base)
    #str+="\tspint q=((spint)1<<{}u); // q is unsaturated radix \n".format(base)
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<{}u)-1);\n".format(base)
    #str+="\tspint mask=(spint)(q-(spint)1);\n"
    str+="\tspint zero=_mm_setzero_si128();\n"
    str+="\tspint one=_mm_set2_epi64(1);\n"
    if fullmonty :
        str+="\tspint ndash=_mm_set2_epi64(0x{:x}u);\n".format(ndash)
        #str+="\tspint ndash=0x{:x}u;\n".format(ndash)
    str=getZSU(str,0)
    gone_neg=False

    if PM :
        str+="spint wl,wh;\n"
        str+="spint sm=_mm_set2_epi64({});\n".format(M)
    
    if fullmonty :
        str+=" spint v0=_mm_and_si128(_mm_madd52lo_epu64(zero,tl,ndash),mask);"
        #str+=" spint v0=((tl*ndash)& mask);"    #" spint v0=(((spint)t*ndash)& mask);"
        if PM :
            gone_neg=True
            str+=" mul(&wl,&wh,sm,_mm_sub_epi64(q,v0)); "
            str+=" add(&tl,&th,wl,wh); "
            #str+=" add(&tl,&th,(spint){}*(q-v0),0); ".format(M)  #" t+=(udpint)(spint)((spint){}*(q-v0)); ".format(M)
        else :
            if ppw[0]==1 :
                str+=" tl=_mm_add_epi64(tl,v0);"
                #str+=" add(&tl,&th,v0,0);"  #" t+=(udpint)v0;"
            else :
                str+=" accum(&tl,&th,v0,p0);"  #" t+=(udpint)v0 * p0;"

    else :
        str+=" spint v0=_mm_and_si128(tl,mask);"
        #str+=" spint v0=(tl & mask);"   #" spint v0=((spint)t & mask);"
    str+=" shiftr(&tl,&th);\n"  #" t>>={};\n".format(base)
    
    str=getZSU(str,1)
    mask_set=False

    for i in range(1,N) :
        if gone_neg :
            if PM :
                str+=" mul(&wl,&wh,sm,mask); "
                str+=" add(&tl,&th,wl,wh);"
                #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)   #" t+=(udpint)(spint)((spint){}*mask);".format(M)
            else :
                if not s_is_declared :
                    str+=" spint s=(spint)mask;"
                    s_is_declared=True
                else :
                    str+=" s=(spint)mask;"
                mask_set=True

        k=1
        str,gone_neg=sqr_process(i,0,str,gone_neg,mask_set)
        while k<i :
            str,gone_neg=sqr_process(i-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" tl=_mm_add_epi64(tl,s);"
            #str+=" add(&tl,&th,s,0);"   #" t+=(udpint)s;"
            mask_set=False   

        if fullmonty :
            str+=" spint v{}=_mm_and_si128(_mm_madd52lo_epu64(zero,tl,ndash),mask); ".format(i)
            #str+=" spint v{}=((tl*ndash) & mask); ".format(i)    #" spint v{}=(((spint)t*ndash) & mask); ".format(i)
            if PM :
                str+=" mul(&wl,&wh,sm,v{}); ".format(i)
                str+=" sub(&tl,&th,wl,wh); "
                #str+=" sub(&tl,&th,(spint){}*v{},0); ".format(M,i)  #" t-=(udpint)(spint)((spint){}*v{}); ".format(M,i)
            else :
                if ppw[0]==1 :
                    str+=" tl=_mm_add_epi64(tl,v{});".format(i)
                    #str+=" add(&tl,&th,v{},0);".format(i)   #" t+=(udpint)v{};".format(i)
                else :
                    str+=" accum(&tl,&th,v{},p0);".format(i)  #" t+=(udpint)v{} * p0;".format(i)
        else :
            str+=" spint v{}=_mm_and_si128(tl,mask);".format(i)
            #str+=" spint v{}=(tl & mask);".format(i)  #" spint v{}=((spint)t & mask);".format(i)
        str+=" shiftr(&tl,&th);\n"    #" t>>={};\n".format(base)
        if i<N-1 :
            str= getZSU(str,i+1)
     
    str=getZSD(str,N)
    if gone_neg :
        if PM :
            str+=" mul(&wl,&wh,sm,mask); "
            str+=" add(&tl,&th,wl,wh);"
            #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)  #" t+=(udpint)(spint)((spint){}*mask);".format(M)
        else :
            if not s_is_declared :
                str+=" spint s=(spint)mask;"
                s_is_declared=True
            else :
                str+=" s=(spint)mask;"
            mask_set=True

    if E :  
        k=0
        while k<N :
            str,gone_neg=sqr_process(N-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" tl=_mm_add_epi64(tl,s);"
            #str+=" add(&tl,&th,s,0);"  #" t+=(udpint)s;"
            mask_set=False   
        if fullmonty :
            str+=" spint v{}=_mm_and_si128(_mm_madd52lo_epu64(zero,tl,ndash),mask); ".format(N)
            #str+=" spint v{}=((tl*ndash) & mask); ".format(N)   #" spint v{}=(((spint)t*ndash) & mask); ".format(N)
            if PM :
                str+=" mul(&wl,&wh,sm,v{}); ".format(N)
                str+=" sub(&tl,&th,wl,wh); "
                #str+=" sub(&tl,&th,(spint){}*v{},0); ".format(M,N)   #" t-=(udpint)(spint)((spint){}*v{}); ".format(M,N)
            else :
                if ppw[0]==1 :
                    str+=" tl=_mm_add_epi64(tl,v{});".format(N)
                    #str+=" add(&tl,&th,v{},0);".format(N)  #" t+=(udpint)v{};".format(N)
                else :
                    str+=" accum(&tl,&th,v{},p0);".format(N)  #" t+=(udpint)v{} * p0;".format(N)
        else :
            str+=" spint v{}=_mm_and_si128(tl,mask); ".format(N)
            #str+=" spint v{}=(tl & mask); ".format(N)   #" spint v{}=((spint)t & mask); ".format(N)
        str+=" shiftr(&tl,&th);\n"      #" t>>={};\n".format(base)

        str=getZSD(str,N+1)
        for i in range(N+1,2*N) :
            if gone_neg :
                if PM :
                    str+=" mul(&wl,&wh,sm,mask); "
                    str+=" add(&tl,&th,wl,wh);"
                    #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)  #" t+=(udpint)(spint)((spint){}*mask);".format(M)
                else :
                    if not s_is_declared :
                        str+=" spint s=(spint)mask;"
                        s_is_declared=True
                    else :
                        str+=" s=(spint)mask;"
                    mask_set=True
                    #str+=" t+=mask;"
            k=i-N 
            while k<=N :
                str,gone_neg=sqr_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" tl=_mm_add_epi64(tl,s);"
                #str+=" add(&tl,&th,s,0);"  #" t+=(udpint)s;"
                mask_set=False  
            str+=" c[{}]=_mm_and_si128(tl,mask); ".format(i-N-1)
            #str+=" c[{}]=(tl & mask); ".format(i-N-1)    #" c[{}]=((spint)t & mask); ".format(i-N-1)
            str+=" shiftr(&tl,&th);\n"   #" t>>={};\n".format(base)

            #if i<=2*N-1 :
            if i<2*N-2 :
                str=getZSD(str,i+1)
            else :
                str+="\t"

        if gone_neg :
            if PM :
                str+=" tl=_mm_add_epi64(tl,_mm_sub_epi64(v{},sm));".format(N)
                #str+=" add(&tl,&th,v{}-(spint){},0);".format(N,M)  #" t+=(udpint)(spint)(v{}-(spint){});".format(N,M)
            else :
                str+=" tl=_mm_add_epi64(tl,_mm_sub_epi64(v{},one));".format(N)
                #str+=" add(&tl,&th,v{}-(spint)1,0);".format(N)     #" t+=(udpint)(spint)(v{}-(spint)1);".format(N)
        else :
            str+=" tl=_mm_add_epi64(tl,v{});".format(N)
            #str+=" add(&tl,&th,v{},0);".format(N)  #" t+=(udpint)v{};".format(N)
        str+=" c[{}] = tl;\n".format(N-1)    #" c[{}] = (spint)t;\n".format(N-1)
    else :
        for i in range(N,2*N-1) :
            k=i-(N-1) 
            while k<=N-1 :
                str,gone_neg=sqr_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" tl=_mm_add_epi64(tl,s);"
                #str+=" add(&tl,&th,s,0);"   #" t+=(udpint)s;"
                mask_set=False   
            if i==2*N-1 :
                break
            str+=" c[{}]=_mm_and_si128(tl,mask); ".format(i-N)  
            #str+=" c[{}]=(tl & mask); ".format(i-N)   #" c[{}]=((spint)t & mask); ".format(i-N)
            str+=" shiftr(&tl,&th);\n"  #" t>>={};\n".format(base)
            if i<=2*N-3 :
                str=getZSD(str,i+1)
                if gone_neg :
                    if PM :
                        str+=" mul(&wl,&wh,sm,mask); "
                        str+=" add(&tl,&th,wl,wh);"
                        #str+=" add(&tl,&th,(spint){}*mask,0);".format(M)  #" t+=(udpint)(spint)((spint){}*mask);".format(M)
                    else :
                        if not s_is_declared :
                            str+=" spint s=(spint)mask;"
                            s_is_declared=True
                        else :
                            str+=" s=(spint)mask;"
                        mask_set=True
        if gone_neg :
            if PM :
                str+="\ttl=_mm_sub_epi64(tl,sm);"
                #str+="\tsub(&tl,&th,{},0);".format(M)   #"\tt-=(udpint){};".format(M)
            else :
                str+="\ttl=_mm_sub_epi64(tl,one);"
                #str+="\tsub(&tl,&th,1u,0);"   #"\tt-=1u;"
        str+="\tc[{}] = tl;\n".format(N-1)    #"\tc[{}] = (spint)t;\n".format(N-1)
    str+="}\n"
    return str

def modcpy() :
    str="//copy\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modcpy{}(const spint *a,spint *c) {{\n".format(DECOR)
    else :
        str+="void modcpy{}(const spint *a,spint *c) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\tc[i]=a[i];\n"
    str+="\t}\n"
    str+="}\n"
    return str

def modnsqr() :
    str="//square n times\n"
    if makestatic :
        str+="static "
    str+="void modnsqr{}(spint *a,int n) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<n;i++) {\n"
    str+="\t\tmodsqr{}(a,a);\n".format(DECOR)
    str+="\t}\n"
    str+="}\n"
    return str

# uses https://github.com/mmcloughlin/addchain to create addition chain
def modpro() :
    cline="addchain search {} > inv.acc".format(PE)
    subprocess.call(cline, shell=True)
    subprocess.call("addchain gen inv.acc > ac.txt", shell=True)
    subprocess.call("rm inv.acc",shell=True)

    f=open('ac.txt')
    lines=f.readlines()
    info=lines[0].split()
    ntmps=len(info)-1
    str="//Calculate progenitor\n"
    if makestatic :
        str+="static "
    str+="void modpro{}(const spint *w,spint *z) {{\n".format(DECOR)
    str+="\tspint x[{}];\n".format(N)
    for i in range(0,ntmps) :
        str+="\tspint t{}[{}];\n".format(i,N)  
    str+="\tmodcpy{}(w,x);\n".format(DECOR)
    for i in range(1,len(lines)) :
        info=lines[i].split()
        if info[0]=="double" :
            str+="\tmodsqr{}({},{});\n".format(DECOR,info[2],info[1])
        if info[0]=="add" :
            str+="\tmodmul{}({},{},{});\n".format(DECOR,info[2],info[3],info[1])
        if info[0]=="shift" :
            if info[1]!=info[2] :
                str+="\tmodcpy{}({},{});\n".format(DECOR,info[2],info[1])
            str+="\tmodnsqr{}({},{});\n".format(DECOR,info[1],int(info[3]))
            #str+="\tfor (i=0;i<{};i++) {{\n".format(int(info[3]))
            #str+="\t\tmodsqr{}({},{});\n".format(DECOR,info[1],info[1])
            #str+="\t}\n"
    str+="}\n"
    f.close()
    subprocess.call("rm ac.txt",shell=True)    
    return str

def modinv() :
    str="//calculate inverse, provide progenitor h if available\n"
    if makestatic :
        str+="static "
    str+="void modinv{}(const spint *x,const spint *h,spint *z) {{\n".format(DECOR)
    if PM1D2>1 :
        str+="\tint i;\n"
    str+="\tspint s[{}];\n".format(N)
    str+="\tspint t[{}];\n".format(N)
    str+="\tif (h==NULL) {\n"
    str+="\t\tmodpro{}(x,t);\n".format(DECOR)
    str+="\t} else {\n"
    str+="\t\tmodcpy{}(h,t);\n".format(DECOR)
    str+="\t}\n"
    str+="\tmodcpy{}(x,s);\n".format(DECOR)
    if PM1D2>1 :
        str+="\tfor (i=0;i<({}-1);i++) {{\n".format(PM1D2)
        str+="\t\tmodsqr{}(s,s);\n".format(DECOR)
        str+="\t\tmodmul{}(s,x,s);\n".format(DECOR)
        str+="\t}\n"
    str+="\tmodnsqr{}(t,{});\n".format(DECOR,PM1D2+1)
    #str+="\tfor (i=0;i<={};i++) {{\n".format(PM1D2)
    #str+="\t\tmodsqr{}(t,t);\n".format(DECOR)
    #str+="\t}\n"
    str+="\tmodmul{}(s,t,z);\n".format(DECOR)
    str+="}\n"
    return str

def modqr() :
    str="//Test for quadratic residue \n"
    if makestatic :
        str+="static "
    str+="spint modqr{}(const spint *h,const spint *x) {{\n".format(DECOR)
    str+="\tspint r[{}];\n".format(N)
    str+="\tif (h==NULL) {\n"
    str+="\t\tmodpro{}(x,r);\n".format(DECOR)
    str+="\t\tmodsqr{}(r,r);\n".format(DECOR)
    str+="\t} else {\n"
    str+="\t\tmodsqr{}(h,r);\n".format(DECOR)
    str+="\t}\n"
    str+="\tmodmul{}(r,x,r);\n".format(DECOR)
    if PM1D2>1 :
        str+="\tmodnsqr{}(r,{});\n".format(DECOR,PM1D2-1)
    str+="\treturn modis1{}(r) | modis0{}(x);\n}}\n".format(DECOR,DECOR)
    return str

#modular square root
def modsqrt() :
    str="//Modular square root, provide progenitor h if available, NULL if not\n"
    if makestatic :
        str+="static "
    str+="void modsqrt{}(const spint *x,const spint *h,spint *r) {{\n".format(DECOR)
    if PM1D2>1 :
        str+="\tint k;\n"
        str+="\tspint t[{}];\n".format(N)
        str+="\tspint b[{}];\n".format(N)
        str+="\tspint v[{}];\n".format(N)
        str+="\tspint z[{}]={{".format(N)
        for i in range(0,N-1) :
            str+=hex(ROI[i])
            str+="u,"
        str+=hex(ROI[N-1])
        str+="u};\n"
    str+="\tspint s[{}];\n".format(N)
    str+="\tspint y[{}];\n".format(N)
    str+="\tif (h==NULL) {\n"
    str+="\t\tmodpro{}(x,y);\n".format(DECOR)
    str+="\t} else {\n"
    str+="\t\tmodcpy{}(h,y);\n".format(DECOR)
    str+="\t}\n" 
    str+="\tmodmul{}(y,x,s);\n".format(DECOR)
    if PM1D2>1 :
        str+="\tmodmul{}(s,y,t);\n".format(DECOR)
        str+="\tnres{}(z,z);\n".format(DECOR)
        str+="\tfor (k={};k>1;k--) {{\n".format(PM1D2)
        str+="\t\tmodcpy{}(t,b);\n".format(DECOR)
        str+="\t\tmodnsqr{}(b,k-2);\n".format(DECOR)

        str+="\t\tspint one=_mm_set2_epi64(1);\n";
        str+="\t\tspint d=_mm_sub_epi64(one,modis1{}(b));\n".format(DECOR)

        #str+="\t\tint d=1-modis1{}(b);\n".format(DECOR)
        str+="\t\tmodmul{}(s,z,v);\n".format(DECOR)
        str+="\t\tmodcmv{}(d,v,s);\n".format(DECOR)
        str+="\t\tmodsqr{}(z,z);\n".format(DECOR)
        str+="\t\tmodmul{}(t,z,v);\n".format(DECOR)
        str+="\t\tmodcmv{}(d,v,t);\n".format(DECOR)
        str+="\t}\n" 
    str+="\tmodcpy{}(s,r);\n".format(DECOR)
    str+="}\n"
    return str

# test for unity
def modis1(n) :
    str="//is unity?\n"
    if makestatic :
        str+="static "
    str+="spint modis1{}(const spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tspint c0;\n"
    str+="\tspint one=_mm_set2_epi64(1);\n"
    str+="\tspint d=_mm_setzero_si128();\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\td=_mm_or_si128(d,c[i]);\n\t}\n"
    #str+="\tc0=(spint)c[0];\n"
    str+="\treturn _mm_and_si128(_mm_and_si128(one,_mm_srli_epi64(_mm_sub_epi64(d,one),{}u)),_mm_srli_epi64(_mm_sub_epi64(_mm_xor_si128(c[0],one),one),{}u));\n}}\n".format(base,base)
    #str+="\treturn ((spint)1 & ((d-(spint)1)>>{}u) & (((c0^(spint)1)-(spint)1)>>{}u));\n}}\n".format(base,base)
    return str

#test for zero
def modis0(n) :
    str="//is zero?\n"
    if makestatic :
        str+="static "
    str+="spint modis0{}(const spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tspint one=_mm_set2_epi64(1);\n"
    str+="\tspint d=_mm_setzero_si128();\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\td=_mm_or_si128(d,c[i]);\n\t}\n" 
    str+="\treturn _mm_and_si128(_mm_srli_epi64(_mm_sub_epi64(d,one),{}u),one);\n}}\n".format(base)
    return str

#set to zero
def modzer() :
    str="//set to zero\n"
    if makestatic :
        str+="static "
    str+="void modzer{}(spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=_mm_setzero_si128();\n"
    str+="\t}\n"
    str+="}\n"
    return str

def modone() :
    str="//set to one\n"
    if makestatic :
        str+="static "
    str+="void modone{}(spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\ta[0]=_mm_set2_epi64(1);\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=_mm_setzero_si128();\n"
    str+="\t}\n"
    str+="\tnres{}(a,a);\n".format(DECOR)
    str+="}\n"
    return str

def modint() :
    str="//set to integer\n"
    if makestatic :
        str+="static "
    str+="void modint{}(int x,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\ta[0]=_mm_set2_epi64(x);\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=_mm_setzero_si128();\n"
    str+="\t}\n"
    str+="\tnres{}(a,a);\n".format(DECOR)
    str+="}\n"
    return str

# convert to Montgomery nresidue form
def nres(n) :
    str="//Convert m to n-residue form, n=nres(m) \n"
    if makestatic :
        str+="static "
    str+="void nres{}(const spint *m,spint *n) {{\n".format(DECOR)
    str+="\tspint c[{}];\n".format(N)
    for i in range(0,N-1) :
        str+="\t\tc[{}]=_mm_set2_epi64({});\n".format(i,hex(cw[i]))
   
    str+="\tmodmul{}(m,c,n);\n".format(DECOR)
    str+="}\n"
    return str 

#convert back from Montgomery to integer form
def redc(n) :
    str="//Convert n back to normal form, m=redc(n) \n"
    if makestatic :
        str+="static "
    str+="void redc{}(const spint *n,spint *m) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tc[0]=_mm_set2_epi64(1);\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\tc[i]=_mm_setzero_si128();\n"
    str+="\t}\n"
    str+="\tmodmul{}(n,c,m);\n".format(DECOR)
    str+="\t(void)modfsb{}(m);\n".format(DECOR)
    str+="}\n"
    return str 

#Reduce double length mod p
#def modred(n) :
#    str="//reduce double length input to n-residue \n"
#    if makestatic :
#        str+="static "
#    str+="void modred{}(const spint *n,spint *b) {{\n".format(DECOR)
#    str+="\tint i;\n"
#    if E :
#        str+="\tconst spint h[{}]={{".format(N)
#        for i in range(0,N-1) :
#            str+=hex(twopn[i])
#            str+="u,"
#        str+=hex(twopn[N-1])
#        str+="u};\n"
#    str+="\tspint t[{}];\n".format(N)
#    str+="\tfor (i=0;i<{};i++) {{\n".format (N)
#    str+="\t\tb[i]=n[i];\n"
#    str+="\t\tt[i]=n[i+{}];\n".format(N)
#    str+="\t}\n"
#    str+="\tnres{}(t,t);\n".format(DECOR)
#    if E:
#        str+="\tmodmul{}(t,h,t);\n".format(DECOR)
#        str+="\tnres{}(b,b);\n".format(DECOR)
#        str+="\tmodadd{}(b,t,b);\n".format(DECOR)
#    else :
#        str+="\tmodadd{}(b,t,b);\n".format(DECOR)
#        str+="\tnres{}(b,b);\n".format(DECOR)
#    str+="}\n"
#    return str 

#conditional swap -  see Santos and Scott eprint 2025
def modcsw() :
    str="//conditional swap g and f if d=1\n"
    str+="//strongly recommend inlining be disabled using compiler specific syntax\n"
    if makestatic :
        str+="static "
    str+="void modcsw{}(spint b,volatile spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    if PSCR :
        str+="\tspint c0,c1,s,t,w,v,aux;\n"
        str+="\tstatic uint64_t R0=0,R1=0;\n"
        str+="\tspint zero=_mm_setzero_si128();\n"
        str+="\tspint one=_mm_set2_epi64(1);\n"
        str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"
        str+="\tR0+=0x3cc3c33c5aa5a55au;\n"
        str+="\tR1+=0x7447e88e1ee1e11eu;\n" 

        str+="\tw=_mm_setc_epi64(R0,R1);\n"
        str+="\tc0=_mm_andnot_si128(b,_mm_add_epi64(w,one));\n"
        str+="\tc1=_mm_add_epi64(b,w);\n" 
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\ts=g[i]; t=f[i];\n"
        str+="\t\tv=_mm_madd52lo_epu64(_mm_madd52lo_epu64(zero,w,t),w,s);\n"
        #str+="\t\tv=w*(t+s);\n"
        str+="\t\tf[i]=aux=_mm_madd52lo_epu64(_mm_madd52lo_epu64(zero,c0,t),c1,s);\n"
        #str+="\t\tf[i] = aux = c0*t+c1*s;\n"
        str+="\t\tf[i]=_mm_and_si128(_mm_sub_epi64(aux,v),mask);\n"
        #str+="\t\tf[i] = aux - v;\n"
        str+="\t\tg[i]=aux=_mm_madd52lo_epu64(_mm_madd52lo_epu64(zero,c0,s),c1,t);\n"
        #str+="\t\tg[i] = aux = c0*s+c1*t;\n"
        str+="\t\tg[i]=_mm_and_si128(_mm_sub_epi64(aux,v),mask);\n"
        #str+="\t\tg[i] = aux - v;\n\t}\n"
        str+="\t}\n"
        str+="}\n"
    else :
        str+="\tspint zero=_mm_setzero_si128();\n"
        str+="\tspint delta,mask=_mm_sub_epi64(zero,b);\n"
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\tdelta=_mm_and_si128(_mm_xor_si128(g[i],f[i]),mask);\n"
        str+="\t\tg[i]=_mm_xor_si128(g[i],mask);\n"
        str+="\t\tf[i]=_mm_xor_si128(f[i],mask);\n\t}\n"
        str+="}\n"
    return str

#conditional move
def modcmv() :
    str="//conditional move g to f if d=1\n"
    str+="//strongly recommend inlining be disabled using compiler specific syntax\n"
    if makestatic :
        str+="static "
    str+="void modcmv{}(spint b,const spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    if PSCR :
        str+="\tspint c0,c1,s,t,v,w,aux;\n"
        str+="\tstatic uint64_t R0=0,R1=0;\n"
        str+="\tspint zero=_mm_setzero_si128();\n"
        str+="\tspint one=_mm_set2_epi64(1);\n"
        str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"
        str+="\tR0+=0x3cc3c33c5aa5a55au;\n"
        str+="\tR1+=0x7447e88e1ee1e11eu;\n" 

        str+="\tw=_mm_setc_epi64(R0,R1);\n"
        str+="\tc0=_mm_andnot_si128(b,_mm_add_epi64(w,one));\n"
        str+="\tc1=_mm_add_epi64(b,w);\n" 
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\ts=g[i]; t=f[i];\n"
        str+="\t\tv=_mm_madd52lo_epu64(_mm_madd52lo_epu64(zero,w,t),w,s);\n"

        str+="\t\tf[i]=aux=_mm_madd52lo_epu64(_mm_madd52lo_epu64(zero,c0,t),c1,s);\n"
        str+="\t\tf[i]=_mm_and_si128(_mm_sub_epi64(aux,v),mask);\n"
        str+="\t}\n"
        str+="}\n"
    else :
        str+="\tspint zero=_mm_setzero_si128();\n"
        str+="\tspint delta,mask=_mm_sub_epi64(zero,b);\n"
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\tdelta=_mm_and_si128(_mm_xor_si128(g[i],f[i]),mask);\n"
        str+="\t\tf[i]=_mm_xor_si128(f[i],mask);\n\t}\n"
        str+="}\n"
    return str

#shift left
def modshl(n) :
    N=getN(n)
    str="//shift left by less than a word\n"
    if makestatic :
        str+="static "
    str+="void modshl{}(unsigned int n,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"
    str+="\ta[{}]=_mm_or_si128(_mm_slli_epi64(a[{}],n),_mm_srli_epi64(a[{}],{}u-n));\n".format(N-1,N-1,N-2,base)
    #str+="\ta[{}]=((a[{}]<<n)) | (a[{}]>>({}u-n));\n".format(N-1,N-1,N-2,base)
    str+="\tfor (i={};i>0;i--) {{\n".format(N-2)
    str+="\t\ta[i]=_mm_or_si128(_mm_and_si128(_mm_slli_epi64(a[i],n),mask),_mm_srli_epi64(a[i-1],{}u-n));\n\t}}\n".format(base)
    #str+="\t\ta[i]=((a[i]<<n)&(spint)0x{:x}) | (a[i-1]>>({}u-n));\n\t}}\n".format(mask,base)
    str+="\ta[0]=_mm_and_si128(_mm_slli_epi64(a[0],n),mask);\n"
    #str+="\ta[0]=(a[0]<<n)&(spint)0x{:x};\n".format(mask)
    str+="}\n"
    return str 

#shift right
def modshr(n) :
    N=getN(n)
    str="//shift right by less than a word. Return shifted out part\n"
    if makestatic :
        str+="static "
    str+="spint modshr{}(unsigned int n,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"
    str+="\tspint mskn=_mm_set2_epi64((1<<n)-1);\n"    
    str+="\tspint r=_mm_and_si128(a[0],mskn);\n"

    #str+="\tspint r=a[0]&(((spint)1<<n)-(spint)1);\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N-1)
    str+="\t\ta[i]=_mm_or_si128(_mm_srli_epi64(a[i],n),_mm_and_si128(_mm_slli_epi64(a[i+1],{}u-n),mask));\n\t}}\n".format(base)
    #str+="\t\ta[i]=(a[i]>>n) | ((a[i+1]<<({}u-n))&(spint)0x{:x});\n\t}}\n".format(base,mask)
    str+="\ta[{}]=_mm_srli_epi64(a[{}],n);\n".format(N-1,N-1)
    #str+="\ta[{}]=a[{}]>>n;\n".format(N-1,N-1)
    str+="\treturn r;\n}\n"
    return str


def mod2r() :
    str="//set a= 2^r\n"
    if makestatic :
        str+="static "
    str+="void mod2r{}(unsigned int r,spint *a) {{\n".format(DECOR)
    str+="\tunsigned int n=r/{}u;\n".format(base)
    str+="\tunsigned int m=r%{}u;\n".format(base)
    str+="\tspint one=_mm_set2_epi64(1);\n"
    str+="\tmodzer{}(a);\n".format(DECOR)
    str+="\tif (r>={}*8) return;\n".format(Nbytes)
    str+="\ta[n]=one; a[n]=_mm_slli_epi64(a[n],m);\n"
    str+="nres{}(a,a);\n}}\n".format(DECOR)
    return str

#export to byte array
def modexp() :
    str="//export to byte array\n"
    if makestatic :
        str+="static "
    str+="void modexp{}(const spint *a,char *b,char *e) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i={};i>=0;i--) {{\n".format(Nbytes-1)
    str+="\t\tb[i]=_mm_extract_epi16(c[0],0)&0xff;\n"
    str+="\t\tif (e!=NULL) e[i]=_mm_extract_epi16(c[0],4)&0xff;\n"
    #str+="\t\telse e[i]=0;\n"
    #str+="\t\tb[i]=c[0]&(spint)0xff;\n"
    str+="\t\t(void)modshr{}(8,c);\n\t}}\n".format(DECOR)
    str+="}\n"
    return str 

#import from byte array
def modimp() :
    str="//import from byte array\n"
    str+="//returns 1 if in range, else 0\n"
    if makestatic :
        str+="static "
    str+="spint modimp{}(const char *b, const char *e, spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint res;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=_mm_setzero_si128();\n\t}\n"
    #str+="\t\ta[i]=0;\n\t}\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(Nbytes)
    str+="\t\tmodshl{}(8,a);\n".format(DECOR)
    str+="\t\tif (e!=NULL) a[0]=_mm_add_epi64(a[0],_mm_setc_epi64((unsigned char)b[i],(unsigned char)e[i]));\n"
    str+="\t\telse a[0]=_mm_add_epi64(a[0],_mm_setc_epi64((unsigned char)b[i],0));\n"
    str+="\t}\n"
    str+="\tres=modfsb{}(a);\n".format(DECOR)
    str+="\tnres{}(a,a);\n".format(DECOR)
    str+="\treturn res;\n"
    str+="}\n"
    return str 

#get sign (parity of value)
def modsign() :
    str="//determine sign\n"
    if makestatic :
        str+="static "
    str+="spint modsign{}(const spint *a) {{\n".format(DECOR)
    str+="\tspint c[{}];\n".format(N)
    str+="\tspint one=_mm_set2_epi64(1);\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\treturn _mm_and_si128(c[0],one);\n"
    str+="}\n"
    return str

#compare for equality
def modcmp() :
    str="//return true if equal\n"
    if makestatic :
        str+="static "
    str+="spint modcmp{}(const spint *a,const spint *b) {{\n".format(DECOR)
    str+="\tspint c[{}],d[{}];\n".format(N,N)
    str+="\tint i;\n"
    str+="\tspint one=_mm_set2_epi64(1);\n"
    str+="\tspint eq=one;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tredc{}(b,d);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\teq=_mm_and_si128(eq, _mm_and_si128(_mm_srli_epi64(_mm_sub_epi64(_mm_xor_si128(c[i],d[i]),one),{}),one));\n\t}}\n".format(base)

    str+="\treturn eq;\n"
    str+="}\n"
    return str

# for timings
def time_modmul(n,ra,rb) :
    N=getN(n)
    mask=(1<<base)-1
    bit=n%base
    smask=(1<<bit)-1
    rap=makebig(ra,base,N)
    rbp=makebig(rb,base,N)
    str="void time_modmul{}() {{\n".format(DECOR)
    str+="\tspint x[{}],y[{}],z[{}];\n".format(N,N,N)
    str+="\tint i,j;\n"

    str+="\tuint64_t start,finish;\n"
    str+="\tclock_t begin;\n"
    str+="\tint elapsed;\n"

    str+="\t"
    for i in range(0,N) :
        str+="x[{}]=_mm_setc_epi64({},{}); ".format(i,hex(rap[i]),hex(rbp[i]))
        #str+="x[{}]={}; ".format(i,hex(rap[i]))
    str+="\n\t"
    for i in range(0,N) :
        str+="y[{}]=_mm_setc_epi64({},{}); ".format(i,hex(rbp[i]),hex(rap[i]))
        #str+="y[{}]={}; ".format(i,hex(rbp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    str+="\tnres{}(y,y);\n".format(DECOR)

    if cyclescounter :
        str+="\tstart=cpucycles();\n"
    if use_rdtsc :
        str+="\tstart=__rdtsc();\n"
    str+="\tbegin=clock();\n"
    str+="\tfor (i=0;i<{};i++)\n".format(100000//scale)
    str+="\t\tfor (j=0;j<200;j++) {\n"
    str+="\t\t\tmodmul{}(x,y,z);\n".format(DECOR)
    str+="\t\t\tmodmul{}(z,x,y);\n".format(DECOR)
    str+="\t\t\tmodmul{}(y,z,x);\n".format(DECOR)
    str+="\t\t\tmodmul{}(x,y,z);\n".format(DECOR)
    str+="\t\t\tmodmul{}(z,x,y);\n".format(DECOR)
    str+="\t\t}\n"
    if cyclescounter :
        str+="\tfinish=cpucycles();\n"
    if use_rdtsc :
        str+="\tfinish=__rdtsc();\n"
    str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10*scale)
    str+="\tredc{}(z,z);\n".format(DECOR)
 
    if cyclescounter or use_rdtsc :
        str+='\tprintf("modmul check 0x%06x Clock cycles= %d Nanosecs= %d\\n",_mm_extract_epi16(z[0],0)&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
    else :
        str+='\tprintf("modmul check 0x%06x Nanosecs= %d\\n",_mm_extract_epi16(z[0],0)&0xFFFFFF,elapsed);\n'
    str+="}\n"
    return str

def time_modsqr(n,r) :
    N=getN(n)
    mask=(1<<base)-1
    bit=n%base
    smask=(1<<bit)-1
    rp=makebig(r,base,N)
    str="void time_modsqr{}() {{\n".format(DECOR)
    str+="\tspint x[{}],z[{}];\n".format(N,N)
    str+="\tint i,j;\n"

    str+="\tuint64_t start,finish;\n"
    str+="\tclock_t begin;\n"
    str+="\tint elapsed;\n"

    str+="\t"
    for i in range(0,N) :
        str+="x[{}]=_mm_setc_epi64({},{}); ".format(i,hex(rp[i]),hex(rp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    if cyclescounter :
        str+="\tstart=cpucycles();\n"
    if use_rdtsc :
        str+="\tstart=__rdtsc();\n" 
    str+="\tbegin=clock();\n"
    str+="\tfor (i=0;i<{};i++)\n".format(100000//scale)
    str+="\t\tfor (j=0;j<500;j++) {\n"
    str+="\t\t\tmodsqr{}(x,z);\n".format(DECOR)
    str+="\t\t\tmodsqr{}(z,x);\n".format(DECOR)
    str+="\t\t}\n"
    if cyclescounter :
        str+="\tfinish=cpucycles();\n"
    if use_rdtsc :
        str+="\tfinish=__rdtsc();\n"
    str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10*scale)
    str+="\tredc{}(z,z);\n".format(DECOR)
 
    if cyclescounter or use_rdtsc :
        str+='\tprintf("modsqr check 0x%06x Clock cycles= %d Nanosecs= %d\\n",_mm_extract_epi16(z[0],0)&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
    else :
        str+='\tprintf("modsqr check 0x%06x Nanosecs= %d\\n",_mm_extract_epi16(z[0],0)&0xFFFFFF,elapsed);\n'
    str+="}\n"
    return str

def time_modinv(n,r) :
    N=getN(n)
    mask=(1<<base)-1
    bit=n%base
    smask=(1<<bit)-1
    rp=makebig(r,base,N)
    str="void time_modinv{}() {{\n".format(DECOR)
    str+="\tspint x[{}],z[{}];\n".format(N,N)
    str+="\tint i,j;\n"
 
    str+="\tuint64_t start,finish;\n"
    str+="\tclock_t begin;\n"
    str+="\tint elapsed;\n"

    str+="\t"
    for i in range(0,N) :
        str+="x[{}]=_mm_setc_epi64({},{}); ".format(i,hex(rp[i]),hex(rp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    if cyclescounter :
        str+="\tstart=cpucycles();\n"
    if use_rdtsc :
        str+="\tstart=__rdtsc();\n"
    str+="\tbegin=clock();\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(50000//scale)
    str+="\t\t\tmodinv{}(x,NULL,z);\n".format(DECOR)
    str+="\t\t\tmodinv{}(z,NULL,x);\n".format(DECOR)
    str+="\t}\n"
    if cyclescounter :
        str+="\tfinish=cpucycles();\n"
    if use_rdtsc :
        str+="\tfinish=__rdtsc();\n"
    str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10000*scale)
    str+="\tredc{}(z,z);\n".format(DECOR)
 
    if cyclescounter or use_rdtsc:
        str+='\tprintf("modinv check 0x%06x Clock cycles= %d Nanosecs= %d\\n",_mm_extract_epi16(z[0],0)&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000//scale)
    else :
        str+='\tprintf("modinv check 0x%06x Microsecs= %d\\n",_mm_extract_epi16(z[0],0)&0xFFFFFF,elapsed);\n'
    str+="}\n"
    return str

def header() :
    print("\n//Automatically generated modular arithmetic C code")
    print("//Command line : python {} {}".format(sys.argv[0], sys.argv[1]))
    print("//Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)\n")
    print("#include <stdio.h>")
    print("#include <stdint.h>\n")
    print("#include <emmintrin.h>")
    print("#include <smmintrin.h>\n")
    print("#include <immintrin.h>\n")
    print("#define sspint __m128i")
    print("#define spint __m128i")

    print("#define Wordlength{} {}".format(DECOR,WL))
    print("#define Nlimbs{} {}".format(DECOR,N))
    print("#define Radix{} {}".format(DECOR,base))
    print("#define Nbits{} {}".format(DECOR,n))
    print("#define Nbytes{} {}\n".format(DECOR,Nbytes))
    print("#define MONTGOMERY")
    if prime[0].isalpha() :
        print("#define",prime.upper(),"\n")
    if trin>0 :
        print("#define MULBYINT")


def functions() :
    print(intrinsics())
    print(prop(n))
    print(flat(n))
    print(modfsb(n))
    print(modadd(n))
    print(modsub(n))
    print(modneg(n))
    print(modmul(n))
    print(modsqr(n))
    print(modcpy())
    print(modnsqr())
    print(modpro())
    print(modinv())
    print(nres(n))
    print(redc(n))
    #print(modred(n))
    print(modis1(n))
    print(modis0(n))
    print(modzer())
    print(modone())
    print(modint())
    print(modmli(n))
    print(modqr())
    print(modcmv())
    print(modcsw())
    print(modsqrt())
    print(modshl(n))
    print(modshr(n))
    print(mod2r())
    print(modexp())
    print(modimp())
    print(modsign())
    print(modcmp())

def main() :
    str="int main() {\n"
    str+='\tprintf("C code - timing some functions - please wait\\n");\n'
    str+="\ttime_modmul();\n"
    str+="\ttime_modsqr();\n"
    str+="\ttime_modinv();\n"
    str+="\treturn 0;\n"
    str+="}\n"
    return str

if len(sys.argv)!=2 :
    print("Syntax error")
    print("Valid syntax - python monty_ifma.py <prime> OR <prime name>")
    print("For example - python monty_ifma.py NIST256")
    print("For example - python monty_ifma.py 0x01fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffa51868783bf2f966b7fcc0148f709a5d03bb5c9b8899c47aebb6fb71e91386409")
    exit(2)

WL=64

prime=sys.argv[1]
p=0
base=52

algorithm=False
mp=0

# Note that number base MUST be less than wordlength WL (unsaturated).
# The radix for an n-bit prime is 2^base where base is usually in the range 51-60 for a 64-bit processor, 26-30 for 
# 32-bit and 12-13 for 16-bit.
# The radix must be low enough so that the sum of L 2^(2*base) double length integers does not encroach on the sign bit, 
# that is bits(L)+2*base<2*WL, where L=n/base (rounded up) is the number of limbs in the number representation. But 
# at the same time radix should be high enough to minimize L. An underused top digit allows for lazy reduction.

# More named primes can be added here

### Start of user editable area


if prime=="NIST256" :
    p=2**256-2**224+2**192+2**96-1

if prime=="NIST384" :
    p=2**384-2**128-2**96+2**32-1

if prime=="X25519" :
    p=2**255-19

if prime=="ED448" :
    p=2**448-2**224-1
    if not generic :
        algorithm=True
        mp=4            # Assuming Edwards curve - see https://eprint.iacr.org/2017/437

if prime=="X448" :
    p=2**448-2**224-1
    if not generic :
        algorithm=True  # if algorithm is known, fix multiple of prime (for modular subtractions) as described in https://eprint.iacr.org/2017/437
        mp=2            # Make sure there is sufficient excess - otherwise change default base. Here assuming Montgomery ladder algorithm RFC7748. 
                        # Now no reduction required after modular additions/subtractions.

if prime=="NIST521" :
    p=2**521-1

if prime=="GM270" :
    p=2**270-2**162-1

if prime=="GM240" :
    p=2**240-2**183-1

if prime=="PM383" :
    p=2**383-187

if prime=="GM360" :
    p=2**360-2**171-1

if prime=="GM480" :
    p=2**480-2**240-1

if prime=="GM378" :
    p=2**378-2**324-1

if prime=="GM384" :
    p=2**384-2**186-1

if prime=="GM512" :
    p=2**512-2**127-1

if prime=="C41417" :
    p=2**414-17

if prime=="PM512" :
    p=2**512-569

if prime=="PM266" :
    p=2**266-3

if prime=="PM336" :
    p=2**336-3

if prime=="NIST224" :
    p=2**224-2**96+1

if prime=="SECP256K1" :
    p=2**256-2**32-977

if prime=="TWEEDLE" :
    p=0x40000000000000000000000000000000038aa127696286c9842cafd400000001

if prime=="SIDH434" :
    p=2**216*3**137-1

if prime=="SIDH503" :
    p=2**250*3**159-1
 
if prime=="SIDH610" :
    p=2**305*3**192-1

if prime=="SIDH751" :
    p=2**372*3**239-1

if prime=="MFP4" :
    p=3*67*(2**246)-1

if prime=="MFP7" :
    p=2**145*(3**9)*(59**3)*(311**3)*(317**3)*(503**3)-1

if prime=="MFP1973" :
    p=0x34e29e286b95d98c33a6a86587407437252c9e49355147ffffffffffffffffff

if prime=="SQISIGN_1" :
    p=5*2**248-1

if prime=="SQISIGN_2" :
    p=65*2**376-1

if prime=="SQISIGN_3" :
    p=27*2**500-1

if prime=="CSIDH512" :
    p=5326738796327623094747867617954605554069371494832722337612446642054009560026576537626892113026381253624626941643949444792662881241621373288942880288065659

### End of user editable area

field=True
noname=False
if p==0 :
    if not prime[0].isdigit() :
        print("This named prime not supported")
        exit(2)
    if prime[0]=='0' and prime[1]=='0' :
        field=False
        p=int(prime)
    else :
        p=eval(prime)
    noname=True    # unnamed prime

fnamec='field.c'
fnameo='field.o'
if not field :
    fnamec='group.c'
    fnameo='group.o'

n=p.bit_length() 
if n<120 or pow(3,p-1,p)!=1 :
#if pow(3,p-1,p)!=1 :
    print("Not a sensible modulus, too small or not a prime")
    exit(2)

PM=False


trin=trinomial(p,base)
if trin>0 :
    print("Prime is a lucky trinomial")


b=2**base
N=getN(n)  # Number of digits in modulus (to given radix)
xcess=N*base-n  # excess bits in representation - should be 2 or greater, or 0 in which case extra "virtual" limb will be added

if N>9 :
    inline=False

# Is it an exploitable Pseudo-Mersenne?
m=2**n-p
if m>1 and bits(m)+base<WL :    # product of m and 2**base must fit in a word for this optimization to work
    PM=True
    print("Exploitable Pseudo Mersenne detected")

# find biggest power of 2 which divides p-1
p1=p-1
PM1D2=0
while (p1%2)==0 :
    PM1D2+=1
    p1>>=1
e=(1<<PM1D2)
PE=(p-1-e)//(2*e)  # exponent for use in inversion, QR check, and for square roots

# get number of bytes for export
Nbytes=n//8
if (n%8)!=0 :
    Nbytes+=1

#get a non-trivial root of unity
qnr=0
roi=0
if PM1D2==1: 
    roi=p-1
if PM1D2==2:
    roi=pow(2,(p-1)//4,p)
if PM1D2>2 : 
    qnr=2
    while pow(qnr,(p-1)//2,p)==1 :
        qnr+=1
    roi=pow(qnr,(p-1)//e,p)

# convert to radix representation
ROI=makebig(roi,base,N)

mod8=p%8
print("Prime is of length",n,"bits and =",mod8,"mod 8. Chosen radix is",base,"bits, using",N,"limbs with excess of",xcess,"bits")
print("Using standard Comba for modmul")

# process prime, check if "virtual" extra limb required
ppw,E=process_prime(p,base,N)

# get M for pseudo-Mersenne 2^n-M
M=0
if PM :
    M=-ppw[0]

minus_ones=0
print("Displaying prime limbs (we like 0, -1, 1 and powers of 2)")
for i in range(0,len(ppw)) :
    if i>0 and ppw[i]==-1 :
        minus_ones+=1
    print(hex(ppw[i]))      

if minus_ones>1 :
    print("Sorry - too many -1s for this script to handle")
    exit(1)

# get Montgomery modulus R
if E :
    print("Extra virtual limb added")
    R=getR(n+base) # add extra "virtual" limb
else :
    R=getR(n)

twopn=makebig(((1<<n)*R)%p,base,N)  # nres(2^n)

#print("twppn= ",twopn)

# should only be an issue for bad user choice of radix
if xcess<2 and (not E) :
    print("Error - Excess is only one bit, consider change of radix")
    exit(1)

# calculate ndash constant
b=2**base
ndash=inverse((R-p)%b,b)

# check if Montgomery friendly
fullmonty=True
if ndash==1 :
    fullmonty=False

# represent Montgomery constant to given radix
c=(R*R)%p              # Montgomery constant
cw=makebig(c,base,N)

#pm1=makebig(p-1,base,N)

maxdigit=2**base-1
maxnum=0

# generate fused WL-bit modular multiplication and squaring code for shaped prime
# where the unsaturated radix is 2^base and base<WL

from contextlib import redirect_stdout

# Note that the accumulated partial products must not exceed the double precision limit. 
# this limit can be 2**(2*WL)
# use unsigned integer types to store limbs

makestatic=False
DECOR=""
modulus=p

import random

makestatic=True
random.seed(42)
ra=random.randint(0,modulus-1)
rb=random.randint(0,modulus-1)
rs=random.randint(0,modulus-1)
ri=random.randint(0,modulus-1)

subprocess.call("rm time.c", shell=True)

with open('time.c', 'w') as f:
    with redirect_stdout(f):
        header()
        if cyclescounter :
            print("#include <cpucycles.h>\n")
        print("#include <time.h>\n")
        if use_rdtsc :
            print("#include <x86intrin.h>\n")
        print(intrinsics())
        print(prop(n))
        print(flat(n))
        print(modfsb(n))
        print(modmul(n))
        print(nres(n))
        print(redc(n))
        print(modint())
        print(modmli(n))
        print(modsqr(n))
        print(modcpy())
        print(modnsqr())
        print(modpro())
        print(modinv())
        print(time_modmul(n,ra,rb))
        print(time_modsqr(n,rs))
        print(time_modinv(n,ri))
        print(main())
f.close()

if cyclescounter :
    subprocess.call(compiler + " -march=native -mtune=native -O3 time.c -lcpucycles -o time", shell=True)
else :
    subprocess.call(compiler + " -march=native -mtune=native -O3 time.c -o time", shell=True)
print("For timings run ./time")

# to determine code size
makestatic=False
with open(fnamec, 'w') as f:
    with redirect_stdout(f):
        header()
        functions()
f.close()

if decoration :
    if PM :
        if noname :
            DECOR="_"+str(n)+str(m)+"_ct"
        else :
            DECOR="_"+prime+"_ct"
    else :
        if noname :
            print("Modulus must have a name - unable to make one for you")
            exit(1)
        DECOR="_"+prime+"_ct"

makestatic=True
with open(fnamec, 'w') as f:
    with redirect_stdout(f):
        header()
        functions()

f.close()

if formatted :
    subprocess.call("clang-format -i "+fnamec, shell=True)  # tidy up the format

if check: 
    subprocess.call("cppcheck --enable=all --addon=misc  --suppress=unusedFunction --suppress=missingIncludeSystem --suppress=checkersReport "+fnamec, shell=True)  # tidy up the format


if field :
    print("field code is in field.c")
else :
    print("group code is in group.c")

sys.exit(base)
