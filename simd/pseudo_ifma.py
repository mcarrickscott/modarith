# Python program to generate reasonably efficient C/C++ modular arithmetic code for pseudo-mersenne primes on a 64-bit processor
# Modulus should be a pseudo-mersenne of the form 2^n-m
# uses unsaturated radix of 2^52
# Special version for AVX-IFMA Intel SIMD extensions
#
# In particular this script generates code for primes like
#
# X25519 = 2^255-19
# PM266 = 2^266-3
# NIST521 = 2^521-1
#
# requires addchain utility in the path - see https://github.com/mmcloughlin/addchain 
#
# Execute this program as: python3 pseudoms64.py X25519
# Production code is output to file field.c
# 
# Mike Scott 19th February 2025
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

import sys
import subprocess

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

def intrinsics() :

    str="// set both lanes to same constant\n"
    str+="static inline spint _mm_set2_epi64(int64_t c) {\n"
    str+="\treturn _mm_set_epi64x((int64_t)c,(int64_t)c);\n"
    str+="}\n" 
    
    str+="// set each lane to a constant\n"
    str+="static inline spint _mm_setc_epi64(int c0,int c1) {\n"
    str+="\treturn _mm_set_epi64x(c1,c0);\n"
    str+="}\n" 

    str+="// t+=a*b\n"
    str+="static inline void accum(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\t*tl=_mm_madd52lo_epu64(*tl,a,b);\n"
    str+="\t*th=_mm_madd52hi_epu64(*th,a,b);\n"
    str+="}\n\n"

    str+="// t+=(lo+hi)*spm. Note (lo+hi) could be > 2^52\n"
    str+="static inline void accumx(spint *tl,spint *th,spint lo,spint hi,spint spm) {\n"
    str+="\tspint s=_mm_add_epi64(lo,hi);\n"
    str+="\t*tl=_mm_madd52lo_epu64(*tl,s,spm);\n"
    str+="\t*th=_mm_madd52hi_epu64(*th,s,spm);\n"
    str+="\t*th=_mm_madd52lo_epu64(*th,_mm_srli_epi64(s,52),spm);\n"
    str+="}\n\n"

    str+="// t=a*b\n"
    str+="static inline void mul(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\tspint r=_mm_set2_epi64(0);\n"
    str+="\t*tl=_mm_madd52lo_epu64(r,a,b);\n"
    str+="\t*th=_mm_madd52hi_epu64(r,a,b);\n"
    str+="}\n\n"

    str+="// t>>52\n"
    str+="static inline void shiftr(spint *tl,spint *th){\n"
    str+="\t*tl=_mm_srai_epi64(*tl,52);\n"  # arithmetic as could be negative
    str+="\t*tl=_mm_add_epi64(*tl,*th);\n"
    str+="\t*th=_mm_set2_epi64(0);\n"
    str+="}\n\n"

    str+="// r=t>>52\n"
    str+="static inline spint shiftout(spint tl,spint th) {\n"
    str+="\tspint r=_mm_srai_epi64(tl,52);\n"
    str+="\tr=_mm_add_epi64(r,th);\n"
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

    str+="// t*=m\n"
    str+="static inline void muli(spint *tl,spint *th,int64_t m) {\n"
    str+="\t*tl=muls(*tl,m);\n"
    str+="\t*th=muls(*th,m);\n"
    str+="}\n"

    str+="// t<<s\n"
    str+="static inline void shiftl(spint *tl,spint *th,char s){\n"
    str+="\t*th = _mm_slli_epi64(*th,s);\n"
    str+="\t*tl = _mm_slli_epi64(*tl,s);\n"
    str+="}\n\n"

    str+="// set each lane to a constant\n"
    str+="static inline __m128i tospint(int c0,int c1) {\n"
    str+="\treturn _mm_set_epi64x(c1,c0);\n"
    str+="}\n" 


    return str

#conditional add of p
def caddp() :
    str=""
    str+="\tbot=_mm_and_si128(bot,carry);\n"
    str+="\ttop=_mm_and_si128(top,carry);\n"
    str+="\tn[0]=_mm_sub_epi64(n[0],bot);\n"
    str+="\tn[{}]=_mm_add_epi64(n[{}],top);\n".format(N-1,N-1)
    return str

#add p
def addp() :
    str=""
    str+="\tn[0]=_mm_sub_epi64(n[0],bot);\n"
    str+="\tn[{}]=_mm_add_epi64(n[{}],top);\n".format(N-1,N-1) 
    return str

#subtract p
def subp() :
    str=""
    str+="\tn[0]=_mm_add_epi64(n[0],bot);\n"
    str+="\tn[{}]=_mm_sub_epi64(n[{}],top);\n".format(N-1,N-1)  
    return str

# Propagate carries. 
def prop(n) :
    str="//propagate carries\n"
    str+="static spint inline prop(spint *n) {\n"
    str+="\tint i;\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"

    #str+="\tsspint carry=(sspint)n[0]>>{}u;\n".format(base)
    str+="\tsspint carry=(sspint)n[0];\n"
    str+="\tcarry=_mm_srai_epi64(carry,{}u);\n".format(base)

    str+="\tn[0]=_mm_and_si128(n[0],mask);\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N-1)

    str+="\t\tcarry=_mm_add_epi64(carry,n[i]);\n"
    str+="\t\tn[i]=_mm_and_si128(carry,mask);\n"
    str+="\t\tcarry=_mm_srai_epi64(carry,{}u);\n".format(base)
    str+="\t}\n"
    str+="\tn[{}]=_mm_add_epi64(n[{}],carry);\n".format(N-1,N-1)
    str+="\treturn (_mm_srai_epi64(n[{}],{}u));\n}}\n".format(N-1,WL-1)
    return str


#propagate carries and add p if negative, propagate carries again
def flat(n) :
    mask=(1<<base)-1
    str="//propagate carries and add p if negative, propagate carries again\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="spint inline flatten(spint *n) {\n"
    else :
        str+="spint flatten(spint *n) {\n"

    str+="\tspint bot=_mm_set2_epi64({}u);\n".format(m) 
    str+="\tspint top=_mm_set2_epi64(0x{:x}u);\n".format(TW)  

    str+="\tspint carry=prop(n);\n"
    str+=caddp()
    str+="\t(void)prop(n);\n"
    str+="\tspint mask=_mm_set2_epi64(1);\n"
    str+="\treturn _mm_and_si128(carry,mask);\n"
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
    str+="\tspint bot=_mm_set2_epi64({}u);\n".format(m) 
    str+="\tspint top=_mm_set2_epi64(0x{:x}u);\n".format(TW) 
    str+=subp()
    str+="\treturn flatten(n);\n"
    str+="}\n"
    return str

#modular addition
def modadd(n,m) :
    str="//Modular addition - reduce less than 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic :
        str+="void inline modadd{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)
    else :
        str+="void modadd{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)
    if not algorithm :
        str+="\tspint carry;\n"
    str+="\tspint bot=_mm_set2_epi64({}u);\n".format(m*2) 
    str+="\tspint top=_mm_set2_epi64(0x{:x}u);\n".format(2*TW) 
    for i in range(0,N) :
        str+="\tn[{}]=_mm_add_epi64(a[{}],b[{}]);\n".format(i,i,i)
    if not algorithm :
        str+=subp()
        str+="\tcarry=prop(n);\n"
        str+=caddp()
    str+="\t(void)prop(n);\n"    
    str+="}\n"
    return str

#modular subtraction
def modsub(n,m) :
    str="//Modular subtraction - reduce less than 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modsub{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)
    else :
        str+="void modsub{}(const spint *a,const spint *b,spint *n) {{\n".format(DECOR)
    if not algorithm :
        str+="\tspint carry;\n"
        str+="\tspint bot=_mm_set2_epi64({}u);\n".format(m*2) 
        str+="\tspint top=_mm_set2_epi64(0x{:x}u);\n".format(2*TW) 
    else :
        str+="\tspint bot=_mm_set2_epi64({}u);\n".format(m*mp)
        str+="\tspint top=_mm_set2_epi64(0x{:x}u);\n".format(mp*TW)
    for i in range(0,N) :
        str+="\tn[{}]=_mm_sub_epi64(a[{}],b[{}]);\n".format(i,i,i)
    if not algorithm :
        str+="\tcarry=prop(n);\n"
        str+=caddp()
    else :
        str+=addp()
    str+="\t(void)prop(n);\n" 
    str+="}\n"
    return str

#modular negation
def modneg(n,m) :
    str="//Modular negation\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modneg{}(const spint *b,spint *n) {{\n".format(DECOR)
    else :
        str+="void modneg{}(const spint *b,spint *n) {{\n".format(DECOR)
    if not algorithm :
        str+="\tspint carry;\n"
    str+="\tspint zero=_mm_set2_epi64(0);\n"
    if not algorithm :
        str+="\tspint bot=_mm_set2_epi64({}u);\n".format(m*2)
        str+="\tspint top=_mm_set2_epi64(0x{:x}u);\n".format(2*TW)       
    else :
        str+="\tspint bot=_mm_set2_epi64({}u);\n".format(m*mp)
        str+="\tspint top=_mm_set2_epi64(0x{:x}u);\n".format(mp*TW) 
    for i in range(0,N) :
        str+="\tn[{}]=_mm_sub_epi64(zero,b[{}]);\n".format(i,i)
    if not algorithm :
        str+="\tcarry=prop(n);\n"
        str+=caddp()
    else :
        str+=addp(mp)
    str+="\t(void)prop(n);\n" 
    str+="}\n"
    return str

#multiplication macro
def getZM(str,row,n,m) :
    N=getN(n)
    xcess=N*base-n
    mm=m*2**xcess
    k=row+1
    L=N-1

    first=True
    while k<N :
        if first :
            str+="\tmul(&ttl,&tth,a[{}],b[{}]);".format(k,L)    #"\ttt=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
            first=False
        else :
            str+=" accum(&ttl,&tth,a[{}],b[{}]);".format(k,L)     #" tt+=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
        L-=1
        k+=1
    if row<N-1:
        str+=" lo=_mm_and_si128(ttl,mask);"
        #str+=" lo=ttl & mask;"    #" lo=(spint)tt & mask;"
        if row==0 :
            str+=" accum(&tl,&th,lo,spm);"     #" t+=(dpint)lo*(dpint)0x{:x};".format(mm)
        else :
            str+=" accumx(&tl,&th,lo,hi,spm);"

        str+=" hi=shiftout(ttl,tth);"  #" hi=(spint)(tt>>{}u);".format(base)
    else :
        first=True
    
    k=0
    while k<=row :
        if first :
            str+="\taccum(&tl,&th,a[{}],b[{}]);".format(k,row-k)    #"\tt+=(dpint)a[{}]*(dpint)b[{}];".format(k,row-k)
            first=False
        else :
            str+=" accum(&tl,&th,a[{}],b[{}]);".format(k,row-k)   #" t+=(dpint)a[{}]*(dpint)b[{}];".format(k,row-k)
        k+=1
    if row==N-1 :
        str+=" accum(&tl,&th,hi,spm);"   #" t+=(dpint)hi*(dpint)0x{:x};".format(mm)
    str+=" spint v{}=_mm_and_si128(tl,mask); shiftr(&tl,&th);\n".format(row)    #" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
    return str

#squaring macro
def getZS(str,row,n,m) :
    N=getN(n)
    xcess=N*base-n
    mm=m*2**xcess
    k=row+1
    L=N-1
    first=True
    dble=False
    if k<L :
        dble=True
    while k<L :
        if first :
            str+="\tmul(&ttl,&tth,a[{}],a[{}]);".format(k,L)    #"\ttt=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
            first=False
        else :
            str+=" accum(&ttl,&tth,a[{}],a[{}]);".format(k,L)    #" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)

        L-=1
        k+=1
    if dble :
        str+=" ttl=_mm_add_epi64(ttl,ttl); tth=_mm_add_epi64(tth,tth);"
    if k==L :
        if first :
            str+="\tmul(&ttl,&tth,a[{}],a[{}]);".format(k,k)     #"\ttt=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
            first=False
        else :
            str+=" accum(&ttl,&tth,a[{}],a[{}]);".format(k,k)   #" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
    first=True
    if row<N-1:
        str+=" lo=_mm_and_si128(ttl,mask); "
        str+=" "
    else: 
        str+="\t"

    first=True
    k=0
    L=row
    
    dble=False
    if k<L :
        dble=True

    while k<L :
        if first :
            str+="mul(&t2l,&t2h,a[{}],a[{}]);".format(k,L)    #"t2=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
            first=False
        else :
            str+=" accum(&t2l,&t2h,a[{}],a[{}]);".format(k,L)   #" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
        k+=1
        L-=1

    if dble :
        str+=" t2l=_mm_add_epi64(t2l,t2l); t2h=_mm_add_epi64(t2h,t2h);"
    if k==L :
        if first :
            str+="mul(&t2l,&t2h,a[{}],a[{}]);".format(k,k)   #"t2=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
            first=False
        else :
            str+=" accum(&t2l,&t2h,a[{}],a[{}]);".format(k,k)  #" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
 
    if row==N-1 : 
        str+=" accum(&tl,&th,hi,spm);"   #" t+=(udpint)hi*(udpint)0x{:x};".format(mm)
    else :
        if row==0 :
            str+=" accum(&t2l,&t2h,lo,spm);"     #" t2+=(udpint)lo*(udpint)0x{:x};".format(mm)
        else :
            str+=" accumx(&t2l,&t2h,lo,hi,spm);"

        str+=" hi=shiftout(ttl,tth);"  #" hi=(spint)(tt>>{}u);".format(base)
    str+=" tl=_mm_add_epi64(tl,t2l); th=_mm_add_epi64(th,t2h); "
    str+=" spint v{}=_mm_and_si128(tl,mask); shiftr(&tl,&th);\n".format(row)     #" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
    return str

# second reduction pass
def second_pass(str,n,m) :
    N=getN(n)
    mask=(1<<base)-1
    xcess=N*base-n
# second reduction pass
    str+="// second reduction pass\n\n"  
    k=0
    if fred :
        str+="\tspint ut=tl;\n"  
        if xcess>0 :
            str+="\tspint smask=_mm_set2_epi64({});\n".format((1<<(base-xcess))-1)
            str+= "\tut=_mm_add_epi64(_mm_slli_epi64(ut,{}),_mm_srli_epi64(v{},{}); v{}=_mm_and_si128(v{},smask);\n".format(xcess,N-1,base-xcess,N-1,N-1)
            #str+= "\tut=(ut<<{})+(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)
        if m>1 :
            str+="\tut=muls(ut,0x{:x});\n".format(m)
            #str+="\tut*=0x{:x};\n".format(m)
        str+="\ts=_mm_add_epi64(v0,_mm_and_si128(ut,mask));\n"
        str+="\tc[0]=_mm_and_si128(s,mask);\n"
        str+= "\tcarry=_mm_add_epi64(_mm_srli_epi64(s,{}),_mm_srli_epi64(ut,{}));\n".format(base,base)
    else : 
        if xcess>0 :
            str+="\tspint smask=_mm_set2_epi64({});\n".format((1<<(base-xcess))-1)
            str+="\tshiftl(&tl,&th,{}); tl=_mm_add_epi64(tl,_mm_srli_epi64(v{},{}u)); v{}=_mm_and_si128(v{},smask);\n".format(xcess,N-1,base-xcess,N-1,N-1)   #"\tut=(ut<<{})+(spint)(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)
        if m>1 :
            str+= "\tmuli(&tl,&th,0x{:x});\n".format(m)  #"\tut*=0x{:x};\n".format(m)
        str+= "\ts=_mm_add_epi64(v0,_mm_and_si128(tl,mask));\n"
        str+= "\tc[0]=_mm_and_si128(s,mask);\n"
        str+= "\tcarry=_mm_add_epi64(_mm_srli_epi64(s,{}),shiftout(tl,th));\n".format(base)    #"\tcarry=(s>>{})+(spint)(ut>>{});\n".format(base,base)
    k=k+1

    str+= "\tc[{}]=_mm_add_epi64(v{},carry);\n".format(k,k)

    for i in range(k+1,N) :
        str+= "\tc[{}]=v{};\n".format(i,i)
    
    str+="\t(void)prop(c);\n"   # fully propagate carries
    str+="}\n"
    return str


# modular multiplication
# Note that allowing inlining gives significant speed-up
def modmul(n,m) :
    N=getN(n)
    mask=(1<<base)-1
    xcess=N*base-n
    mm=m*2**xcess
    str="// Modular multiplication, c=a*b mod 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modmul{}(const spint *a,const spint *b,spint *c) {{\n".format(DECOR)
    else :
        str+="void modmul{}(const spint *a,const spint *b,spint *c) {{\n".format(DECOR)

    str+="\tspint tl=_mm_set2_epi64(0);\n"
    str+="\tspint th=_mm_set2_epi64(0);\n"
    str+="\tspint spm=_mm_set2_epi64({});\n".format(mm)

    str+="\tspint ttl,tth;\n"
    str+="\tspint lo,hi;\n"
    #str+="\tuspint carry,s,mask=((uspint)1<<{}u)-(uspint)1;\n".format(base)
    
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"

    for row in range(0,N) :
        str=getZM(str,row,n,m)

    str=second_pass(str,n,m)

    return str

# modular squaring
# Note that allowing inlining gives significant speed-up
def modsqr(n,m) :
    N=getN(n)
    mask=(1<<base)-1
    xcess=N*base-n
    mm=m*2**xcess
    str="// Modular squaring, c=a*a mod 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modsqr{}(const spint *a,spint *c) {{\n".format(DECOR)
    else :
        str+="void modsqr{}(const spint *a,spint *c) {{\n".format(DECOR)
    
    str+="\tspint tl=_mm_set2_epi64(0);\n"
    str+="\tspint th=_mm_set2_epi64(0);\n"
    str+="\tspint spm=_mm_set2_epi64({});\n".format(mm)
    str+="\tspint ttl,tth;\n"
    str+="\tspint t2l,t2h;\n"
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"

    str+="\tspint lo,hi;\n"

    for row in range(0,N) :
        str=getZS(str,row,n,m)

    str=second_pass(str,n,m)

    return str

# multiply by an integer
def modmli(n,m) :
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
    str+="\tspint tl=_mm_set2_epi64(0);\n"
    str+="\tspint th=_mm_set2_epi64(0);\n"

    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint bw=_mm_set2_epi64(b);\n"
    str+="\tspint mask=_mm_set2_epi64(((int64_t)1<<52)-1);\n"

    for i in range(0,N) :
        str+="\taccum(&tl,&th,a[{}],bw);".format(i)
        str+="\tspint v{}=_mm_and_si128(tl,mask); shiftr(&tl,&th);\n".format(i)
    str=second_pass(str,n,m)

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
    str="//Calculate progenitor - use optimal addition chain\n"
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
    str+="}\n"
    subprocess.call("rm ac.txt",shell=True)    
    return str

#modular inversion
def modinv() :
    str="//calculate inverse, provide progenitor h if available, NULL if not\n"
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
        str+="\tfor (i=0;i<{};i++) {{\n".format(PM1D2-1)
        str+="\t\tmodsqr{}(s,s);\n".format(DECOR)
        str+="\t\tmodmul{}(s,x,s);\n".format(DECOR)
        str+="\t}\n"

    str+="\tmodnsqr{}(t,{});\n".format(DECOR,PM1D2+1)
    str+="\tmodmul{}(s,t,z);\n".format(DECOR)
    str+="}\n"
    return str

#test for quadratic residue
def modqr() :
    str="//Test for quadratic residue, provide progenitor h if available, NULL if not\n"
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
            str+=","
        str+=hex(ROI[N-1])
        str+="};\n"
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
    str+="\tspint d=_mm_set2_epi64(0);\n"
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
    str+="\tspint d=_mm_set2_epi64(0);\n"
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
    str+="\t\ta[i]=_mm_set2_epi64(0);\n"
    str+="\t}\n"
    str+="}\n"
    return str

#set to unity
def modone() :
    str="//set to one\n"
    if makestatic :
        str+="static "
    str+="void modone{}(spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\t\ta[0]=_mm_set2_epi64(1);\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=_mm_set2_epi64(0);\n"
    str+="\t}\n"
    str+="\tnres{}(a,a);\n".format(DECOR)
    str+="}\n"
    return str

#set to a small integer
def modint() :
    str="//set to integer\n"
    if makestatic :
        str+="static "
    str+="void modint{}(int x,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\ta[0]=_mm_set2_epi64(x);\n"  
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=_mm_set2_epi64(0);\n"
    str+="\t}\n"
    str+="\tnres{}(a,a);\n".format(DECOR)
    str+="}\n"
    return str

# convert to internal form
def nres(n) :
    str="//Convert m to internal form, n=nres(m) \n"
    if makestatic :
        str+="static "
    str+="void nres{}(const spint *m,spint *n) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\tn[i]=m[i];\n"
    str+="\t}\n"
    str+="}\n"
    return str 

#convert back to integer form (and do final subtract)
def redc(n) :
    str="//Convert n back to normal form, m=redc(n) \n"
    if makestatic :
        str+="static "
    str+="void redc{}(const spint *n,spint *m) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\tm[i]=n[i];\n"
    str+="\t}\n"
    str+="\t(void)modfsb{}(m);\n".format(DECOR)
    str+="}\n"
    return str 

#conditional swap -  see Loiseau et al. 2021
def modcsw() :
    str="//conditional swap g and f if d=1\n"
    str+="//strongly recommend inlining be disabled using compiler specific syntax\n"
    if makestatic :
        str+="static "
    str+="void modcsw{}(spint b,volatile spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c0,c1,s,t,w,v,aux;\n"
    str+="\tstatic uint64_t R0=0,R1=0;\n"
    str+="\tspint zero=_mm_set2_epi64(0);\n"
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
    return str

#conditional move
def modcmv() :
    str="//conditional move g to f if d=1\n"
    str+="//strongly recommend inlining be disabled using compiler specific syntax\n"
    if makestatic :
        str+="static "
    str+="void modcmv{}(spint b,const spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c0,c1,s,t,v,w,aux;\n"
    str+="\tstatic uint64_t R0=0,R1=0;\n"
    str+="\tspint zero=_mm_set2_epi64(0);\n"
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
    str+="\ta[n]=one; a[n]=_mm_slli_epi64(a[n],m);\n}\n"
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
    str+="\t\ta[i]=_mm_set2_epi64(0);\n\t}\n"
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
    print("\n//Automatically generated modular arithmetic C code for pseudo-Mersenne primes")
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
    print("#define MERSENNE")
    print("#define MULBYINT")
    if prime[0].isalpha() :
        print("#define",prime,"\n")

def functions() :
    print(intrinsics())
    print(prop(n))
    print(flat(n))
    print(modfsb(n))
    print(modadd(n,m))
    print(modsub(n,m))
    print(modneg(n,m))
    print(modmli(n,m))
    print(modmul(n,m))
    print(modsqr(n,m))
    print(modcpy())
    print(modnsqr())
    print(modpro())
    print(modinv())
    print(nres(n))
    print(redc(n))
    print(modis1(n))
    print(modis0(n,))
    print(modzer())
    print(modone())
    print(modint())
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
    print("Valid syntax - python pseudo_ifma.py <prime> OR <prime name>")
    print("For example - python pseudo_ifma.py X25519")
    print("For example - python pseudo_ifma.py 2**255-19")
    exit(2)

WL=64
prime=sys.argv[1]

n=0
base=52
m=0

# Note that number base is fixed as 2^52

### Start of user editable area

# More named primes can be added here
p=0
algorithm=False
mp=0

if prime=="PM266" :
    p=2**266-3

if prime=="NUMS256W" :
    p=2**256-189

if prime=="NUMS256E" :
    p=2**256-189

if prime=="NIST521" :
    p=2**521-1

if prime=="ED521" :
    p=2**521-1

if prime=="ED25519" :
    p=2**255-19
    if not generic :
        algorithm=True
        mp=4            # Assuming Edwards curve - see https://eprint.iacr.org/2017/437

if prime=="X25519" :
    p=2**255-19
    if not generic :
        algorithm=True  # if algorithm is known, fix multiple of prime mp (for modular subtractions) as described in https://eprint.iacr.org/2017/437
        mp=2            # Make sure there is sufficient excess - otherwise change default radix. Here assuming Montgomery ladder algorithm  RFC7748. 
                        # Now no reduction required after modular additions/subtractions.
    
if prime=="C2065" :
    p=2**206-5

if prime=="PM336" :
    p=2**336-3

if prime=="PM383" :
    p=2**383-187

if prime=="C41417" :
    p=2**414-17

if prime=="PM512" :
    p=2**512-569


if prime=="SECP256K1" :
    p=2**256-2**32-977

### End of user editable area

noname=False
if p==0 :
    if prime[0].isdigit() :
        p=eval(prime)
        noname=True    # unnamed prime
    else :
        print("This named modulus not supported (not a pseudo-Mersenne?)")
        exit(2)

n=p.bit_length() 
if n<120 or pow(3,p-1,p)!=1 :
    print("Not a sensible modulus, too small or not a prime")
    exit(2)

m=2**n-p
b=2**base

# find biggest power of 2 which divides p-1
p1=p-1                           
PM1D2=0
while (p1%2)==0 :
    PM1D2+=1
    p1>>=1
e=(1<<PM1D2)
PE=(p-1-e)//(2*e)                # exponent for use in inversion, QR check, and for square roots

if m>=b :
    print("Not an exploitable pseudo-Mersenne - Use Montgomery method instead")
    exit(1)

# get number of limbs
N=getN(n)
xcess=N*base-n
mm=m*2**xcess

if N>9 :
    inline=False

# check for excess too large
if mm >= 2**(WL-1) : #or m*((N+2)*2**(2*WL)) >= 2**(2*base+WL-1-xcess) :
    print("Unfortunate choice of radix - excess",xcess,"too large - try using smaller radix")
    exit(1)

if (n%base)==0 :
    TW=b
else :
    TW=(1<<(n%base)) # Top word of prime

#get number of bytes for export
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
print("Using alternate method")

# faster reduction
fred=False
if bits(N+1)+base+bits(mm)<WL :
    fred=True
    print("Tighter reduction")


# generate WL bit modular multiplication and squaring code for Mersenne prime 2^n-m
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
        print(modmul(n,m))
        print(modsqr(n,m))
        print(modcpy())
        print(modnsqr())
        print(modpro())
        print(modinv())
        print(nres(n))
        print(redc(n))
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
    #subprocess.call("rm time.c", shell=True)

# to determine code size
makestatic=False
with open('field.c', 'w') as f:
    with redirect_stdout(f):
        header()
        functions()
f.close()

subprocess.call(compiler+" -march=native -O3 -c field.c",shell=True)
subprocess.call("size field.o > size.txt",shell=True)

f=open('size.txt')
lines=f.readlines()
info=lines[1].split()

print("Code size using -O3 = ",info[0])

subprocess.call(compiler+" -march=native -Os -c field.c",shell=True)
subprocess.call("size field.o > size.txt",shell=True)

f=open('size.txt')
lines=f.readlines()
info=lines[1].split()

print("Code size using -Os = ",info[0])
subprocess.call("rm size.txt",shell=True)     
subprocess.call("rm field.o",shell=True)  

if decoration :
    if noname :
        DECOR="_"+str(n)+str(m)+"_ct"
    else :
        DECOR="_"+prime+"_ct"

makestatic=True
with open('field.c', 'w') as f:
    with redirect_stdout(f):
        header()
        functions()
f.close()

if formatted :
    subprocess.call("clang-format -i field.c", shell=True)  # tidy up the format

if check:
    subprocess.call("cppcheck --enable=all --addon=misc  --suppress=unusedFunction --suppress=missingIncludeSystem --suppress=checkersReport field.c", shell=True) 

print("Field code is in field.c")

sys.exit(base)
