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
# Execute this program as: python3 pseudo_ifma.py X25519
# Production code is output to file field.c
# 
# Mike Scott 24th November 2025
# TII
#

# Some default settings

NLANES=2 # Number of available SIMD 64-bit lanes in AVX-IFMA architecture (could be 2, 4 or 8)
if NLANES<2 :
    NLANES=2
if NLANES>8 :
    NLANES=8

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
scale=1 # set to 10 or 100 for faster timing loops. Default to 1
PSCR=False # Power Side Channel Resistant conditional moves and swaps

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

#conditional add of p
def caddp() :
    str=""
    str+="\tbot=MR_AND(bot,carry);\n"
    str+="\ttop=MR_AND(top,carry);\n"
    str+="\tn[0]=MR_SUB64U(n[0],bot);\n"
    str+="\tn[{}]=MR_ADD64U(n[{}],top);\n".format(N-1,N-1)
    return str

#add p
def addp() :
    str=""
    str+="\tn[0]=MR_SUB64U(n[0],bot);\n"
    str+="\tn[{}]=MR_ADD64U(n[{}],top);\n".format(N-1,N-1) 
    return str

#subtract p
def subp() :
    str=""
    str+="\tn[0]=MR_ADD64U(n[0],bot);\n"
    str+="\tn[{}]=MR_SUB64U(n[{}],top);\n".format(N-1,N-1)  
    return str

# Propagate carries. 
def prop(n) :
    str="//propagate carries\n"
    str+="static spint inline prop(spint *n) {\n"
    str+="\tint i;\n"
    str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"

    #str+="\tsspint carry=(sspint)n[0]>>{}u;\n".format(base)
    str+="\tsspint carry=(sspint)n[0];\n"
    str+="\tcarry=MR_SRA64S(carry,{}u);\n".format(base)

    str+="\tn[0]=MR_AND(n[0],mask);\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N-1)

    str+="\t\tcarry=MR_ADD64U(carry,n[i]);\n"
    str+="\t\tn[i]=MR_AND(carry,mask);\n"
    str+="\t\tcarry=MR_SRA64S(carry,{}u);\n".format(base)
    str+="\t}\n"
    str+="\tn[{}]=MR_ADD64U(n[{}],carry);\n".format(N-1,N-1)
    str+="\treturn (MR_SRA64S(n[{}],{}u));\n}}\n".format(N-1,WL-1)
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

    str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m) 
    str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(TW)  

    str+="\tspint carry=prop(n);\n"
    str+=caddp()
    str+="\t(void)prop(n);\n"
    str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
    str+="\treturn MR_AND(carry,mask);\n"
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
    str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m) 
    str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(TW) 
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
    str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m*2) 
    str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(2*TW) 
    for i in range(0,N) :
        str+="\tn[{}]=MR_ADD64U(a[{}],b[{}]);\n".format(i,i,i)
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
        str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m*2) 
        str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(2*TW) 
    else :
        str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m*mp)
        str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(mp*TW)
    for i in range(0,N) :
        str+="\tn[{}]=MR_SUB64U(a[{}],b[{}]);\n".format(i,i,i)
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
    str+="\tspint zero=MR_ZERO();\n"
    if not algorithm :
        str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m*2)
        str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(2*TW)       
    else :
        str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m*mp)
        str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(mp*TW) 
    for i in range(0,N) :
        str+="\tn[{}]=MR_SUB64U(zero,b[{}]);\n".format(i,i)
    if not algorithm :
        str+="\tcarry=prop(n);\n"
        str+=caddp()
    else :
        str+=addp()
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
            str+="\tttl=MR_ZERO();\n"
            str+="\ttth=MR_ZERO();\n"
            str+="\tttl=MR_MULADD_LO(ttl,a[{}],b[{}]); tth=MR_MULADD_HI(tth,a[{}],b[{}]);".format(k,L,k,L)
            #str+="\tmul(&ttl,&tth,a[{}],b[{}]);".format(k,L)    #"\ttt=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
            first=False
        else :
            str+=" ttl=MR_MULADD_LO(ttl,a[{}],b[{}]); tth=MR_MULADD_HI(tth,a[{}],b[{}]);".format(k,L,k,L)
            #str+=" accum(&ttl,&tth,a[{}],b[{}]);".format(k,L)     #" tt+=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
        L-=1
        k+=1
    if row<N-1:
        str+=" lo=MR_AND(ttl,mask);"
        #str+=" lo=ttl & mask;"    #" lo=(spint)tt & mask;"
        if row==0 :
            str+=" tl=MR_MULADD_LO_64(tl,lo,spm);"     #" t+=(dpint)lo*(dpint)0x{:x};".format(mm)
        else :
            str+=" tl=MR_MULADD_LO_64(tl,MR_ADD64U(lo,hi),spm);"

        str+=" hi=MR_SHR52(ttl,tth);"  #" hi=(spint)(tt>>{}u);".format(base)
    else :
        first=True
    
    k=0
    while k<=row :
        if first :
            str+="\ttl=MR_MULADD_LO(tl,a[{}],b[{}]); th=MR_MULADD_HI(th,a[{}],b[{}]);".format(k,row-k,k,row-k)
            #str+="\taccum(&tl,&th,a[{}],b[{}]);".format(k,row-k)    #"\tt+=(dpint)a[{}]*(dpint)b[{}];".format(k,row-k)
            first=False
        else :
            str+="tl=MR_MULADD_LO(tl,a[{}],b[{}]); th=MR_MULADD_HI(th,a[{}],b[{}]);".format(k,row-k,k,row-k) 
            #str+=" accum(&tl,&th,a[{}],b[{}]);".format(k,row-k)   #" t+=(dpint)a[{}]*(dpint)b[{}];".format(k,row-k)
        k+=1
    if row==N-1 :
        str+=" tl=MR_MULADD_LO_64(tl,hi,spm);"   #" t+=(dpint)hi*(dpint)0x{:x};".format(mm)
    str+=" spint v{}=MR_AND(tl,mask); tl=MR_SHR52(tl,th); th=MR_ZERO();\n".format(row)    #" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
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
            str+="\tttl=MR_ZERO(); "
            str+="tth=MR_ZERO(); "
            str+="ttl=MR_MULADD_LO(ttl,a[{}],a[{}]); tth=MR_MULADD_HI(tth,a[{}],a[{}]);".format(k,L,k,L)
            first=False
        else :
            str+=" ttl=MR_MULADD_LO(ttl,a[{}],a[{}]); tth=MR_MULADD_HI(tth,a[{}],a[{}]);".format(k,L,k,L)
            #str+=" accum(&ttl,&tth,a[{}],a[{}]);".format(k,L)    #" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)

        L-=1
        k+=1
    if dble :
        str+=" ttl=MR_ADD64U(ttl,ttl); tth=MR_ADD64U(tth,tth);"
    if k==L :
        if first :
            str+="\tttl=MR_ZERO(); "
            str+="tth=MR_ZERO(); "
            str+="ttl=MR_MULADD_LO(ttl,a[{}],a[{}]); tth=MR_MULADD_HI(tth,a[{}],a[{}]);".format(k,k,k,k)
            #str+="\tmul(&ttl,&tth,a[{}],a[{}]);".format(k,k)     #"\ttt=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
            first=False
        else :
            str+=" ttl=MR_MULADD_LO(ttl,a[{}],a[{}]); tth=MR_MULADD_HI(tth,a[{}],a[{}]); ".format(k,k,k,k) 
            #str+=" accum(&ttl,&tth,a[{}],a[{}]);".format(k,k)   #" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
    first=True
    if row<N-1:
        str+=" lo=MR_AND(ttl,mask); "
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
            str+="\tt2l=MR_ZERO(); "
            str+="t2h=MR_ZERO(); "
            str+=" t2l=MR_MULADD_LO(t2l,a[{}],a[{}]); t2h=MR_MULADD_HI(t2h,a[{}],a[{}]);".format(k,L,k,L)
            #str+="mul(&t2l,&t2h,a[{}],a[{}]);".format(k,L)    #"t2=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
            first=False
        else :
            str+=" t2l=MR_MULADD_LO(t2l,a[{}],a[{}]); t2h=MR_MULADD_HI(t2h,a[{}],a[{}]);".format(k,L,k,L)  
            #str+=" accum(&t2l,&t2h,a[{}],a[{}]);".format(k,L)   #" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
        k+=1
        L-=1

    if dble :
        str+=" t2l=MR_ADD64U(t2l,t2l); t2h=MR_ADD64U(t2h,t2h);"
    if k==L :
        if first :
            str+="\tt2l=MR_ZERO(); "
            str+="t2h=MR_ZERO(); "
            str+=" t2l=MR_MULADD_LO(t2l,a[{}],a[{}]); t2h=MR_MULADD_HI(t2h,a[{}],a[{}]);".format(k,k,k,k)
            #str+="mul(&t2l,&t2h,a[{}],a[{}]);".format(k,k)   #"t2=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
            first=False
        else :
            str+=" t2l=MR_MULADD_LO(t2l,a[{}],a[{}]); t2h=MR_MULADD_HI(t2h,a[{}],a[{}]);".format(k,k,k,k)
            #str+=" accum(&t2l,&t2h,a[{}],a[{}]);".format(k,k)  #" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
 
    if row==N-1 : 
        str+=" tl=MR_MULADD_LO_64(tl,hi,spm);"   #" t+=(udpint)hi*(udpint)0x{:x};".format(mm)
    else :
        if row==0 :
            str+=" t2l=MR_MULADD_LO_64(t2l,lo,spm);"     #" t2+=(udpint)lo*(udpint)0x{:x};".format(mm)
        else :
            str+=" t2l=MR_MULADD_LO_64(t2l,MR_ADD64U(lo,hi),spm);"

        str+=" hi=MR_SHR52(ttl,tth);"  #" hi=(spint)(tt>>{}u);".format(base)
    str+=" tl=MR_ADD64U(tl,t2l); th=MR_ADD64U(th,t2h); "
    str+=" spint v{}=MR_AND(tl,mask); tl=MR_SHR52(tl,th); th=MR_ZERO();\n".format(row)     #" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
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
            str+="\tspint smask=MR_SET_ALL_LANES_TO_CONSTANT({});\n".format((1<<(base-xcess))-1)
            str+= "\tut=MR_ADD64U(MR_SHL64U(ut,{}),MR_SHR64U(v{},{}); v{}=MR_AND(v{},smask);\n".format(xcess,N-1,base-xcess,N-1,N-1)
            #str+= "\tut=(ut<<{})+(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)
        if m>1 :
            str+="\tut=MR_MULLO(ut,MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}));\n".format(m)
            #str+="\tut*=0x{:x};\n".format(m)
        str+="\ts=MR_ADD64U(v0,MR_AND(ut,mask));\n"
        str+="\tc[0]=MR_AND(s,mask);\n"
        str+= "\tcarry=MR_ADD64U(MR_SHR64U(s,{}),MR_SHR64U(ut,{}));\n".format(base,base)
    else : 
        if xcess>0 :
            str+="\tspint smask=MR_SET_ALL_LANES_TO_CONSTANT({});\n".format((1<<(base-xcess))-1)
            str+="\ttl=MR_SHL64U(tl,{}); th=MR_SHL64U(th,{}); tl=MR_ADD64U(tl,MR_SHR64U(v{},{}u)); v{}=MR_AND(v{},smask);\n".format(xcess,xcess,N-1,base-xcess,N-1,N-1)   #"\tut=(ut<<{})+(spint)(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)
        if m>1 :
            str+= "\ttl=MR_MULLO(tl,MR_SET_ALL_LANES_TO_CONSTANT(0x{:x})); th=MR_MULLO(th,MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}));\n".format(m,m)
            #str+= "\tmuli(&tl,&th,0x{:x});\n".format(m)  #"\tut*=0x{:x};\n".format(m)
        str+= "\ts=MR_ADD64U(v0,MR_AND(tl,mask));\n"
        str+= "\tc[0]=MR_AND(s,mask);\n"
        str+= "\tcarry=MR_ADD64U(MR_SHR64U(s,{}),MR_SHR52(tl,th));\n".format(base)    #"\tcarry=(s>>{})+(spint)(ut>>{});\n".format(base,base)
    k=k+1

    str+= "\tc[{}]=MR_ADD64U(v{},carry);\n".format(k,k)

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

    str+="\tspint tl=MR_ZERO();\n"
    str+="\tspint th=MR_ZERO();\n"
    str+="\tspint spm=MR_SET_ALL_LANES_TO_CONSTANT({});\n".format(mm)

    str+="\tspint ttl,tth;\n"
    str+="\tspint lo,hi;\n"
    #str+="\tuspint carry,s,mask=((uspint)1<<{}u)-(uspint)1;\n".format(base)
    
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"

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
    
    str+="\tspint tl=MR_ZERO();\n"
    str+="\tspint th=MR_ZERO();\n"
    str+="\tspint spm=MR_SET_ALL_LANES_TO_CONSTANT({});\n".format(mm)
    str+="\tspint ttl,tth;\n"
    str+="\tspint t2l,t2h;\n"
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"

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
        str+="void inline modmli{}(const spint *a,spint bw,spint *c) {{\n".format(DECOR)
    else :
        str+="void modmli{}(const spint *a,spint bw,spint *c) {{\n".format(DECOR)
    str+="\tspint tl=MR_ZERO();\n"
    str+="\tspint th=MR_ZERO();\n"

    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    #str+="\tspint bw=MR_SET_ALL_LANES_TO_CONSTANT(b);\n"
    str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"

    for i in range(0,N) :
        str+="\ttl=MR_MULADD_LO(tl,a[{}],bw); th=MR_MULADD_HI(th,a[{}],bw);".format(i,i) 
        #str+="\taccum(&tl,&th,a[{}],bw);".format(i)
        str+="\tspint v{}=MR_AND(tl,mask); tl=MR_SHR52(tl,th); th=MR_ZERO();\n".format(i)
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
    f.close()
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

        str+="\t\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n";
        str+="\t\tspint d=MR_SUB64U(one,modis1{}(b));\n".format(DECOR)

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
    str+="\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
    str+="\tspint d=MR_ZERO();\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\td=MR_OR(d,c[i]);\n\t}\n"
    #str+="\tc0=(spint)c[0];\n"
    str+="\treturn MR_AND(MR_AND(one,MR_SHR64U(MR_SUB64U(d,one),{}u)),MR_SHR64U(MR_SUB64U(MR_XOR(c[0],one),one),{}u));\n}}\n".format(base,base)
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
    str+="\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
    str+="\tspint d=MR_ZERO();\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\td=MR_OR(d,c[i]);\n\t}\n" 
    str+="\treturn MR_AND(MR_SHR64U(MR_SUB64U(d,one),{}u),one);\n}}\n".format(base)
    return str

#set to zero
def modzer() :
    str="//set to zero\n"
    if makestatic :
        str+="static "
    str+="void modzer{}(spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=MR_ZERO();\n"
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
    str+="\t\ta[0]=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=MR_ZERO();\n"
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
    str+="\ta[0]=MR_SET_ALL_LANES_TO_CONSTANT(x);\n"  
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=MR_ZERO();\n"
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
        str+="\tstatic uint64_t R[2]={0,0};\n"
        str+="\tspint zero=MR_ZERO();\n"
        str+="\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
        str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"
        str+="\tR[0]+=0x3cc3c33c5aa5a55au;\n"
        str+="\tR[1]+=0x7447e88e1ee1e11eu;\n" 

        str+="\tw=MR_SET_EACH_LANE_TO_CONSTANT(R);\n"
        str+="\tc0=MR_ANDNOT(b,MR_ADD64U(w,one));\n"
        str+="\tc1=MR_ADD64U(b,w);\n" 
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\ts=g[i]; t=f[i];\n"
        str+="\t\tv=MR_MUL_DOUBLE_ADD(w,t,w,s);\n"
        #str+="\t\tv=w*(t+s);\n"
        str+="\t\tf[i]=aux=MR_MUL_DOUBLE_ADD(c0,t,c1,s);\n" 
        #str+="\t\tf[i] = aux = c0*t+c1*s;\n"
        str+="\t\tf[i]=MR_AND(MR_SUB64U(aux,v),mask);\n"
        #str+="\t\tf[i] = aux - v;\n"
        str+="\t\tg[i]=aux=MR_MUL_DOUBLE_ADD(c0,s,c1,t);\n"
        #str+="\t\tg[i] = aux = c0*s+c1*t;\n"
        str+="\t\tg[i]=MR_AND(MR_SUB64U(aux,v),mask);\n"
        #str+="\t\tg[i] = aux - v;\n\t}\n"
        str+="\t}\n"
        str+="}\n"
    else :
        str+="\tspint zero=MR_ZERO();\n"
        str+="\tspint delta,mask=MR_SUB64U(zero,b);\n"
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\tdelta=MR_AND(MR_XOR(g[i],f[i]),mask);\n"
        str+="\t\tg[i]=MR_XOR(g[i],delta);\n"
        str+="\t\tf[i]=MR_XOR(f[i],delta);\n\t}\n"
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
        str+="\tstatic uint64_t R[2]={0,0};\n"
        str+="\tspint zero=MR_ZERO();\n"
        str+="\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
        str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"
        str+="\tR[0]+=0x3cc3c33c5aa5a55au;\n"
        str+="\tR[1]+=0x7447e88e1ee1e11eu;\n" 

        str+="\tw=MR_SET_EACH_LANE_TO_CONSTANT(R);\n"
        str+="\tc0=MR_ANDNOT(b,MR_ADD64U(w,one));\n"
        str+="\tc1=MR_ADD64U(b,w);\n" 
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\ts=g[i]; t=f[i];\n"
        str+="\t\tv=MR_MUL_DOUBLE_ADD(w,t,w,s);\n"

        str+="\t\tf[i]=aux=MR_MUL_DOUBLE_ADD(c0,t,c1,s);\n"
        str+="\t\tf[i]=MR_AND(MR_SUB64U(aux,v),mask);\n"
        str+="\t}\n"
        str+="}\n"
    else :
        str+="\tspint zero=MR_ZERO();\n"
        str+="\tspint delta,mask=MR_SUB64U(zero,b);\n"
        str+="\tfor (i=0;i<{};i++) {{\n".format(N)
        str+="\t\tdelta=MR_AND(MR_XOR(g[i],f[i]),mask);\n"
        str+="\t\tf[i]=MR_XOR(f[i],delta);\n\t}\n"
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
    str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"
    str+="\ta[{}]=MR_ADD64U(MR_SHL64U(a[{}],n),MR_SHR64U(a[{}],{}u-n));\n".format(N-1,N-1,N-2,base)
    #str+="\ta[{}]=((a[{}]<<n)) | (a[{}]>>({}u-n));\n".format(N-1,N-1,N-2,base)
    str+="\tfor (i={};i>0;i--) {{\n".format(N-2)
    str+="\t\ta[i]=MR_ADD64U(MR_AND(MR_SHL64U(a[i],n),mask),MR_SHR64U(a[i-1],{}u-n));\n\t}}\n".format(base)
    #str+="\t\ta[i]=((a[i]<<n)&(spint)0x{:x}) | (a[i-1]>>({}u-n));\n\t}}\n".format(mask,base)
    str+="\ta[0]=MR_AND(MR_SHL64U(a[0],n),mask);\n"
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
    str+="\tspint mask=MR_SET_ALL_LANES_TO_CONSTANT(((int64_t)1<<52)-1);\n"
    str+="\tspint mskn=MR_SET_ALL_LANES_TO_CONSTANT((1<<n)-1);\n"    
    str+="\tspint r=MR_AND(a[0],mskn);\n"

    #str+="\tspint r=a[0]&(((spint)1<<n)-(spint)1);\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N-1)
    str+="\t\ta[i]=MR_ADD64U(MR_SHR64U(a[i],n),MR_AND(MR_SHL64U(a[i+1],{}u-n),mask));\n\t}}\n".format(base)
    #str+="\t\ta[i]=(a[i]>>n) | ((a[i+1]<<({}u-n))&(spint)0x{:x});\n\t}}\n".format(base,mask)
    str+="\ta[{}]=MR_SHR64U(a[{}],n);\n".format(N-1,N-1)
    #str+="\ta[{}]=a[{}]>>n;\n".format(N-1,N-1)
    str+="\treturn r;\n}\n"
    return str

#divide by 2
def modhaf(n) :
    N=getN(n)
    str="//divide by 2. Shift right 1 bit (or add p and shift right one bit)\n"
    if makestatic :
        str+="static "
    str+="void modhaf{}(spint *n) {{\n".format(DECOR)
    str+="\tspint lsb,one,d;\n"    
    str+="\tspint t[{}];\n".format(N)
    str+="\tspint bot=MR_SET_ALL_LANES_TO_CONSTANT({}u);\n".format(m) 
    str+="\tspint top=MR_SET_ALL_LANES_TO_CONSTANT(0x{:x}u);\n".format(TW)  
    str+="\t(void)prop(n);\n"
    str+="\tmodcpy{}(n,t);\n".format(DECOR)
    str+="\tlsb=modshr{}(1,t);\n".format(DECOR)
    str+=addp();
    str+="\t(void)prop(n);\n" 
    str+="\tmodshr{}(1,n);\n".format(DECOR)
    str+="\tone=MR_SET_ALL_LANES_TO_CONSTANT(1);\n";
    str+="\td=MR_SUB64U(one,lsb);\n".format(DECOR)
    str+="\tmodcmv{}(d,t,n);\n".format(DECOR)
    str+="}\n"
    return str

def mod2r() :
    str="//set a= 2^r\n"
    if makestatic :
        str+="static "
    str+="void mod2r{}(unsigned int r,spint *a) {{\n".format(DECOR)
    str+="\tunsigned int n=r/{}u;\n".format(base)
    str+="\tunsigned int m=r%{}u;\n".format(base)
    str+="\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
    str+="\tmodzer{}(a);\n".format(DECOR)
    str+="\tif (r>={}*8) return;\n".format(Nbytes)
    str+="\ta[n]=one; a[n]=MR_SHL64U(a[n],m);\n}\n"
    return str

#export to byte array
def modexp() :
    str="//export to byte array\n"
    if makestatic :
        str+="static "
    str+="void modexp{}(const spint *a,char *b[]) {{\n".format(DECOR)
    str+="\tint i,j;\n"
    str+="\tchar *ptr;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i={};i>=0;i--) {{\n".format(Nbytes-1)

    for j in range(0,NLANES) :
        str+="\t\tptr=b[{}];\n".format(j)
        str+="\t\tif (ptr!=NULL) ptr[i]=MR_EXTRACT(c[0],{});\n".format(j)
    str+="\t\t(void)modshr{}(8,c);\n\t}}\n".format(DECOR)
    str+="}\n"
    return str 

#import from byte array
def modimp() :
    str="//import from byte array\n"
    str+="//returns 1 if in range, else 0\n"
    if makestatic :
        str+="static "
    str+="spint modimp{}(const char *b[], spint *a) {{\n".format(DECOR)
    str+="\tint i,j;\n"
    str+="\tspint res;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=MR_ZERO();\n\t}\n"
    #str+="\t\ta[i]=0;\n\t}\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(Nbytes)
    str+="\t\tmodshl{}(8,a);\n".format(DECOR)
    str+="\t\tuint32_t bc[8];\n"
    str+="\t\tfor (j=0;j<{};j++) {{\n".format(NLANES)
    str+="\t\t\tconst char *ptr=b[j];\n"
    str+="\t\t\tif (ptr!=NULL) bc[j]=(unsigned char)ptr[i];\n"
    str+="\t\t\telse bc[j]=0;\n"
    str+="\t\t}\n"
    str+="\t\ta[0]=MR_ADD64U(a[0],MR_SET_EACH_LANE_TO_CONSTANT(bc));\n"
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
    str+="\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\treturn MR_AND(c[0],one);\n"
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
    str+="\tspint one=MR_SET_ALL_LANES_TO_CONSTANT(1);\n"
    str+="\tspint eq=one;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tredc{}(b,d);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\teq=MR_AND(eq, MR_AND(MR_SHR64U(MR_SUB64U(MR_XOR(c[i],d[i]),one),{}),one));\n\t}}\n".format(base)

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
        str+="x[{}]=MR_SET_ALL_LANES_TO_CONSTANT({}); ".format(i,hex(rap[i]))
        #str+="x[{}]={}; ".format(i,hex(rap[i]))
    str+="\n\t"
    for i in range(0,N) :
        str+="y[{}]=MR_SET_ALL_LANES_TO_CONSTANT({}); ".format(i,hex(rbp[i]))
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
        str+='\tprintf("modmul check 0x%06x Clock cycles= %d Nanosecs= %d\\n",MR_EXTRACT(z[0],0)&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
    else :
        str+='\tprintf("modmul check 0x%06x Nanosecs= %d\\n",MR_EXTRACT(z[0],0)&0xFFFFFF,elapsed);\n'
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
        str+="x[{}]=MR_SET_ALL_LANES_TO_CONSTANT({}); ".format(i,hex(rp[i]))
        #str+="x[{}]={}; ".format(i,hex(rp[i]))
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
        str+='\tprintf("modsqr check 0x%06x Clock cycles= %d Nanosecs= %d\\n",MR_EXTRACT(z[0],0)&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
    else :
        str+='\tprintf("modsqr check 0x%06x Nanosecs= %d\\n",MR_EXTRACT(z[0],0)&0xFFFFFF,elapsed);\n'
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
        str+="x[{}]=MR_SET_ALL_LANES_TO_CONSTANT({}); ".format(i,hex(rp[i]))
        #str+="x[{}]={}; ".format(i,hex(rp[i]))
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
        str+='\tprintf("modinv check 0x%06x Clock cycles= %d Nanosecs= %d\\n",MR_EXTRACT(z[0],0)&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000//scale)
    else :
        str+='\tprintf("modinv check 0x%06x Microsecs= %d\\n",MR_EXTRACT(z[0],0)&0xFFFFFF,elapsed);\n'
    str+="}\n"
    return str

def header() :
    print("\n//Automatically generated modular arithmetic C code for pseudo-Mersenne primes")
    print("//Command line : python {} {}".format(sys.argv[0], sys.argv[1]))
    print("//Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)\n")
    print("#include <stdio.h>")
    print("#include <stdint.h>\n")
    print("#define AVX128_IFMA 0")
    print("#define AVX256_IFMA 1")
    print("#define AVX512_IFMA 2")
    print("#include <emmintrin.h>")
    print("#include <smmintrin.h>\n")
    print("#include <immintrin.h>\n")

    if NLANES==2 :
        print("#define sspint __m128i")
        print("#define spint __m128i")
        print("#define SIMD_ENGINE AVX128_IFMA")
        print("#define SIMD_LANES 2\n")
    if NLANES==4 :
        print("#define sspint __m256i")
        print("#define spint __m256i")
        print("#define SIMD_ENGINE AVX256_IFMA")
        print("#define SIMD_LANES 4\n")
    if NLANES==8 :
        print("#define sspint __m512i")
        print("#define spint __m512i")
        print("#define SIMD_ENGINE AVX512_IFMA")
        print("#define SIMD_LANES 8\n")
    print('#include "ifma.h"')

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
    print(modhaf(n))
    print(mod2r())
    print(modexp())
    print(modimp())
    print(modsign())
    print(modcmp())


def main() :
    str="int main() {\n"
    str+='\tprintf("C Code generated from command line : python {} {}\\n");\n'.format(sys.argv[0], sys.argv[1])
    str+='\tprintf("Timing some functions - please wait\\n");\n'

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

print("Field code is in field.c")

sys.exit(base)
