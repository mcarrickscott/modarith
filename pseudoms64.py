# Python program to generate reasonably efficient C/C++ modular arithmetic code for pseudo-mersenne primes on a 16, 32 or 64-bit processor
# Modulus should be a pseudo-mersenne of the form 2^n-m
# uses unsaturated radix
# Special version for Microsoft 64-bit C compiler (no 128 bit integers!)
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

use_rdtsc=False # use rdtsc directly, x86 only, for better comparison with other implementations
decoration=False # decorate function names to avoid name clashes
formatted=True # pretty up the final output
inline=True # consider encouraging inlining
generic=True # set to False if algorithm is known in advance, in which case modadd and modsub can be faster - see https://eprint.iacr.org/2017/437. Set False for RFC7748 implementation.
allow_asr=True # Allow Arithmetic Shift Right. Maybe set to False to silence MISRA warnings
check=False # run cppcheck on the output
scale=1 # set to 10 or 100 for faster timing loops. Default to 1
fully_propagate=False # recommended set to True for Production code

import sys
import subprocess

# Determine optimal radix
def getbase(n) :
    limbs=int(n/WL)
    limit=2*WL
    while (True) :
        limbs=limbs+1
        if limbs==1 :   # must be at least 2 limbs
            limbs=2
        base=int(WL/2)
        while limbs*base<n :
            base=base+1
        if base>WL-3 :  # OK, base too large (moddadd/sub needs 3 bits), another limb required
            continue
        if limbs*(2**base-1)**2 < 2**limit :
            break
    return base


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
    str="#include <intrin.h>\n\n"

    str+="#define ARCH_X86_64 // remove for other architectures like ARM64 - Note umulh() intrinsic is still required \n"
    str+="#ifndef ARCH_X86_64\n"
    
    str+="// t+=a*b\n"
    str+="static inline void accum(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\tspint wl, wh;\n"
    str+="\twl=a*b;\n"
    str+="\twh=__umulh(a,b);\n"
    str+="\t*tl+=wl;\n"
    str+="\t*th+=(*tl<wl);\n"
    str+="\t*th+=wh;\n"
    str+="}\n\n"

    str+="// t+=(a+b)*c\n"
    str+="static inline void accumx(spint *tl,spint *th,spint a,spint bl,spint bh,spint c) {\n"
    str+="\tspint wl, wh, l, h;\n"
    str+="\tl=a;\n"
    str+="\tl+=bl;\n"
    str+="\th=bh+(l<bl);\n"
    str+="\twl=l*c;\n"
    str+="\twh=__umulh(l,c)+h*c;\n"
    str+="\t*tl+=wl;\n"
    str+="\t*th+=(*tl<wl);\n"
    str+="\t*th+=wh;\n"
    str+="}\n\n"

    str+="// t>>s\n"
    str+="static inline void shiftr(spint *tl,spint *th,char s){\n"
    str+="\t*tl=((*tl)>>s)|((*th)<<(64-s));\n"
    str+="\t*th=(*th)>>s;\n"
    str+="}\n\n"

    str+="// t<<s\n"
    str+="static inline void shiftl(spint *tl,spint *th,char s){\n"
    str+="\t*th = ((*th)<<s)|((*tl)>>(64-s));\n"
    str+="\t*tl = (*tl) << s;\n"
    str+="}\n\n"

    str+="// t+= (v<<s)\n"
    str+="static inline void accumsl(spint *tl,spint *th,spint v,char s) {\n"
    str+="\tspint wl, wh;\n"
    str+="\twh = v>>(64-s);\n"
    str+="\twl = v << s;\n"      
    str+="\t*tl+=wl;\n"
    str+="\t*th+=(*tl<wl);\n"
    str+="\t*th+=wh;\n"
    str+="}\n\n"

    str+="// r=t>>s\n"
    str+="static inline spint shiftout(spint tl,spint th,char s) {\n"
    str+="\treturn (tl>>s)|(th<<(64-s));\n"
    str+="}\n\n"

    str+="// t+=s\n"
    str+="static inline void add(spint *tl,spint *th,spint sl,spint sh) {\n"
    str+="\t*tl+=sl;\n"
    str+="\t*th+=(*tl<sl);\n"
    str+="\t*th+=sh;\n"
    str+="}\n\n"

    str+="// t-=s\n"
    str+="static inline void sub(spint *tl,spint *th,spint sl,spint sh) {\n"
    str+="\tpint l=*tl;\n"
    str+="\t*tl-=sl;\n"
    str+="\t*th-=((*tl)>l);\n"
    str+="\t*th-=sh;\n"
    str+="}\n\n"

    str+="// t=a*b\n"
    str+="static inline void mul(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\t*tl=a*b;\n"
    str+="\t*th=__umulh(a,b);\n"
    str+="}\n\n"

    str+="// t*=m\n"
    str+="static inline void muli(spint *tl,spint *th,spint m) {\n"
    str+="\tspint w=__umulh(*tl,m);\n"
    str+="\t*tl = *tl*m;\n"
    str+="\t*th = w + (*th) * m;\n"
    str+="}\n"

    str+="#else\n"

    str+="// t+=a*b\n"
    str+="static inline void accum(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\tunsigned char carry;\n"
    str+="\tspint wl,wh;\n"
    str+="\twl=_mulx_u64(a,b,&wh);\n"
    str+="\tcarry=_addcarryx_u64(0,wl,*tl,tl);\n"
    str+="\t_addcarryx_u64(carry,wh,*th,th);\n"
    str+="}\n\n"

    str+="// t+=(a+b)*c\n"
    str+="static inline void accumx(spint *tl,spint *th,spint a,spint bl,spint bh,spint c) {\n"
    str+="\tunsigned char carry;\n"
    str+="\tspint wl,wh,l,h;\n"
    str+="\tcarry=_addcarry_u64(0,a,bl,&l);\n"
    str+="\t_addcarry_u64(carry,0,bh,&h);\n"
    str+="\twl=_mulx_u64(l,c,&wh);\n"
    str+="\twh+=h*c;\n"
    str+="\tcarry=_addcarryx_u64(0,wl,*tl,tl);\n" 
    str+="\t_addcarryx_u64(carry,wh,*th,th);\n"
    str+="}\n\n"

    str+="// t>>s\n"
    str+="static inline void shiftr(spint *tl,spint *th,char s){\n"
    str+="\t*tl=__shiftright128(*tl,*th,s);\n"
    str+="\t*th=(*th>>s);\n"
    str+="}\n\n"

    str+="// t<<s\n"
    str+="static inline void shiftl(spint *tl,spint *th,char s){\n"
    str+="\t*th=__shiftleft128(*tl,*th,s);\n"
    str+="\t*tl=(*tl<<s);\n"
    str+="}\n\n"

    str+="// t+= (v<<s)\n"
    str+="static inline void accumsl(spint *tl,spint *th,spint v,char s) {\n"
    str+="\tunsigned char carry;\n"
    str+="\tspint l,h;\n"
    str+="\th=__shiftleft128(v,0,s);\n"
    str+="\tl=(v<<s);\n"
    str+="\tcarry=_addcarry_u64(0,l,*tl,tl);\n"
    str+="\t_addcarry_u64(carry,h,*th,th);\n"
    str+="}\n\n"

    str+="// r=t>>s\n"
    str+="static inline spint shiftout(spint tl,spint th,char s) {\n"
    str+="\treturn __shiftright128(tl,th,s);\n"
    str+="}\n\n"

    str+="// t+=s\n"
    str+="static inline void add(spint *tl,spint *th,spint sl,spint sh) {\n"
    str+="\tunsigned char carry;\n"
    str+="\tcarry=_addcarryx_u64(0,*tl,sl,tl);\n"
    str+="\t_addcarryx_u64(carry,*th,sh,th);\n"
    str+="}\n\n"

    str+="// t-=s\n"
    str+="static inline void sub(spint *tl,spint *th,spint sl,spint sh) {\n"
    str+="\tunsigned char carry;\n"
    str+="\tcarry=_subborrow_u64(0,*tl,sl,tl);\n"
    str+="\t_subborrow_u64(carry,*th,sh,th);\n"
    str+="}\n\n"

    str+="// t=a*b\n"
    str+="static inline void mul(spint *tl,spint *th,spint a,spint b) {\n"
    str+="\t*tl=_mulx_u64(a,b,th);\n"
    str+="}\n\n"
    
    str+="// t*=m\n"
    str+="static inline void muli(spint *tl,spint *th,spint m) {\n"
    str+="\tunsigned char carry;\n"
    str+="\tspint w;\n"
    str+="\t*tl=_mulx_u64(*tl,m,&w);\n"
    str+="\t*th=w+(*th)*m;\n"
    str+="}\n\n"

    str+="#endif\n\n"
    return str

#conditional add of p
def caddp(x) :
    str=""
    str+="\tn[0]-=((spint){}u)&carry;\n".format(m*x)
    str+="\tn[{}]+=((spint)0x{:x}u)&carry;\n ".format(N-1,x*TW)
    return str

#add p
def addp(x) :
    str=""
    str+="\tn[0]-=(spint){}u;\n".format(m*x)
    str+="\tn[{}]+=(spint)0x{:x}u;\n".format(N-1,x*TW)
    return str

#subtract p
def subp(x) :
    str=""
    str+="\tn[0]+=(spint){}u;\n".format(m*x)
    str+="\tn[{}]-=(spint)0x{:x}u;\n".format(N-1,x*TW)
    return str

# Propagate carries. Modified to please MISRA requirements and clang bug
def prop(n) :
    str="//propagate carries\n"
    str+="static spint inline prop(spint *n) {\n"
    str+="\tint i;\n"
    str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)
    if not allow_asr :
        str+="\tspint cst=0x{:x}u;\n".format(((1<<base)-1)<<(WL-base))
        str+="\tspint carry=n[0]>>{}u;\n".format(base)
        str+="\tcarry+=((-(carry>>{}u))&cst);\n".format(WL-base-1)
    else :
        #str+="\tsspint carry=(sspint)n[0]>>{}u;\n".format(base)
        str+="\tsspint carry=(sspint)n[0];\n"
        str+="\tcarry>>={}u;\n".format(base)

    str+="\tn[0]&=mask;\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N-1)
    if not allow_asr :
        str+="\t\tcarry += n[i];\n"
        str+="\t\tn[i] = carry & mask;\n"
        str+="\t\tcarry>>={}u;\n".format(base)
        str+="\t\tcarry+=((-(carry>>{}u))&cst);\n \n".format(WL-base-1)
    else :
        str+="\t\tcarry+=(sspint)n[i];\n"
        str+="\t\tn[i] = (spint)carry & mask;\n"
        str+="\t\tcarry>>={}u;\n".format(base)
    str+="\t}\n"
    str+="\tn[{}]+=(spint)carry;\n".format(N-1)
    str+="\treturn -((n[{}]>>1)>>{}u);\n}}\n".format(N-1,WL-2)
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
    str+="\tspint carry=prop(n);\n"
    str+=caddp(1)
    str+="\t(void)prop(n);\n"
    str+="\treturn (carry&1);\n"
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
    str+=subp(1)
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
    for i in range(0,N) :
        str+="\tn[{}]=a[{}]+b[{}];\n".format(i,i,i)
    if not algorithm :
        str+=subp(2)
        str+="\tcarry=prop(n);\n"
        str+=caddp(2)
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
    for i in range(0,N) :
        str+="\tn[{}]=a[{}]-b[{}];\n".format(i,i,i)
    if not algorithm :
        str+="\tcarry=prop(n);\n"
        str+=caddp(2)
    else :
        str+=addp(mp)
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
    for i in range(0,N) :
        str+="\tn[{}]=(spint)0-b[{}];\n".format(i,i)
    if not algorithm :
        str+="\tcarry=prop(n);\n"
        str+=caddp(2)
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
        if EPM :
            if first :
                str+="\t"
            else :
                str+=" "
            str+="accum(&tl,&th,ma{},b[{}]);".format(k,L)   #"t+=(dpint)ma{}*(dpint)b[{}];".format(k,L)
        else :
            if first :
                str+="\tmul(&ttl,&tth,a[{}],b[{}]);".format(k,L)    #"\ttt=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
                first=False
            else :
                str+=" accum(&ttl,&tth,a[{}],b[{}]);".format(k,L)     #" tt+=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
        L-=1
        k+=1
    if row<N-1:
        if overflow :
            str+=" lo=ttl & mask;"    #" lo=(spint)tt & mask;"
            if row==0 :
                str+=" accum(&tl,&th,lo,0x{:x});".format(mm)     #" t+=(dpint)lo*(dpint)0x{:x};".format(mm)
            else :
                if bad_overflow :
                    str+=" accumx(&tl,&th,lo,hil,hih,0x{:x});".format(mm)   #" t+=(hi+(dpint)lo)*(dpint)0x{:x};".format(mm)
                else :
                    str+=" accum(&tl,&th,lo+hi,0x{:x});".format(mm)   #" t+=(dpint)(spint)(lo+hi)*(dpint)0x{:x};".format(mm)
            if bad_overflow :
                str+=" hil=ttl; hih=tth; shiftr(&hil,&hih,{}u);".format(base)    #" hi=tt>>{}u;".format(base)
            else :
                str+=" hi=shiftout(ttl,tth,{}u);".format(base)  #" hi=(spint)(tt>>{}u);".format(base)
        else :
            if not EPM :
                str+=" muli(&ttl,&tth,0x{:x});".format(mm)  #" tt*=0x{:x};".format(mm)
                str+=" add(&tl,&th,ttl,tth);"   #" t+=tt;"
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
    if row==N-1 and overflow :
        str+=" accum(&tl,&th,hi,0x{:x});".format(mm)   #" t+=(dpint)hi*(dpint)0x{:x};".format(mm)
    str+=" spint v{}=tl & mask; shiftr(&tl,&th,{}u);\n".format(row,base)    #" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
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
        if EPM :
            if dble :   
                if first :
                    str+="\t"
                else :
                    str+=" "                
                str+="accum(&tl,&th,ma{},ta{});".format(k,L)    #"t+=(udpint)ma{}*(udpint)ta{};".format(k,L)
            else :
                if first :
                    str+="\taccum(&tl,&th,ma{},a[{}]);".format(k,L)   #"\tt+=(udpint)ma{}*(udpint)a[{}];".format(k,L)
                    first=False
                else :
                    str+=" accum(&tl,&th,ma{},a[{}]);".format(k,L)   #" t+=(udpint)ma{}*(udpint)a[{}];".format(k,L)
        else :
            if first :
                str+="\tmul(&ttl,&tth,a[{}],a[{}]);".format(k,L)    #"\ttt=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
                first=False
            else :
                str+=" accum(&ttl,&tth,a[{}],a[{}]);".format(k,L)    #" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)

        L-=1
        k+=1
    if dble :
        if not EPM :
            str+=" add(&ttl,&tth,ttl,tth);"      #" tt*=2;"
    if k==L :
        if EPM :
            if first :
                str+="\t"
            else :
                str+=" "  
            str+="accum(&tl,&th,ma{},a[{}]);".format(k,k)    #"t+=(udpint)ma{}*(udpint)a[{}];".format(k,k)
        else :
            if first :
                str+="\tmul(&ttl,&tth,a[{}],a[{}]);".format(k,k)     #"\ttt=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
                first=False
            else :
                str+=" accum(&ttl,&tth,a[{}],a[{}]);".format(k,k)   #" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
    first=True
    if row<N-1:
        if overflow :
            str+=" lo=tl & mask;"
        else :
            if not EPM :
                str+=" muli(&ttl,&tth,0x{:x});".format(mm)  #" tt*=0x{:x};".format(mm)
                str+=" add(&tl,&th,ttl,tth);"   #" t+=tt;"
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
        if EPM and dble :
            str+="accum(&tl,&th,a[{}],ta{});".format(k,L)   #"t+=(udpint)a[{}]*(udpint)ta{};".format(k,L)
        else :
            if first :
                str+="mul(&t2l,&t2h,a[{}],a[{}]);".format(k,L)    #"t2=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
                first=False
            else :
                str+=" accum(&t2l,&t2h,a[{}],a[{}]);".format(k,L)   #" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
        k+=1
        L-=1

    if dble :
        if not EPM :
            str+=" add(&t2l,&t2h,t2l,t2h);"  #" t2*=2;"
    if k==L :
        if EPM :
            str+=" accum(&tl,&th,a[{}],a[{}]);".format(k,k)   #" t+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
        else :
            if first :
                str+="mul(&t2l,&t2h,a[{}],a[{}]);".format(k,k)   #"t2=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
                first=False
            else :
                str+=" accum(&t2l,&t2h,a[{}],a[{}]);".format(k,k)  #" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
 

    if overflow :
        if row==N-1 : 
            str+=" accum(&tl,&th,hi,0x{:x});".format(mm)   #" t+=(udpint)hi*(udpint)0x{:x};".format(mm)
        else :
            if row==0 :
                str+=" accum(&t2l,&t2h,lo,0x{:x});".format(mm)     #" t2+=(udpint)lo*(udpint)0x{:x};".format(mm)
            else :
                if bad_overflow :
                    str+=" accumx(&t2l,&t2h,lo,hil,hih,0x{:x});".format(mm)    #" t2+=(hi+(udpint)lo)*(udpint)0x{:x};".format(mm)
                else :
                    str+=" accum(&t2l,&t2h,lo+hi,0x{:x});".format(mm)    #" t2+=(udpint)(spint)(lo+hi)*(udpint)0x{:x};".format(mm)
            if bad_overflow :
                str+=" hil=ttl; hih=tth; shiftr(&hil,&hih,{}u);".format(base)  #" hi=tt>>{}u;".format(base)
            else :
                str+=" hi=shiftout(ttl,tth,{}u);".format(base)  #" hi=(spint)(tt>>{}u);".format(base)
    if not EPM :
        str+=" add(&tl,&th,t2l,t2h);"    #" t+=t2;"
    str+=" spint v{}=tl&mask; shiftr(&tl,&th,{}u);\n".format(row,base)     #" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
    return str

# second reduction pass
def second_pass(str,n,m) :
    N=getN(n)
    mask=(1<<base)-1
    xcess=N*base-n
# second reduction pass
    str+="// second reduction pass\n\n"  
    if fred :
        str+="\tspint ut=tl;\n"  
        if xcess>0 :
            smask=(1<<(base-xcess))-1
            str+= "\tut=(ut<<{})+(spint)(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)

        k=0
        if m>1 :
            str+= "\tut*=0x{:x};\n".format(m)
        str+= "\ts=v0+(ut & mask);\n"
        str+= "\tc[0]=(s&mask);\n"

        if carry_on :
            str+="\tut=(s>>{})+(ut>>{});\n".format(base,base)
            str+="\ts=v1+(ut & mask);\n"
            str+= "\tc[1]=(s&mask);\n"
            k+=1

        str+= "\tcarry=(s>>{})+(ut>>{});\n".format(base,base)
    else : 
        if xcess>0 :
            smask=(1<<(base-xcess))-1
            str+= "\tshiftl(&tl,&th,{}); add(&tl,&th,(v{}>>{}u),0); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)   #"\tut=(ut<<{})+(spint)(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)

        k=0
        if m>1 :
            str+= "\tmuli(&tl,&th,0x{:x});\n".format(m)  #"\tut*=0x{:x};\n".format(m)
        str+= "\ts=v0+(tl & mask);\n"
        str+= "\tc[0]=(spint)(s&mask);\n"

        if carry_on :
            str+="\tshiftr(&tl,&th,{}); add(&tl,&th,s>>{},0);\n".format(base,base)       #"\tut=(udpint)(s>>{})+(ut>>{});\n".format(base,base)
            str+="\ts=v1+(tl & mask);\n"
            str+= "\tc[1]=(spint)(s&mask);\n"
            k+=1

        str+= "\tcarry=(s>>{})+shiftout(tl,th,{});\n".format(base,base)    #"\tcarry=(s>>{})+(spint)(ut>>{});\n".format(base,base)
    k=k+1


    str+= "\tc[{}]=v{}+carry;\n".format(k,k)

    for i in range(k+1,N) :
        str+= "\tc[{}]=v{};\n".format(i,i)
    if fully_propagate :
        str+="\t(void)prop(c);\n"
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

    str+="\tspint tl=0, th=0;\n"

    if  EPM  :
        for i in range(1,N) :
            str+="\tspint ma{}=a[{}]*(spint)0x{:x};\n".format(i,i,mm)
    else :
        str+="\tspint ttl,tth;\n"

    if overflow :
        str+="\tspint lo;\n"
        if bad_overflow :
            str+="\tspint hil,hih;\n"   # could overflow single type
        else :
            str+="\tspint hi;\n"
    #str+="\tuspint carry,s,mask=((uspint)1<<{}u)-(uspint)1;\n".format(base)
    
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

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
    str="// Modular squaring, c=a*a mod 2p\n"
    if makestatic :
        str+="static "
    if inline and makestatic:
        str+="void inline modsqr{}(const spint *a,spint *c) {{\n".format(DECOR)
    else :
        str+="void modsqr{}(const spint *a,spint *c) {{\n".format(DECOR)
    
    str+="\tspint tl=0,th=0;\n"

    if  EPM  :
        for i in range(1,N) :
            str+="\tspint ta{}=a[{}]*(spint)2;\n".format(i,i)
        for i in range(1,N) :
            str+="\tspint ma{}=a[{}]*(spint)0x{:x};\n".format(i,i,mm)
    else :
        str+="\tspint ttl,tth;\n"
        str+="\tspint t2l,t2h;\n"
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

    if overflow :
        str+="\tspint lo;\n"
        if bad_overflow :
            str+="\tspint hil,hih;\n"
        else :
            str+="\tspint hi;\n"

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
    str+="\tspint tl=0,th=0;\n"

    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

    for i in range(0,N) :
        str+="\taccum(&tl,&th,a[{}],b);".format(i)
        str+="\tspint v{}=tl & mask; shiftr(&tl,&th,{}u);\n".format(i,base)

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
    str+="int modqr{}(const spint *h,const spint *x) {{\n".format(DECOR)
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
        str+="\t\tint d=1-modis1{}(b);\n".format(DECOR)
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
    str+="int modis1{}(const spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tspint c0;\n"
    str+="\tspint d=0;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\td|=c[i];\n\t}\n"
    str+="\tc0=(spint)c[0];\n"
    str+="\treturn ((spint)1 & ((d-(spint)1)>>{}u) & (((c0^(spint)1)-(spint)1)>>{}u));\n}}\n".format(base,base)
    return str

#test for zero
def modis0(n) :
    str="//is zero?\n"
    if makestatic :
        str+="static "
    str+="int modis0{}(const spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tspint d=0;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\td|=c[i];\n\t}\n"
    str+="\treturn ((spint)1 & ((d-(spint)1)>>{}u));\n}}\n".format(base)
    return str

#set to zero
def modzer() :
    str="//set to zero\n"
    if makestatic :
        str+="static "
    str+="void modzer{}(spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=0;\n"
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
    str+="\ta[0]=1;\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=0;\n"
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
    str+="\ta[0]=(spint)x;\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=0;\n"
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
    str+="void modcsw{}(int b,volatile spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c0,c1,s,t,w,v,aux;\n"
    str+="\tstatic spint R=0;\n"
    str+="\tR+=0x3cc3c33c5aa5a55au;\n"
    str+="\tw=R;\n"
    str+="\t\tc0=(~b)&(w+1);\n"
    str+="\t\tc1=b+w;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ts=g[i]; t=f[i];\n"
    str+="\t\tv=w*(t+s);\n"
    str+="\t\tf[i] = aux = c0*t+c1*s;\n"
    str+="\t\tf[i] = aux - v;\n"
    str+="\t\tg[i] = aux = c0*s+c1*t;\n"
    str+="\t\tg[i] = aux - v;\n\t}\n"
    str+="}\n"
    return str

#conditional move
def modcmv() :
    str="//conditional move g to f if d=1\n"
    str+="//strongly recommend inlining be disabled using compiler specific syntax\n"
    if makestatic :
        str+="static "
    str+="void modcmv{}(int b,const spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c0,c1,s,t,w,aux;\n"
    str+="\tstatic spint R=0;\n"
    str+="\tR+=0x3cc3c33c5aa5a55au;\n"
    str+="\tw=R;\n"
    str+="\t\tc0=(~b)&(w+1);\n"
    str+="\t\tc1=b+w;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ts=g[i]; t=f[i];\n"
    str+="\t\tf[i] = aux = c0*t+c1*s;\n"
    str+="\t\tf[i] = aux - w*(t+s);\n\t}\n"
    str+="}\n"
    return str


#shift left
def modshl(n) :
    N=getN(n)
    mask=(1<<base)-1
    str="//shift left by less than a word\n"
    if makestatic :
        str+="static "
    str+="void modshl{}(unsigned int n,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\ta[{}]=((a[{}]<<n)) | (a[{}]>>({}u-n));\n".format(N-1,N-1,N-2,base)
    str+="\tfor (i={};i>0;i--) {{\n".format(N-2)
    str+="\t\ta[i]=((a[i]<<n)&(spint)0x{:x}) | (a[i-1]>>({}u-n));\n\t}}\n".format(mask,base)
    str+="\ta[0]=(a[0]<<n)&(spint)0x{:x};\n".format(mask)
    str+="}\n"
    return str 

#shift right
def modshr(n) :
    N=getN(n)
    mask=(1<<base)-1
    str="//shift right by less than a word. Return shifted out part\n"
    if makestatic :
        str+="static "
    str+="int modshr{}(unsigned int n,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint r=a[0]&(((spint)1<<n)-(spint)1);\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N-1)
    str+="\t\ta[i]=(a[i]>>n) | ((a[i+1]<<({}u-n))&(spint)0x{:x});\n\t}}\n".format(base,mask)
    str+="\ta[{}]=a[{}]>>n;\n".format(N-1,N-1)
    str+="\treturn r;\n}\n"
    return str

def mod2r() :
    str="//set a= 2^r\n"
    if makestatic :
        str+="static "
    str+="void mod2r{}(unsigned int r,spint *a) {{\n".format(DECOR)
    str+="\tunsigned int n=r/{}u;\n".format(base)
    str+="\tunsigned int m=r%{}u;\n".format(base)
    str+="\tmodzer{}(a);\n".format(DECOR)
    str+="\tif (r>={}*8) return;\n".format(Nbytes)
    str+="\ta[n]=1; a[n]<<=m;\n}\n"
    return str

#export to byte array
def modexp() :
    str="//export to byte array\n"
    if makestatic :
        str+="static "
    str+="void modexp{}(const spint *a,char *b) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i={};i>=0;i--) {{\n".format(Nbytes-1)
    str+="\t\tb[i]=c[0]&(spint)0xff;\n"
    str+="\t\t(void)modshr{}(8,c);\n\t}}\n".format(DECOR)
    str+="}\n"
    return str 

#import from byte array
def modimp() :
    str="//import from byte array\n"
    str+="//returns 1 if in range, else 0\n"
    if makestatic :
        str+="static "
    str+="int modimp{}(const char *b, spint *a) {{\n".format(DECOR)
    str+="\tint i,res;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=0;\n\t}\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(Nbytes)
    str+="\t\tmodshl{}(8,a);\n".format(DECOR)
    str+="\t\ta[0]+=(spint)(unsigned char)b[i];\n\t}\n"
    str+="\tres=(int)modfsb{}(a);\n".format(DECOR)
    str+="\tnres{}(a,a);\n".format(DECOR)
    str+="\treturn res;\n"
    str+="}\n"
    return str 

#get sign (parity of value)
def modsign() :
    str="//determine sign\n"
    if makestatic :
        str+="static "
    str+="int modsign{}(const spint *a) {{\n".format(DECOR)
    str+="\tspint c[{}];\n".format(N)
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\treturn c[0]%2;\n"
    str+="}\n"
    return str

#compare for equality
def modcmp() :
    str="//return true if equal\n"
    if makestatic :
        str+="static "
    str+="int modcmp{}(const spint *a,const spint *b) {{\n".format(DECOR)
    str+="\tspint c[{}],d[{}];\n".format(N,N)
    str+="\tint i,eq=1;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tredc{}(b,d);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\teq&=(((c[i]^d[i])-1)>>{})&1;\n\t}}\n".format(base)
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
        str+="x[{}]={}; ".format(i,hex(rap[i]))
    str+="\n\t"
    for i in range(0,N) :
        str+="y[{}]={}; ".format(i,hex(rbp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    str+="\tnres{}(y,y);\n".format(DECOR)

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
    if use_rdtsc :
        str+="\tfinish=__rdtsc();\n"
    str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10*scale)
    str+="\tredc{}(z,z);\n".format(DECOR)
 
    if use_rdtsc :
        str+='\tprintf("modmul check 0x%06x Clock cycles= %d Nanosecs= %d\\n",(int)z[0]&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
    else :
        str+='\tprintf("modmul check 0x%06x Nanosecs= %d\\n",(int)z[0]&0xFFFFFF,elapsed);\n'
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
        str+="x[{}]={}; ".format(i,hex(rp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    if use_rdtsc :
        str+="\tstart=__rdtsc();\n" 
    str+="\tbegin=clock();\n"
    str+="\tfor (i=0;i<{};i++)\n".format(100000//scale)
    str+="\t\tfor (j=0;j<500;j++) {\n"
    str+="\t\t\tmodsqr{}(x,z);\n".format(DECOR)
    str+="\t\t\tmodsqr{}(z,x);\n".format(DECOR)
    str+="\t\t}\n"
    if use_rdtsc :
        str+="\tfinish=__rdtsc();\n"
    str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10*scale)
    str+="\tredc{}(z,z);\n".format(DECOR)
 
    if use_rdtsc :
        str+='\tprintf("modsqr check 0x%06x Clock cycles= %d Nanosecs= %d\\n",(int)z[0]&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
    else :
        str+='\tprintf("modsqr check 0x%06x Nanosecs= %d\\n",(int)z[0]&0xFFFFFF,elapsed);\n'
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
        str+="x[{}]={}; ".format(i,hex(rp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    if use_rdtsc :
        str+="\tstart=__rdtsc();\n"
    str+="\tbegin=clock();\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(50000//scale)
    str+="\t\t\tmodinv{}(x,NULL,z);\n".format(DECOR)
    str+="\t\t\tmodinv{}(z,NULL,x);\n".format(DECOR)
    str+="\t}\n"
    if use_rdtsc :
        str+="\tfinish=__rdtsc();\n"
    str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10000*scale)
    str+="\tredc{}(z,z);\n".format(DECOR)
 
    if use_rdtsc:
        str+='\tprintf("modinv check 0x%06x Clock cycles= %d Nanosecs= %d\\n",(int)z[0]&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000//scale)
    else :
        str+='\tprintf("modinv check 0x%06x Microsecs= %d\\n",(int)z[0]&0xFFFFFF,elapsed);\n'
    str+="}\n"
    return str

def header() :
    print("\n//Automatically generated modular arithmetic C code for pseudo-Mersenne primes")
    print("//Command line : python {} {}".format(sys.argv[0], sys.argv[1]))
    print("//Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)\n")
    print("#include <stdio.h>")
    print("#include <stdint.h>\n")
    print("#define sspint int{}_t".format(WL))
    print("#define spint uint{}_t".format(WL))
 
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
    print("Valid syntax - python pseudoms64.py <prime> OR <prime name>")
    print("For example - python pseudoms64.py X25519")
    print("For example - python pseudoms64.py 2**255-19")
    exit(2)

WL=64
prime=sys.argv[1]

n=0
base=0
m=0

# Note that number base MUST be less than wordlength WL (unsaturated).
# The radix for an n-bit prime is 2^base where base is usually in the range 51-60 for a 64-bit processor, 26-30 for 
# 32-bit and 12-13 for 16 bit
# The radix must be low enough so that the sum of L 2^(2*base) double length integers does not encroach on the sign bit, 
# that is bits(L)+2*base<2*WL, where L=n/base (rounded up) is the number of limbs in the number representation. But 
# at the same time radix should be high enough to minimize L. An underused top digit allows for lazy reduction.


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
    base=58

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

if base==0 :
    base=getbase(n)   # use default radix
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

overflow=False
bad_overflow=False
if (b-1)*(b-1)*mm*N >= 2**(2*WL) :
    overflow=True
    print("Possibility of overflow... using alternate method")

    if (N-1)*(b-1)**2 >= 2**(2*WL-3) :
        bad_overflow=True
if bad_overflow :
    print("Overflow requires extra resource")

# faster reduction
fred=False
if bits(N+1)+base+bits(mm)<WL :
    fred=True
    print("Tighter reduction")


# check if multiplication can be done single precision
EPM=False
if not overflow :
    if mm*(b-1)<2**WL : 
        EPM=True
        print("Fully Exploitable Pseudo-Mersenne detected")

# Do we need to propagate carries a bit further?
carry_on=False
if m*(2**(2*WL-base+xcess)+2**(base-xcess)) >= 2**(2*base) :
    print("Must propagate carries one limb further in second pass ")
    carry_on=True

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
        print("#include <time.h>\n")

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

print("Timing code is in file time.c")

# to determine code size
makestatic=False
with open('field.c', 'w') as f:
    with redirect_stdout(f):
        header()
        functions()
f.close()

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