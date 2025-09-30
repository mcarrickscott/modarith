# Python program to generate reasonably efficient C/C++ modular arithmetic code for pseudo-mersenne primes using ARM NEON SIMD extensions
# Modulus should be a pseudo-mersenne of the form 2^n-m
# uses unsaturated radix
#
# This version generates code that uses NEON intrinsics for 32-bit ARM processors
# Can execute two operations in parallel
#
# In particular this script generates code for primes like
#
# X25519 = 2^255-19
# PM266 = 2^266-3
# NIST521 = 2^521-1
#
# requires addchain utility in the path - see https://github.com/mmcloughlin/addchain 
#
# Execute this program as: python3 pseudo_neon.py X25519
# Production code is output to file field.c
# 
# Mike Scott 26th August 2025
# TII
#

# Some default settings

embedded=False  # If True then functions to start and stop a performance counter (cycles or micro-seconds) must be made available
                # If True no timing executable is created
cyclesorsecs=True     # if embedded, count cycles otherwise seconds
arduino=False   # set True if embedded and using arduino for timings
if arduino :    # If arduino, count microseconds
    cyclesorsecs=False 
compiler="gcc" # gcc, clang or icx (inlining can sometimes cause icx to hang)
cyclescounter=True # use Bernstein's cpucycle counter, otherwise just provide timings
decoration=False # decorate function names to avoid name clashes
formatted=True # pretty up the final output
inline=True # consider encouraging inlining
generic=True # set to False if algorithm is known in advance, in which case modadd and modsub can be faster - see https://eprint.iacr.org/2017/437. Set False for RFC7748 implementation.
check=False # run cppcheck on the output
scale=1 # set to 10 or 100 for faster timing loops. Default to 1
fully_propagate=False # recommended set to True for Production code

import sys
import subprocess

# Determine optimal radix
def getbase(n) :
    limbs=int(n/32)
    limit=64
    while (True) :
        limbs=limbs+1
        if limbs==1 :   # must be at least 2 limbs
            limbs=2
        base=int(16)
        while limbs*base<n :
            base=base+1
        if base>29 :  # OK, base too large (moddadd/sub needs 3 bits), another limb required
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
        pw.append(tp%b)        #str+="\tsspint carry=(sspint)n[0]>>{}u;\n".format(base)
        tp=tp>>base
        i=i+1
    return pw


def intrinsics() :
    str="// t=a*c where t, a is 64 bits and c is a small constant\n"
    str+="static inline dpint vmulct_n_u64(dpint a,int b) {\n"
    str+="\tudpint pp1=vmull_n_u32(vmovn_u64(a),b);\n"
    str+="\tudpint pp2=vmull_n_u32(vmovn_u64(vshrq_n_u64(a,32)),b);\n"
    str+="\treturn vaddq_u64(pp1,vshlq_n_u64(pp2,32));\n"  
    str+="}\n" 
    str+="// set each lane to a constant\n"
    str+="static inline spint tospint(int c0,int c1) {\n"
    str+="\tuint32_t s[2];\n"
    str+="\ts[0]=c0; s[1]=c1;\n"
    str+="\treturn vld1_u32(s);\n"
    str+="}\n" 
    return str

#conditional add of p
def caddp() :
    str=""
    str+="\tbot=vand_u32(bot,carry);\n"
    str+="\ttop=vand_u32(top,carry);\n"
    str+="\tn[0]=vsub_u32(n[0],bot);\n"
    str+="\tn[{}]=vadd_u32(n[{}],top);\n".format(N-1,N-1)

    #str+="\tn[0]-=((spint){}u)&carry;\n".format(m*x)
    #str+="\tn[{}]+=((spint)0x{:x}u)&carry;\n ".format(N-1,x*TW)
    return str

#add p
def addp() :
    str=""
    str+="\tn[0]=vsub_u32(n[0],bot);\n"
    str+="\tn[{}]=vadd_u32(n[{}],top);\n".format(N-1,N-1)    

    #str+="\tn[0]-=(spint){}u;\n".format(m*x)
    #str+="\tn[{}]+=(spint)0x{:x}u;\n".format(N-1,x*TW)
    return str

#subtract p
def subp() :
    str=""
    str+="\tn[0]=vadd_u32(n[0],bot);\n"
    str+="\tn[{}]=vsub_u32(n[{}],top);\n".format(N-1,N-1)    

    #str+="\tn[0]+=(spint){}u;\n".format(m*x)
    #str+="\tn[{}]-=(spint)0x{:x}u;\n".format(N-1,x*TW)
    return str

# Propagate carries.
def prop(n) :
    str="//propagate carries\n"
    str+="static spint inline prop(spint *n) {\n"
    str+="\tint i;\n"
    str+="\tspint mask=vdup_n_u32((1<<{}u)-1);\n".format(base)
    #str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

    str+="\tsspint carry=vreinterpret_s32_u32 (n[0]);\n"
    str+="\tcarry=vshr_n_s32 (carry, {}u);\n".format(base)

    #str+="\tsspint carry=(sspint)n[0];\n"
    #str+="\tcarry>>={}u;\n".format(base)

    str+="\tn[0]=vand_u32 (n[0], mask);\n"
    #str+="\tn[0]&=mask;\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N-1)

    str+="\t\tcarry=vadd_s32(carry,vreinterpret_s32_u32 (n[i]));\n"
    str+="\t\tn[i]=vand_u32(vreinterpret_u32_s32(carry),mask);\n"
    str+="\t\tcarry=vshr_n_s32 (carry, {}u);\n".format(base)

    #str+="\t\tcarry+=(sspint)n[i];\n"
    #str+="\t\tn[i] = (spint)carry & mask;\n"
    #str+="\t\tcarry>>={}u;\n".format(base)
    str+="\t}\n"
    str+="\tn[{}]=vadd_u32(n[{}],vreinterpret_u32_s32(carry));\n".format(N-1,N-1)

    #str+="\tn[{}]+=(spint)carry;\n".format(N-1)

    str+="\treturn vreinterpret_u32_s32(vshr_n_s32(vreinterpret_s32_u32(n[{}]),{}u));\n}}\n".format(N-1,31)
    #str+="\treturn -((n[{}]>>1)>>{}u);\n}}\n".format(N-1,32-2)
    return str


#propagate carries and add p if negative, propagate carries again
def flat(n) :
    mask=(1<<base)-1
    str="//propagate carries and add p if negative, propagate carries again\n"
    if makestatic :
        str+="static "
    if inline and makestatic :
        str+="spint inline flatten(spint *n) {\n"
    else :
        str+="spint flatten(spint *n) {\n"
    str+="\tspint bot=vdup_n_u32({}u);\n".format(m)
    str+="\tspint top=vdup_n_u32(0x{:x}u);\n".format(TW)            
    str+="\tspint carry=prop(n);\n"
    str+=caddp()
    str+="\t(void)prop(n);\n"
    str+="\tspint mask=vdup_n_u32(1);\n"
    str+="\treturn vand_u32(carry,mask);\n"

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
    str+="\tspint bot=vdup_n_u32({}u);\n".format(m)
    str+="\tspint top=vdup_n_u32(0x{:x}u);\n".format(TW)           
    str+=subp()
    str+="\treturn flatten(n);\n"
    str+="}\n"
    return str

#modular addition)"
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
    str+="\tspint bot=vdup_n_u32({}u);\n".format(m*2)
    str+="\tspint top=vdup_n_u32(0x{:x}u);\n".format(2*TW)        
    for i in range(0,N) :
        str+="\tn[{}]=vadd_u32(a[{}],b[{}]);\n".format(i,i,i)

        #str+="\tn[{}]=a[{}]+b[{}];\n".format(i,i,i)
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
        str+="\tspint bot=vdup_n_u32({}u);\n".format(m*2)
        str+="\tspint top=vdup_n_u32(0x{:x}u);\n".format(2*TW)       
    else :
        str+="\tspint bot=vdup_n_u32({}u);\n".format(m*mp)
        str+="\tspint top=vdup_n_u32(0x{:x}u);\n".format(mp*TW)              
    for i in range(0,N) :
        str+="\tn[{}]=vsub_u32(a[{}],b[{}]);\n".format(i,i,i)
        #str+="\tn[{}]=a[{}]-b[{}];\n".formaflattent(i,i,i)
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
    str+="\tspint zero=vdup_n_u32(0);\n"
    if not algorithm :
        str+="\tspint bot=vdup_n_u32({}u);\n".format(m*2)
        str+="\tspint top=vdup_n_u32(0x{:x}u);\n".format(2*TW)       
    else :
        str+="\tspint bot=vdup_n_u32({}u);\n".format(m*mp)
        str+="\tspint top=vdup_n_u32(0x{:x}u);\n".format(mp*TW)            
    for i in range(0,N) :
        str+="\tn[{}]=vsub_u32(zero,b[{}]);\n".format(i,i)
        #str+="\tn[{}]=(spint)0-b[{}];\n".format(i,i)
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
        if EPM :
            if first :
                str+="\t"
            else :
                str+=" "
            str+="t=vmlal_u32(t,ma{},b[{}]);".format(k,L)
            #str+="t+=(dpint)ma{}*(dpint)b[{}];".format(k,L)
        else :
            if first :
                str+="\ttt=vmull_u32(a[{}],b[{}]);".format(k,L)
                #str+="\ttt=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
                first=False
            else :
                str+=" tt=vmlal_u32(tt,a[{}],b[{}]);".format(k,L)
                #str+=" tt+=(dpint)a[{}]*(dpint)b[{}];".format(k,L)
        L-=1
        k+=1
    if row<N-1:
        if overflow :
            str+=" lo=vand_u32(vmovn_u64(tt),mask);"
            #str+=" lo=(spint)tt & mask;"
            if row==0 :
                str+=" t=vmlal_n_u32(t,lo,0x{:x});".format(mm)
                #str+=" t+=(dpint)lo*(dpint)0x{:x};".format(mm)
            else :
                if bad_overflow :
                    str+=" t=vadd_u64(t,vmulct_n_u6(vaddw_u32(hi,lo),0x{:x}));".format(mm)
                    #str+=" t+=(hi+(dpint)lo)*(dpint)0x{:x};".format(mm)
                else :
                    str+=" t=vmlal_n_u32(t,vadd_u32(lo,hi),0x{:x});".format(mm)
                    #str+=" t+=(dpint)(spint)(lo+hi)*(dpint)0x{:x};".format(mm)
            if bad_overflow :
                str+= "hi=vshrq_n_u64(tt,{}u);".format(base)
                #str+=" hi=tt>>{}u;".format(base)
            else :
                str+=" hi=vmovn_u64(vshrq_n_u64(tt,{}u));".format(base)
                #str+=" hi=(spint)(tt>>{}u);".format(base)
        else :
            if not EPM :
                str+=" tt=vmulct_n_u64(tt,0x{:x});".format(mm)
                str+=" t=vaddq_u64(t,tt);"
                #str+=" tt*=0x{:x};".format(mm)
                #str+=" t+=tt;"
    else :
        first=True
    
    k=0
    while k<=row :
        if first :
            str+="\tt=vmlal_u32(t,a[{}],b[{}]);".format(k,row-k)
            #str+="\tt+=(dpint)a[{}]*(dpint)b[{}];".format(k,row-k)
            first=False
        else :
            str+=" t=vmlal_u32(t,a[{}],b[{}]);".format(k,row-k)
            #str+=" t+=(dpint)a[{}]*(dpint)b[{}];".format(k,row-k)
        k+=1
    if row==N-1 and overflow :
        if bad_overflow :
            str+=" t=vaddq_u64(t,vmulct_n_u64(hi,0x{:x}));".format(mm)
        else:
            str+=" t=vmlal_n_u32(t,hi,0x{:x});".format(mm)
        #str+=" t+=(dpint)hi*(dpint)0x{:x};".format(mm)
    str+=" spint v{}=vand_u32(vmovn_u64(t),mask); t=vshrq_n_u64(t,{}u);\n".format(row,base)
    #str+=" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
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
                str+="t=vmlal_u32(t,ma{},ta{});".format(k,L)
                #str+="t+=(udpint)ma{}*(udpint)ta{};".format(k,L)
            else :
                if first :
                    str+="\tt=vmlal_u32(t,ma{},a[{}]);".format(k,L)
                    #str+="\tt+=(udpint)ma{}*(udpint)a[{}];".format(k,L)
                    first=False
                else :
                    str+=" t=vmlal_u32(t,ma{},a[{}]);".format(k,L)
                    #str+=" t+=(udpint)ma{}*(udpint)a[{}];".format(k,L)
        else :
            if first :
                str+="\ttt=vmull_u32(a[{}],a[{}]);".format(k,L)
                #str+="\ttt=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
                first=False
            else :
                str+=" tt=vmlal_u32(tt,a[{}],a[{}]);".format(k,L)
                #str+=" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)

        L-=1
        k+=1
    if dble :
        if not EPM :
            str+=" tt=vaddq_u64(tt,tt);\n"
            #str+=" tt*=2;"            
    if k==L :
        if EPM :
            if first :
                str+="\t"
            else :
                str+=" "  
            str+="t=vmlal_u32(t,ma{},a[{}]);".format(k,k)
            #str+="t+=(udpint)ma{}*(udpint)a[{}];".format(k,k)
        else :
            if first :
                str+="\ttt=vmull_u32(a[{}],a[{}]);".format(k,k)
                #str+="\ttt=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
                first=False
            else :
                str+=" tt=vmlal_u32(tt,a[{}],a[{}]);".format(k,k)
                #str+=" tt+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
    first=True
    if row<N-1:
        if overflow :
            str+=" lo=vand_u32(vmovn_u64(tt),mask);"
            #str+=" lo=(spint)tt & mask;"
        else :
            if not EPM :
                str+=" tt=vmulct_n_u64(tt,0x{:x});".format(mm)
                str+=" t=vaddq_u64(t,tt);"                
                #str+=" tt*=0x{:x};".format(mm)
                #str+=" t+=tt;"
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
            str+="t=vmlal_u32(t,a[{}],ta{});\n".format(k,L)
            #str+="t+=(udpint)a[{}]*(udpint)ta{};".format(k,L)
        else :
            if first :
                str+="t2=vmull_u32(a[{}],a[{}]);".format(k,L)
                #str+="t2=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
                first=False
            else :
                str+=" t2=vmlal_u32(t2,a[{}],a[{}]);".format(k,L)
                #str+=" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,L)
        k+=1
        L-=1

    if dble :
        if not EPM :
            str+=" t2=vaddq_u64(t2,t2);"
            #str+=" t2*=2;"
    if k==L :
        if EPM :
            str+=" t=vmlal_u32(t,a[{}],a[{}]);".format(k,k)
            #str+=" t+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
        else :
            if first :
                str+="t2=vmull_u32(a[{}],a[{}]);".format(k,k)
                #str+="t2=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
                first=False
            else :
                str+=" t2=vmlal_u32(t2,a[{}],a[{}]);".format(k,k)
                #str+=" t2+=(udpint)a[{}]*(udpint)a[{}];".format(k,k)
 

    if overflow :
        if row==N-1 : 
            if bad_overflow :
                str+= "t=vaddq_u64(t,vmulct_n_u64(hi,0x{:x}));".format(mm)
                #str+=" t+=(udpint)hi*(udpint)0x{:x};".format(mm)
            else :
                str+=" t=vmlal_n_u32(t,hi,0x{:x});".format(mm)
                #str+=" t+=(udpint)hi*(udpint)0x{:x};".format(mm)
        else :
            if row==0 :
                str+=" t2=vmlal_n_u32(t2,lo,0x{:x});".format(mm)
                #str+=" t2+=(udpint)lo*(udpint)0x{:x};".format(mm)
            else :
                if bad_overflow :
                    str+=" t2=vadd_u64(t2,vmulct_n_u6(vaddw_u32(hi,lo),0x{:x}));".format(mm)
                    #str+=" t2+=(hi+(udpint)lo)*(udpint)0x{:x};".format(mm)
                else :
                    str+=" t2=vmlal_n_u32(t2,vadd_u32(lo,hi),0x{:x});".format(mm)
                    #str+=" t2+=(udpint)(spint)(lo+hi)*(udpint)0x{:x};".format(mm)
            if bad_overflow :
                str+=" hi=vshrq_n_u64(tt,{}u);".format(base)
                #str+=" hi=tt>>{}u;".format(base)
            else :
                str+=" hi=vmovn_u64(vshrq_n_u64(tt,{}u));".format(base)
                #str+=" hi=(spint)(tt>>{}u);".format(base)
    if not EPM :
        str+=" t=vaddq_u64(t,t2);"
        #str+=" t+=t2;"
    str+=" spint v{}=vand_u32(vmovn_u64(t),mask); t=vshrq_n_u64(t,{}u);\n".format(row,base)
    #str+=" spint v{}=(spint)t & mask; t=t>>{}u;\n".format(row,base)
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
        str+="\tspint ut=vmovn_u64(t);\n"
        #str+="\tspint ut=(spint)t;\n"  
        if xcess>0 :
            str+="\tspint smask=vdup_n_u32({});\n".format((1<<(base-xcess))-1)
            str+= "\tut=vadd_u32(vshl_n_u32(ut,{}),vshr_n_u32(v{},{}); v{}=vand_u32(v{},smask);\n".format(xcess,N-1,base-xcess,N-1,N-1)
            #str+= "\tut=(ut<<{})+(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)

        if m>1 :
            str+= "\tut=vmul_n_u32(ut,{});\n".format(m)
            #str+= "\tut*=0x{:x};\n".format(m)
            
        str+="\ts=vadd_u32(v0,vand_u32(ut,mask));\n"
        str+="\tc[0]=vand_u32(s,mask);\n"
        #str+= "\ts=v0+(ut & mask);\n"
        #str+= "\tc[0]=(s&mask);\n"

        if carry_on :
            str+="\tut=vadd_u32(vshr_n_u32(s,{}),vshr_n_u32(ut,{}));\n".format(base,base)
            str+="\ts=vadd_u32(v1,vand_u32(ut,mask));\n"
            str+="\tc[1]=vand_u32(s,mask);\n"
            
            #str+="\tut=(s>>{})+(ut>>{});\n".format(base,base)
            #str+="\ts=v1+(ut & mask);\n"
            #str+= "\tc[1]=(s&mask);\n"
            k+=1

        str+= "\tcarry=vadd_u32(vshr_n_u32(s,{}),vshr_n_u32(ut,{}));\n".format(base,base)
        #str+= "\tcarry=(s>>{})+(ut>>{});\n".format(base,base)
            
    else :
        str+="\tudpint ut=(udpint)t;\n"    
        if xcess>0 :
            str+="\tspint smask=vdup_n_u32({});\n".format((1<<(base-xcess))-1)
            str+="\tut=vaddw_u32(vshlq_n_u64(ut,{}),vshr_n_u32(v{},{})); v{}=vand_u32(v{},smask);\n".format(xcess,N-1,base-xcess,N-1,N-1) 
            #str+= "\tut=(ut<<{})+(udpint)(v{}>>{}u); v{}&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)

        if m>1 :
            str+="\tut=vmulct_n_u64(ut,{});\n".format(m)
            #str+= "\tut*=0x{:x};\n".format(m)
        str+="\ts=vadd_u32(v0,vand_u32(vmovn_u64(ut),mask));\n"   
        str+="\tc[0]=vand_u32(s,mask);\n"         
        #str+= "\ts=v0+((spint)ut & mask);\n"
        #str+= "\tc[0]=(s&mask);\n"

        if carry_on :
            str+="\tut=vaddw_u32(vshrq_n_u64(ut,{}),vshr_n_u32(s,{}));\n".format(base,base)
            str+="\ts=vadd_u32(v1,vand_u32(vmovn_u64(ut),mask));\n"
            str+="\tc[1]=vand_u32(s,mask);\n"
            #str+="\tut=(udpint)(s>>{})+(ut>>{});\n".format(base,base)
            #str+="\ts=v1+((spint)ut & mask);\n"
            #str+= "\tc[1]=(s&mask);\n"
            k+=1

        str+= "\tcarry=vadd_u32(vshr_n_u32(s,{}),vmovn_u64(vshrq_n_u64(ut,{})));\n".format(base,base)
        #str+= "\tcarry=(s>>{})+(spint)(ut>>{});\n".format(base,base)

    k=k+1
    str+="\tc[{}]=vadd_u32(v{},carry);\n".format(k,k)
    #str+= "\tc[{}]=v{}+carry;\n".format(k,k)

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

    str+="\tdpint t=vdupq_n_u64(0);\n"
    #str+="\tdpint t=0;\n"

    if  EPM  :
        for i in range(1,N) :
            str+="\tspint ma{}=vmul_n_u32(a[{}],0x{:x});\n".format(i,i,mm);
            #str+="\tspint ma{}=a[{}]*(spint)0x{:x};\n".format(i,i,mm)
    else :
        str+="\tdpint tt;\n"

    if overflow :
        str+="\tspint lo;\n"
        if bad_overflow :
            str+="\tdpint hi;\n"   # could overflow single type
        else :
            str+="\tspint hi;\n"
    #str+="\tuspint carry,s,mask=((uspint)1<<{}u)-(uspint)1;\n".format(base)
    
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=vdup_n_u32((1<<{}u)-1);\n".format(base)    
    #str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

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
    
    str+="\tudpint t=vdupq_n_u64(0);\n"
    #str+="\tudpint t=0;\n"

    if  EPM  :
        for i in range(1,N) :
            str+="\tspint ta{}=vmul_n_u32(a[{}],2);\n".format(i,i)
            #str+="\tspint ta{}=a[{}]*(spint)2;\n".format(i,i)
        for i in range(1,N) :
            str+="\tspint ma{}=vmul_n_u32(a[{}],0x{:x});\n".format(i,i,mm)
            #str+="\tspint ma{}=a[{}]*(spint)0x{:x};\n".format(i,i,mm)
    else :
        str+="\tudpint tt;\n"
        str+="\tudpint t2;\n"
    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint mask=vdup_n_u32((1<<{}u)-1);\n".format(base) 
    #str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

    if overflow :
        str+="\tspint lo;\n"
        if bad_overflow :
            str+="\tudpint hi;\n"
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
        
    str+="\tudpint t=vdupq_n_u64(0);\n"
    #str+="\tudpint t=0;\n"

    str+="\tspint carry;\n"
    str+="\tspint s;\n"
    str+="\tspint bw=vdup_n_u32(b);\n"
    str+="\tspint mask=vdup_n_u32((1<<{}u)-1);\n".format(base)
    #str+="\tspint mask=((spint)1<<{}u)-(spint)1;\n".format(base)

    for i in range(0,N) :
        str+="\tt=vmlal_u32(t,a[{}],bw); ".format(i)
        str+="spint v{}=vand_u32(vmovn_u64(t),mask); t=vshrq_n_u64(t,{}u);\n".format(i,base)
        
        #str+="\tt+=(udpint)a[{}]*(udpint)b; ".format(i)
        #str+="spint v{}=(spint)t & mask; t=t>>{}u;\n".format(i,base)

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
        str+="\t\tspint one=vdup_n_u32(1);\n";
        str+="\t\tspint d=vsub_u32(one,modis1{}(b));\n".format(DECOR)
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
    str+="\tspint one=vdup_n_u32(1);\n"
    str+="\tspint d=vdup_n_u32(0);\n"
    #str+="\tspint d=0;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\td=vorr_u32(d,c[i]);\n\t}\n"
    #str+="\t\td|=c[i];\n\t}\n"
    #str+="\tc0=(spint)c[0];\n"
    str+="\treturn vand_u32(vand_u32(one,vshr_n_u32(vsub_u32(d,one),{}u)) ,    vshr_n_u32(vsub_u32(veor_u32(c[0],one)   ,one),{}u)    );\n}}\n".format(base,base)   
    #str+="\treturn (   (spint)1 & ((d-(spint)1)>>{}u)      &    (         (   (c0^(spint)1)  -  (spint)1    )    >>{}u)      );\n}}\n".format(base,base)
    return str

#test for zero
def modis0(n) :
    str="//is zero?\n"
    if makestatic :
        str+="static "
    str+="spint modis0{}(const spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tspint one=vdup_n_u32(1);\n"    
    str+="\tspint d=vdup_n_u32(0);\n"    
    #str+="\tspint d=0;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\td=vorr_u32(d,c[i]);\n\t}\n"    
    #str+="\t\td|=c[i];\n\t}\n"
    str+="\treturn vand_u32(  vshr_n_u32(vsub_u32(d,one),{}u)    ,one);\n}}\n".format(base)
    #str+="\treturn ((spint)1 & ((d-(spint)1)>>{}u));\n}}\n".format(base)
    return str

#set to zero
def modzer() :
    str="//set to zero\n"
    if makestatic :
        str+="static "
    str+="void modzer{}(spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=vdup_n_u32(0);\n"
    #str+="\t\ta[i]=0;\n"
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
    str+="\ta[0]=vdup_n_u32(1);\n"
    #str+="\ta[0]=1;\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=vdup_n_u32(0);\n"
    #str+="\t\ta[i]=0;\n"
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
    str+="\ta[0]=vdup_n_u32(x);\n"    
    #str+="\ta[0]=(spint)x;\n"
    str+="\tfor (i=1;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=vdup_n_u32(0);\n"    
    #str+="\t\ta[i]=0;\n"
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
    str+="void __attribute__ ((noinline)) modcsw{}(spint b,volatile spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c0,c1,s,t,w,v,aux;\n"
    
    str+="\tstatic uint32_t R[2]={0,0};\n"
    str+="\tspint one=vdup_n_u32(1);\n"
    str+="\tR[0]+=0x5aa5a55au;\n"
    str+="\tR[1]+=0x7447e88eu;\n"    
    
    #str+="\tstatic spint R=0;\n"
    #str+="\tR+=0x5aa5a55au;\n"
    str+="\tw=vld1_u32(R);\n"
    #str+="\tw=R;\n"
    
    str+="\tc0=vand_u32(vmvn_u32(b),vadd_u32(w,one));\n"
    str+="\tc1=vadd_u32(b,w);\n"    
    #str+="\t\tc0=(~b)&(w+1);\n"
    #str+="\t\tc1=b+w;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ts=g[i]; t=f[i];\n"
    str+="\t\tv=vmul_u32(w,vadd_u32(t,s));\n"
    #str+="\t\tv=w*(t+s);\n"
    str+="\t\tf[i]=aux=vadd_u32(vmul_u32(c0,t),vmul_u32(c1,s));\n"    
    str+="\t\tf[i]=vsub_u32(aux,v);\n"
    #str+="\t\tf[i] = aux = c0*t+c1*s;\n"
    #str+="\t\tf[i] = aux - v;\n"
    str+="\t\tg[i]=aux=vadd_u32(vmul_u32(c0,s),vmul_u32(c1,t));\n"    
    str+="\t\tg[i]=vsub_u32(aux,v);\n\t}\n"
    #str+="\t\tg[i] = aux = c0*s+c1*t;\n"
    #str+="\t\tg[i] = aux - v;\n\t}\n"
    str+="}\n"
    return str

#conditional move
def modcmv() :
    str="//conditional move g to f if d=1\n"
    str+="//strongly recommend inlining be disabled using compiler specific syntax\n"
    if makestatic :
        str+="static "
    str+="void __attribute__ ((noinline)) modcmv{}(spint b,const spint *g,volatile spint *f) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint c0,c1,s,t,w,aux;\n"
    str+="\tstatic uint32_t R[2]={0,0};\n"
    str+="\tspint one=vdup_n_u32(1);\n"
    str+="\tR[0]+=0x5aa5a55au;\n"
    str+="\tR[1]+=0x7447e88eu;\n"
    str+="\tw=vld1_u32(R);\n"
    #str+="\tw=R;\n"
    str+="\tc0=vand_u32(vmvn_u32(b),vadd_u32(w,one));\n"
    str+="\tc1=vadd_u32(b,w);\n"
    #str+="\t\tc0=(~b)&(w+1);\n"
    #str+="\t\tc1=b+w;\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ts=g[i]; t=f[i];\n"
    str+="\t\tf[i]=aux=vadd_u32(vmul_u32(c0,t),vmul_u32(c1,s));\n"
    str+="\t\tf[i]=vsub_u32(aux,vmul_u32(w,vadd_u32(t,s)));\n\t}\n"
    #str+="\t\tf[i] = aux = c0*t+c1*s;\n"
    #str+="\t\tf[i] = aux - w*(t+s);\n\t}\n"
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
    str+="\tspint mask=vdup_n_u32(0x{:x});\n".format(mask)
    str+="\ta[{}]=vorr_u32(vshl_n_u32(a[{}],n),vshr_n_u32(a[{}],{}u-n));\n".format(N-1,N-1,N-2,base)
    #str+="\ta[{}]=((a[{}]<<n)) | (a[{}]>>({}u-n));\n".format(N-1,N-1,N-2,base)
    str+="\tfor (i={};i>0;i--) {{\n".format(N-2)
    str+="\t\ta[i]=vorr_u32(vand_u32(vshl_n_u32(a[i],n),mask)    ,vshr_n_u32(a[i-1],{}u-n)    );\n\t}}\n".format(base)
    #str+="\t\ta[i]=((a[i]<<n)&(spint)0x{:x}) | (a[i-1]>>({}u-n));\n\t}}\n".format(mask,base)
    str+="\ta[0]=vand_u32(vshl_n_u32(a[0],n),mask);\n"
    #str+="\ta[0]=(a[0]<<n)&(spint)0x{:x};\n".format(mask)
    str+="}\n"
    return str 

#shift right
def modshr(n) :
    N=getN(n)
    mask=(1<<base)-1
    str="//shift right by less than a word. Return shifted out part\n"
    if makestatic :
        str+="static "
    str+="spint modshr{}(unsigned int n,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint mask=vdup_n_u32(0x{:x});\n".format(mask)
    str+="\tspint mskn=vdup_n_u32((1<<n)-1);\n"    
    str+="\tspint r=vand_u32(a[0],mskn);\n"
    #str+="\tspint r=a[0]&(((spint)1<<n)-(spint)1);\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(N-1)
    str+="\t\ta[i]=vorr_u32(vshr_n_u32(a[i],n),vand_u32(vshl_n_u32(a[i+1],{}u-n),mask));\n\t}}\n".format(base)
    #str+="\t\ta[i]=(a[i]>>n) | ((a[i+1]<<({}u-n))&(spint)0x{:x});\n\t}}\n".format(base,mask)
    str+="\ta[{}]=vshr_n_u32(a[{}],n);\n".format(N-1,N-1)
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
    str+="\tspint one=vdup_n_u32(1);\n"
    str+="\tmodzer{}(a);\n".format(DECOR)
    str+="\tif (r>={}*8) return;\n".format(Nbytes)
    str+="\ta[n]=one; a[n]=vshl_n_u32(a[n],m);\n}\n"
    #str+="\ta[n]=1; a[n]<<=m;\n}\n"
    return str

#export to byte array
def modexp() :
    str="//export to byte arrays\n"
    if makestatic :
        str+="static "
    str+="void modexp{}(const spint *a,char *b,char *e) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tuint32_t s[2];\n"
    str+="\tspint mask=vdup_n_u32(0xff);\n"
    str+="\tspint c[{}];\n".format(N)
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tfor (i={};i>=0;i--) {{\n".format(Nbytes-1)
    
    str+="\t\tvst1_u32(s,vand_u32(c[0],mask));\n"
    str+="\t\tb[i]=(char)s[0];\n"
    str+="\t\tif (e!=NULL) e[i]=(char)s[1];\n"
    #str+="\t\tb[i]=c[0]&(spint)0xff;\n"
    
    str+="\t\t(void)modshr{}(8,c);\n\t}}\n".format(DECOR)
    str+="}\n"
    return str 

#import from byte array
def modimp() :
    str="//import from byte arrays\n"
    str+="//returns 1 if in range, else 0\n"
    if makestatic :
        str+="static "
    str+="spint modimp{}(const char *b, const char *e,spint *a) {{\n".format(DECOR)
    str+="\tint i;\n"
    str+="\tspint res;\n"
    str+="\tuint32_t s[2];\n"
        
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    str+="\t\ta[i]=vdup_n_u32(0);\n\t}\n"
    #str+="\t\ta[i]=0;\n\t}\n"
    str+="\tfor (i=0;i<{};i++) {{\n".format(Nbytes)
    str+="\t\tmodshl{}(8,a);\n".format(DECOR)
    
    str+="\t\ts[0]=(uint32_t)b[i]; \n"
    str+="\t\tif (e!=NULL) s[1]=(uint32_t)e[i];\n"
    str+="\t\telse s[1]=0;\n"
    str+="\t\ta[0]=vadd_u32(a[0],vld1_u32(s));\n\t}\n"
    #str+="\t\ta[0]+=(spint)(unsigned char)b[i];\n\t}\n"
    
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
    str+="\tspint one=vdup_n_u32(1);\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\treturn vand_u32(c[0],one);\n"
    #str+="\treturn c[0]%2;\n"
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
    str+="\tspint one=vdup_n_u32(1);\n"
    str+="\tspint eq=one;\n"
    str+="\tredc{}(a,c);\n".format(DECOR)
    str+="\tredc{}(b,d);\n".format(DECOR)
    str+="\tfor (i=0;i<{};i++) {{\n".format(N)
    
    str+="\t\teq=vand_u32(eq, vand_u32(vshr_n_u32(vsub_u32(veor_u32(c[i],d[i]),one),{}) ,one));\n\t}}\n".format(base)
    #str+="\t\teq&=(((c[i]^d[i])-1)>>{})&1;\n\t}}\n".format(base)
    
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
    str+="\tuint32_t s[2];\n"   
    if not embedded :
        str+="\tuint64_t start,finish;\n"
        str+="\tclock_t begin;\n"
        str+="\tint elapsed;\n"
    else :
        str+="\tlong start,finish;\n"
    str+="\t"
    for i in range(0,N) :
        str+="s[0]={}; ".format(hex(rap[i]))
        str+="s[1]={}; ".format(hex(rbp[i]))
        str+="x[{}]=vld1_u32(s); ".format(i)
        #str+="x[{}]={}; ".format(i,hex(rap[i]))
    str+="\n\t"
    for i in range(0,N) :
        str+="s[0]={}; ".format(hex(rbp[i]))
        str+="s[1]={}; ".format(hex(rap[i]))
        str+="y[{}]=vld1_u32(s); ".format(i)        
        #str+="y[{}]={}; ".format(i,hex(rbp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    str+="\tnres{}(y,y);\n".format(DECOR)

    if embedded :
        if arduino :
            str+="\tstart=micros();\n"
        else :
            str+="\t//provide code to start counter, start=?;\n"
        str+="\tfor (j=0;j<200;j++) {\n"
        str+="\t\tmodmul{}(x,y,z);\n".format(DECOR)
        str+="\t\tmodmul{}(z,x,y);\n".format(DECOR)
        str+="\t\tmodmul{}(y,z,x);\n".format(DECOR)
        str+="\t\tmodmul{}(x,y,z);\n".format(DECOR)
        str+="\t\tmodmul{}(z,x,y);\n".format(DECOR)
        str+="\t}\n"
        if arduino :
            str+="\tfinish=micros();\n"
        else :
            str+="\t//provide code to stop counter, finish=?;\n"
        str+="\tredc{}(z,z);\n".format(DECOR)
        str+="\tvst1_u32(s,z[0]);\n"
        if arduino :
            str+='\tSerial.print("modmul usecs= "); Serial.println((finish-start)/1000);\n'
        else :
            if cyclesorsecs :
                str+='\tprintf("modmul check %x Clock cycles= %d\\n",(int)s[0]&0xFFFFFF,(int)((finish-start)/1000));\n'
            else :
                str+='\tprintf("modmul check %x Nanosecs= %d\\n",(int)s[0]&0xFFFFFF,(int)((finish-start)));\n'
        str+="}\n"
    else :
        if cyclescounter :
            str+="\tstart=cpucycles();\n"
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
        str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10*scale)
        str+="\tredc{}(z,z);\n".format(DECOR)
        str+="\tvst1_u32(s,z[0]);\n"
        if cyclescounter :
            str+='\tprintf("modmul check 0x%06x Clock cycles= %d Nanosecs= %d\\n",(int)s[0]&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
        else :
            str+='\tprintf("modmul check 0x%06x Nanosecs= %d\\n",(int)s[0]&0xFFFFFF,elapsed);\n'
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
    str+="\tuint32_t s[2];\n"   
    if not embedded :
        str+="\tuint64_t start,finish;\n"
        str+="\tclock_t begin;\n"
        str+="\tint elapsed;\n"
    else :
        str+="\tlong start,finish;\n"

    str+="\t"
    for i in range(0,N) :
        str+="s[0]={}; ".format(hex(rp[i]))
        str+="s[1]={}; ".format(hex(rp[i]))
        str+="x[{}]=vld1_u32(s); ".format(i)        
        #str+="x[{}]={}; ".format(i,hex(rp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    if embedded :
        if arduino :
            str+="\tstart=micros();\n"
        else :
            str+="\t//provide code to start counter, start=?;\n"
        str+="\tfor (j=0;j<500;j++) {\n"
        str+="\t\tmodsqr{}(x,z);\n".format(DECOR)
        str+="\t\tmodsqr{}(z,x);\n".format(DECOR)
        str+="\t}\n"
        if arduino :
            str+="\tfinish=micros();\n"
        else :
            str+="\t//provide code to stop counter, finish=?;\n"
        str+="\tredc{}(z,z);\n".format(DECOR)
        str+="\tvst1_u32(s,z[0]);\n"        
        if arduino :
            str+='\tSerial.print("modsqr usecs= "); Serial.println((finish-start)/1000);\n'
        else :
            if cyclesorsecs :
                str+='\tprintf("modsqr check %x Clock cycles= %d\\n",(int)s[0]&0xFFFFFF,(int)((finish-start)/1000));\n'
            else :
                str+='\tprintf("modsqr check %x Nanosecs= %d\\n",(int)s[0]&0xFFFFFF,(int)((finish-start)));\n'
        str+="}\n"
    else :
        if cyclescounter :
            str+="\tstart=cpucycles();\n"
        str+="\tbegin=clock();\n"
        str+="\tfor (i=0;i<{};i++)\n".format(100000//scale)
        str+="\t\tfor (j=0;j<500;j++) {\n"
        str+="\t\t\tmodsqr{}(x,z);\n".format(DECOR)
        str+="\t\t\tmodsqr{}(z,x);\n".format(DECOR)
        str+="\t\t}\n"
        if cyclescounter :
            str+="\tfinish=cpucycles();\n"
        str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10*scale)
        str+="\tredc{}(z,z);\n".format(DECOR)
        str+="\tvst1_u32(s,z[0]);\n"         
        if cyclescounter :
            str+='\tprintf("modsqr check 0x%06x Clock cycles= %d Nanosecs= %d\\n",(int)s[0]&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000000//scale)
        else :
            str+='\tprintf("modsqr check 0x%06x Nanosecs= %d\\n",(int)s[0]&0xFFFFFF,elapsed);\n'
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
    str+="\tuint32_t s[2];\n"      
    if not embedded :
        str+="\tuint64_t start,finish;\n"
        str+="\tclock_t begin;\n"
        str+="\tint elapsed;\n"
    else :
        str+="\tlong start,finish;\n"

    str+="\t"
    for i in range(0,N) :
        str+="s[0]={}; ".format(hex(rp[i]))
        str+="s[1]={}; ".format(hex(rp[i]))
        str+="x[{}]=vld1_u32(s); ".format(i)             
        #str+="x[{}]={}; ".format(i,hex(rp[i]))
    str+="\n"

    str+="\tnres{}(x,x);\n".format(DECOR)
    if embedded :
        if arduino :
            str+="\tstart=micros();\n"
        else :
            str+="\t//provide code to start counter, start=?;\n"
        str+="\tmodinv{}(x,NULL,z);\n".format(DECOR)
        if arduino :
            str+="\tfinish=micros();\n"
        else :
            str+="\t//provide code to stop counter, finish=?;\n"
        str+="\tredc{}(z,z);\n".format(DECOR)
        str+="\tvst1_u32(s,z[0]);\n"          
        if arduino :
            str+='\tSerial.print("modinv usecs= "); Serial.println(finish-start);\n'
        else :
            if cyclesorsecs :
                str+='\tprintf("modinv check %x Clock cycles= %d\\n",(int)s[0]&0xFFFFFF,(int)(finish-start));\n'
            else :
                str+='\tprintf("modinv check %x Microsecs= %d\\n",(int)s[0]&0xFFFFFF,(int)(finish-start));\n'
        str+="}\n"
    else :
        if cyclescounter :
            str+="\tstart=cpucycles();\n"
        str+="\tbegin=clock();\n"
        str+="\tfor (i=0;i<{};i++) {{\n".format(50000//scale)
        str+="\t\t\tmodinv{}(x,NULL,z);\n".format(DECOR)
        str+="\t\t\tmodinv{}(z,NULL,x);\n".format(DECOR)
        str+="\t}\n"
        if cyclescounter :
            str+="\tfinish=cpucycles();\n"
        str+="\telapsed = {}*(clock() - begin) / CLOCKS_PER_SEC;\n".format(10000*scale)
        str+="\tredc{}(z,z);\n".format(DECOR)
        str+="\tvst1_u32(s,z[0]);\n"          
        if cyclescounter :
            str+='\tprintf("modinv check 0x%06x Clock cycles= %d Nanosecs= %d\\n",(int)s[0]&0xFFFFFF,(int)((finish-start)/{}ULL),elapsed);\n'.format(100000//scale)
        else :
            str+='\tprintf("modinv check 0x%06x Microsecs= %d\\n",(int)s[0]&0xFFFFFF,elapsed);\n'
        str+="}\n"
    return str

def header() :
    print("\n//Automatically generated modular arithmetic C code for pseudo-Mersenne primes")
    print("//Command line : python {} {}".format(sys.argv[0], sys.argv[1]))
    print("//Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)\n")
    print("#include <stdio.h>")
    print("#include <stdint.h>\n")
    print("#include <arm_neon.h>")
    print("#define sspint int32x2_t")
    print("#define spint uint32x2_t")

    print("#define udpint uint64x2_t")

    print("#define dpint uint64x2_t\n")
    print("#define Wordlength{} 32".format(DECOR))
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
    if cyclescounter :
        str+='\tprintf("%s %s\\n",cpucycles_implementation(),cpucycles_version());\n'
    str+="\ttime_modmul();\n"
    str+="\ttime_modsqr();\n"
    str+="\ttime_modinv();\n"
    str+="\treturn 0;\n"
    str+="}\n"
    return str

if len(sys.argv)!=2 :
    print("Syntax error")
    print("Valid syntax - python pseudo_neon.py <prime> OR <prime name>")
    print("For example - python pseudo_neon.py X25519")
    print("For example - python pseudo_neon.py 2**255-19")
    exit(2)

prime=sys.argv[1]

n=0
base=0
m=0

# Note that number base MUST be less than wordlength 32 (unsaturated).
# The radix for an n-bit prime is 2^base where base is usually in the range 51-60 for a 64-bit processor, 26-30 for 
# 32-bit and 12-13 for 16 bitmakestatic
# The radix must be low enough so that the sum of L 2^(2*base) double length integers does not encroach on the sign bit, 
# that is bits(L)+2*base<2*32, where L=n/base (rounded up) is the number of limbs in the number representation. But 
# at the same time radix should be high enough to minimize L. An underused top digit allows for lazy reduction.


### Start of user editable area

inline=False

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
    base=29


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
if mm >= 2**(31) : #or m*((N+2)*2**(64)) >= 2**(2*base+31-xcess) :
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
print("Compiler is "+compiler)
 
print("Using standard Comba for modmul")

overflow=False
bad_overflow=False
if (b-1)*(b-1)*mm*N >= 2**(64) :
    overflow=True
    print("Possibility of overflow... using alternate method")

    if (N-1)*(b-1)**2 >= 2**(2*32-3) :
        bad_overflow=True
if bad_overflow :
    print("Overflow requires extra resource")

# faster reduction
fred=False
if bits(N+1)+base+bits(mm)<32 :
    fred=True
    print("Tighter reduction")


# check if multiplication can be done single precision
EPM=False
if not overflow :
    if mm*(b-1)<2**32 : 
        EPM=True
        print("Fully Exploitable Pseudo-Mersenne detected")

# Do we need to propagate carries a bit further?
carry_on=False
if m*(2**(64-base+xcess)+2**(base-xcess)) >= 2**(2*base) :
    print("Must propagate carries one limb further in second pass ")
    carry_on=True

# generate 32 bit modular multiplication and squaring code for Mersennmakestatice prime 2^n-m
# where the unsaturated radix is 2^base and base<32 

from contextlib import redirect_stdout

DECOR=""
modulus=p
import random

makestatic=True
random.seed(42)
ra=random.randint(0,modulus-1)
rb=random.randint(0,modulus-1)
rs=random.randint(0,modulus-1)
ri=random.randint(0,modulus-1)
#if embedded :
subprocess.call("rm time.c", shell=True)

with open('time.c', 'w') as f:
    with redirect_stdout(f):
        header()
        if not embedded :
            if cyclescounter :
                print("#include <cpucycles.h>\n")
            print("#include <time.h>\n")
#        if use_rdtsc :
#            print("#include <x86intrin.h>\n")
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
        if not embedded :
            print(main())

f.close()

#maybe -march=rv64gc for RISC-V
#maybe -fPIC
if not embedded :   # Create timing program for this processor
    if cyclescounter :
        subprocess.call(compiler + " -march=native -mtune=native -O3 time.c -lcpucycles -o time", shell=True)
    else :
        subprocess.call(compiler + " -march=native -mtune=native -O3 time.c -o time", shell=True)
    print("For timings run ./time")
    #subprocess.call("rm time.c", shell=True)
else :
    print("Timing code is in file time.c")

# to determine code size
makestatic=False
with open('field.c', 'w') as f:
    with redirect_stdout(f):
        header()
        functions()
f.close()

subprocess.call(compiler+" -O3 -c -fmax-errors=100 field.c",shell=True)
subprocess.call("size field.o > size.txt",shell=True)

f=open('size.txt')
lines=f.readlines()
info=lines[1].split()

print("Code size using -O3 = ",info[0])

subprocess.call(compiler+" -Os -c field.c",shell=True)
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
