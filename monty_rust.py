# Python program to generate reasonably efficient C/C++ modular arithmetic code for any prime, on a 16, 32 or 64-bit processor
# IMPORTANT:- Uses Montgomery representation
# uses unsaturated radix
#
# In particular this script generates code for the NIST primes
#
# NIST256 = 2^256-2^224+2^192+2^96-1
# NIST384 = 2^384-2^128-2^96+2^32-1
#
# and the "Goldilocks" prime
#
# X448=2^448-2^224-1
#
# requires addchain utility in the path - see https://github.com/mmcloughlin/addchain 
#
# How to use. 
# (1) First execute this program: python3 monty_rust.py 64 NIST256. Output code is written to file field.rs or group.rs
# (2) All constants and inputs must be converted to Montgomery nresidue form by calling nres()
# (3) All final outputs must be converted back to integer form by calling redc()
#
# The modulus can either be entered directly, or derived from the hard-wired name of an elliptic curve
# Normally this will be the field prime modulus
# However if the modulus is entered in decimal with two leading 0s, it is instead assumed to be a group 
# order (a large prime factor of the number of points on a curve)
#
# Note that even though a modulus is represented to the base 2^29 rather than 2^32, it still retains some shape
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

cyclesorsecs=True # use a cpu cycle counter, otherwise just provide timings
karatsuba=False # default setting
formatted=True # pretty up the final output
inline=True # consider encouraging inlining
generic=True # set to False if algorithm is known in advance, in which case modadd and modsub can be faster - see https://eprint.iacr.org/2017/437
scale=1 # set to 10 or 100 for faster timing loops. Default to 1
makepublic=False # Make functions public
fully_propagate=False # recommended set to True for Production code

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

# Determine default optimal radix
def getbase(n) :
    limbs=int(n/WL)
    limit=2*WL
    if karatsuba :
        limit=2*WL-1
    while (True) :
        limbs=limbs+1
        if limbs==1 :   # must be at least 2 limbs
            limbs=2
        base=int(WL/2)
        while limbs*base<n or base-(n%base)<2 :  # insist on excess of 2 or more
            base=base+1
        if n%base==0 :
            _,E=process_prime(p,base,limbs)      # if excess=0, and if no extra limb required, bump base
            if not E :
                base=base+1
        while limbs*base<n or base-(n%base)<2 :  # insist on excess of 2 or more
            base=base+1
        if base>WL-3 :    # OK, base too large, another limb required
            continue
        if limbs*(2**base-1)**2 < 2**limit :   # worst case scenario for non-shaped primes
            break
    return base

# get number of limbs
def getN(n,base) :
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
def getR(n,base) :
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

#conditional add of p
def caddp(x) :
    str=""
    for i in range(0,N) :
        if ppw[i]==0 :
            continue
        if ppw[i]==-1:
            str+="\tn[{}]-={}&(carry as SPINT);\n".format(i,x)
            continue
        if ppw[i]<0:
            str+="\tn[{}]-=(0x{:x} as SPINT)&(carry as SPINT);\n".format(i,-ppw[i]*x)
            continue
        str+="\tn[{}]+=(0x{:x} as SPINT)&(carry as SPINT);\n".format(i,ppw[i]*x)
    if E:
        str+="\tn[{}]+=({}*q)&(carry as SPINT);\n".format(N-1,x)
    return str

#add p
def addp(x) :
    str=""
    for i in range(0,N) :
        if ppw[i]==0 :
            continue
        if ppw[i]==-1:
            str+="\tn[{}]-={};\n".format(i,x)
            continue
        if ppw[i]<0:
            str+="\tn[{}]-=0x{:x} as SPINT;\n".format(i,-ppw[i]*x)
            continue
        str+="\tn[{}]+=0x{:x} as SPINT;\n".format(i,ppw[i]*x)
    if E:
        str+="\tn[{}]+={}*q;\n".format(N-1,x)
    return str


#subtract p
def subp(x) :
    str=""
    for i in range(0,N) :
        if ppw[i]==0 :
            continue
        if ppw[i]==-1:
            str+="\tn[{}]+={};\n".format(i,x)
            continue
        if ppw[i]<0:
            str+="\tn[{}]+=0x{:x} as SPINT;\n".format(i,-ppw[i]*x)
            continue
        str+="\tn[{}]-=0x{:x} as SPINT;\n".format(i,ppw[i]*x)
    if E:
        str+="\tn[{}]-={}*q;\n".format(N-1,x)
    return str

#propagate carries
def prop(n,base) :
    str="//propagate carries\n"
    str+="#[allow(unused_variables)]\n"
    str+="#[inline]\n"
    str+="fn prop(n: &mut [SPINT]) -> SPINT{\n"
    str+="\tlet mask=((1 as SPINT)<<{})-1;\n".format(base)
    str+="\tlet mut carry=(n[0] as SSPINT)>>{};\n".format(base)
    str+="\tn[0]&=mask;\n"
    str+="\tfor i in 1..{} {{\n".format(N-1)
    str+="\t\tcarry += n[i] as SSPINT;\n"
    str+="\t\tn[i] = (carry as SPINT) & mask;\n"
    str+="\t\tcarry>>={};\n".format(base)
    str+="\t}\n"
    str+="\tn[{}]+=carry as SPINT;\n".format(N-1)
    str+="\treturn ((n[{}] as SSPINT)>>{}) as SPINT;\n}}\n".format(N-1,WL-1)
    return str

#propagate carries and add p if negative, propagate carries again
def flat(n,base) :
    mask=(1<<base)-1
    str="//propagate carries and add p if negative, propagate carries again\n"
    str+="#[allow(unused_variables)]\n"
    str+="#[inline]\n"
    str+="fn flatten(n: &mut [SPINT]) -> bool {\n"
    str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    str+="\tlet carry=prop(n);\n"
    str+=caddp(1)
    str+="\tprop(n);\n"
    str+="\treturn (carry&1) == 1;\n}\n"
    return str

#final subtract
def modfsb(n,base) : 
    str="//Montgomery final subtract\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modfsb(n: &mut[SPINT]) -> bool {\n"
    str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    str+=subp(1)
    str+="\treturn flatten(n);\n}\n"
    return str

#modular addition
def modadd(n,base) :
    str="//Modular addition - reduce less than 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modadd(b: &[SPINT],n: &mut [SPINT]) {\n"

    if not algorithm :
        if E:
            str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    for i in range(0,N) :
        str+="\tn[{}]=n[{}]+b[{}];\n".format(i,i,i)
    if not algorithm :
        str+=subp(2)
        str+="\tlet carry=prop(n);\n"
        str+=caddp(2)
    str+="\tprop(n);\n"    
    str+="\treturn;\n}\n"
    return str

#modular subtraction
def modsub(n,base) :
    str="//Modular subtraction - reduce less than 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modsub(b: &[SPINT],n: &mut [SPINT]) {\n"
    if not algorithm :
        if E:
            str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    else :
        str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    for i in range(0,N) :
        str+="\tn[{}]=n[{}]-b[{}];\n".format(i,i,i)
    if not algorithm :
        str+="\tlet carry=prop(n);\n"
        str+=caddp(2)
    else :
        str+=addp(mp)
    str+="\tprop(n);\n"    
    str+="\treturn;\n}\n"
    return str

#modular negation
def modneg(n,base) :
    str="//Modular negation\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modneg(n: &mut [SPINT]) {\n"
    if not algorithm :
        if E:
            str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    else :
        str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    for i in range(0,N) :
        str+="\tn[{}]=(0 as SPINT)-n[{}];\n".format(i,i)
    if not algorithm :
        str+="\tlet carry=prop(n);\n"
        str+=caddp(2)
    else :
        str+=addp(mp)
    str+="\tprop(n);\n"    
    str+="\treturn;\n}\n"
    return str

# add column of partial products from multiplication on way up
def getZMU(str,i) :
    first=True
    global maxnum

    if karatsuba :
        if i==0 :
            str+="\tlet mut u=d0; t += u;"
            maxnum+=maxdigit*maxdigit
        else :
            str+="\tu+=d{}; t+=u;".format(i)
            for m in range(i,int(i/2),-1) :
                str+=" t+=(((c[{}] as SSPINT)-(c[{}] as SSPINT)) as DPINT) * (((b[{}] as SSPINT)-(b[{}] as SSPINT)) as DPINT); ".format(m,i - m, i - m, m)
                maxnum+=maxdigit*maxdigit
        return str

    k=0
    while (k<=i) :
        if first :
            str+="\tt+=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,i-k)
            first=False
        else :
            str+=" t+=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,i-k)
        k+=1
        maxnum+=maxdigit*maxdigit
    return str

# add column of partial products from multiplication on way down
def getZMD(str,i) :

    if karatsuba :
        str+="\tu-=d{}; t+=u; ".format(i - N)
        for m in range(N-1,int(i/2),-1) :
            str+="t+=(((c[{}] as SSPINT)-(c[{}] as SSPINT)) as DPINT) * (((b[{}] as SSPINT)-(b[{}] as SSPINT)) as DPINT); ".format(m, i - m, i - m, m)
        return str

    first=True
    k=i-(N-1)
    while (k<=N-1) :
        if first :
            first=False
            str+="\tt+=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,i-k)
        else :
            str+=" t+=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,i-k)
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
            str+="\ttot=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,i-k)
            first=False
        else :
            str+=" tot+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,i-k)
        k+=1
        j-=1
    if hap:
        str+=" tot*=2;"
    if i%2==0:
        if first :
            str+="\ttot=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(int(i/2),int(i/2))
        else :
            str+=" tot+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(int(i/2),int(i/2))
    if i==0 :
        str+=" t=tot;"
    else :
        str+=" t+=tot; "
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
            str+="\ttot=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,i-k)
            first=False
        else :
            str+=" tot+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,i-k)
        k+=1
        j-=1
    if hap :
        str+=" tot*=2;"
    if i%2==0:
        if first :
            str+="\ttot=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(int(i/2),int(i/2))
        else :
            str+=" tot+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(int(i/2),int(i/2))
    str+=" t+=tot; "
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
            str+=" t+=(v{} as DPINT)<<{}; ".format(j,e)
            maxnum+=maxdigit*2**e
        else :
            str+=" t+=(v{} as DPINT)*(p{} as DPINT); ".format(j,i)
            maxnum+=maxdigit*ppw[i]
    if n == 1 :
        if mask_set :
            str+=" s+=v{}; ".format(j)
        else :
            str+=" t+=v{} as DPINT; ".format(j)
            maxnum+=maxdigit
    if n == -1 :
        if mask_set :
            if not gone_neg :
                str+=" s+=q-v{};".format(j)
            else :
                str+=" s-=v{}; ".format(j)
        else: 
            if not gone_neg :
                str+=" t+=(q-v{}) as DPINT;".format(j)
                maxnum+=maxdigit
            else :
                str+=" t-=v{} as DPINT; ".format(j)
        gone_neg=True
    return str,gone_neg

def sqr_process(i,j,str,gone_neg,mask_set) :
    global maxnum
    n=ppw[i]   # do nothing if n=0
    if abs(n)>1 :
        e=ispowerof2(n)
        if e > 0  :
            str+=" t+=(v{} as UDPINT)<<{}; ".format(j,e)
            maxnum+=maxdigit*2**e
        else :
            str+=" t+=(v{} as UDPINT)*(p{} as UDPINT); ".format(j,i)
            maxnum+=maxdigit*ppw[i]
    if n == 1 :
        if mask_set :
            str+=" s+=v{}; ".format(j)
        else :
            str+=" t+=v{} as UDPINT; ".format(j)
            maxnum+=maxdigit
    if n == -1 :
        if mask_set :
            if not gone_neg :
                str+=" s+=q-v{};".format(j)
            else :
                str+=" s-=v{}; ".format(j)
        else: 
            if not gone_neg :
                str+=" t+=(q-v{}) as UDPINT;".format(j)
                maxnum+=maxdigit
            else :
                str+=" t-=v{} as UDPINT; ".format(j)
        gone_neg=True
    return str,gone_neg


def modmli(n,base) :
    N=getN(n,base)
    mask=(1<<base)-1
    xcess=N*base-n
    str="// Modular multiplication by an integer, c=c*b mod 2p\n"
    str+="#[allow(dead_code)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modmli(b: usize,c: &mut [SPINT]) {\n"
    if trin>0 :
        str+="\tlet mut t=0 as UDPINT;\n"
        str+="\tlet mask=((1 as SPINT)<<{})-1;\n".format(base)

        for i in range(0,N) :
            str+="\tt+=(c[{}] as UDPINT)*(b as UDPINT); ".format(i)
            str+="c[{}]=(t as SPINT) & mask; t=t>>{};\n".format(i,base)

        str+="// reduction pass\n\n"  
        str+="\tlet s=t as SPINT;\n"  

        if xcess>0 :
            smask=(1<<(base-xcess))-1
            str+= "\ts=(s<<{})+(c[{}]>>{}); c[{}]&=0x{:x};\n".format(xcess,N-1,base-xcess,N-1,smask)

        str+="\tc[0]+=s;\n"
        str+="\tc[{}]+=s;\n".format(trin)
    else :

# Barrett/Dhem method
#        str+="/*\n"
        for i in range(0,N) :
            d=ppw[i]
            if abs(d)>1 :
                if d>1 :
                    if i==0 or ispowerof2(d)<0 :
                        str+="\tlet p{}={} as SPINT;\n".format(i,hex(d))
                if d<1 :
                    if i==0 :
                        str+="\tlet p{}={} as SPINT;\n".format(i,hex(-d))

        str+="\tlet mask=((1 as SPINT)<<{})-1;\n".format(base)
        str+="\tlet mut t=0 as UDPINT;\n"
        str+="\tlet r=0x{:x} as SPINT;\n".format((2**(n+base))//p)
        for i in range(0,N-1) :
            str+="\tt+=(c[{}] as UDPINT)*(b as UDPINT); ".format(i)
            str+="c[{}]=(t as SPINT) & mask; t=t>>{};\n".format(i,base)
        str+="\tt+=(c[{}] as UDPINT)*(b as UDPINT); ".format(N-1)
        str+="c[{}]=t as SPINT;\n".format(N-1)
        str+="\t\n//Barrett-Dhem reduction\n"
        str+="\tlet h = (t>>{}) as SPINT;\n".format((n-WL)%base)
        str+="\tlet q=(((h as UDPINT)*(r as UDPINT))>>{}) as SPINT;\n".format(WL)  

# required to propagate carries?
        propc=False
        for i in range(1,N-1) :
            if ppw[i]!=0 :
                propc=True
        if ppw[0]>0 :
            propc=True

        for i in range(0,N) :
            d=ppw[i]
            if d==0 :
                continue
            if d==-1 :
                str+="\tc[{}]+=q;\n".format(i)
                continue
            if d==1 :
                str+="\tc[{}]-=q;\n".format(i)
            e=ispowerof2(d)
            if e>0 :
                if i<N-1 :
                    str+="\tt=(q as UDPINT)<<{}; c[{}]-=(t as SPINT)&mask; c[{}]-=(t>>{}) as SPINT;\n".format(e,i,i+1,base)
                else :
                    str+="\tc[{}]-=q<<{};\n".format(i,e)
            else :
                if d<0 :
                    str+="\tt=(q as UDPINT)*(p{} as UDPINT); c[{}]+=(t as SPINT)&mask; c[{}]+=(t>>{}) as SPINT;\n".format(i,i,i+1,base)
                else :
                    if i<N-1 :
                        str+="\tt=(q as UDPINT)*(p{} as UDPINT); c[{}]-=(t as SPINT)&mask; c[{}]-=(t>>{}) as SPINT;\n".format(i,i,i+1,base)
                    else :
                        str+="\tc[{}]-=q*p{};\n".format(i,i)
                        #str+="\tc[{}]=(c[{}]-(q*p{}))&mask;\n".format(i,i,i)
        if propc :
            str+="\tprop(c);\n"
        else :
            if fully_propagate :
                str+="\tprop(c);\n"
#        str+="*/\n"
#        str+="\tlet mut t: [SPINT; {}] = [0; {}];\n".format(N,N)
#        str+="\tmodint(b,&mut t);\n"
#        str+="\tmodmul(&t,c);\n"

    str+="\treturn;\n}\n"
    return str

#def dummymodmli(n,base) :
#    str="// Modular multiplication by an integer, c=c*b mod 2p\n"
#    str+="#[allow(dead_code)]\n"
#    if inline :
#        str+="#[inline]\n"
#    if makepublic :
#        str+="pub "
#    str+="fn modmli(_b: usize,_c: &mut [SPINT]) {\n"
#    str+="\treturn;\n}\n"
#    return str

# modular multiplication, modulo p. Exploits 0 digits in p.
# Note that allowing inlining gives significant speed-up
def modmul(n,base) :
    global maxnum
    mask=(1<<base)-1

    str="// Modular multiplication, c=c*b mod 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modmul(b: &[SPINT],c: &mut [SPINT]) {\n"

    str+="\tlet mut t = 0 as DPINT;\n"
    str+="\tlet mut s:SPINT;\n"
    if karatsuba :
        for i in range(0,N) :
            str+="\tlet d{}=(c[{}] as DPINT)*(b[{}] as DPINT);\n".format(i, i, i)

    for i in range(0,N) :
        if i==0 and PM :
            continue
        n=ppw[i]
        if abs(n)>1 :
            if i==0 or ispowerof2(n)<0 :
                str+="\tlet p{}={} as SPINT;\n".format(i,hex(ppw[i]))
    str+="\tlet q=(1 as SPINT)<<{}; // q is unsaturated radix \n".format(base)
    str+="\tlet mask=(q-1) as SPINT;\n"
    if fullmonty :
        str+="\tlet ndash=0x{:x} as SPINT;\n".format(ndash)
    maxnum=0
    str=getZMU(str,0)
    gone_neg=False 
    
    if fullmonty :
        str+=" let v0=(((t as SPINT)*ndash)&mask) as SPINT;"
        if PM :
            gone_neg=True
            str+=" t+=({}*(q-v0)) as DPINT; ".format(M)
            maxnum+=M*maxdigit
        else :
            if ppw[0]==1 :
                str+=" t+=v0 as DPINT;"
            else :
                str+=" t+=(v0 as DPINT) * (p0 as DPINT);"
            maxnum+=ppw[0]*maxdigit
    else :
        str+=" let v0=((t as SPINT) & mask) as SPINT;"
    str+=" t>>={};\n".format(base)
    maxnum=2**(2*WL-base)
    str=getZMU(str,1)
    mask_set=False

    for i in range(1,N) :
        if gone_neg :
            if PM :
                str+=" t+=({}*mask) as DPINT;".format(M)
                maxnum+=M*maxdigit
            else :
                str+=" s=mask as SPINT;"
                mask_set=True
                #str+=" t+=mask;"
        k=1
        str,gone_neg=mul_process(i,0,str,gone_neg,mask_set)
        while k<i :
            str,gone_neg=mul_process(i-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" t+=s as DPINT;"
            maxnum+=maxdigit
            mask_set=False

        if fullmonty :
            str+=" let v{}=(((t as SPINT)*ndash) & mask) as SPINT; ".format(i)
            if PM :
                str+=" t-=({}*v{}) as DPINT; ".format(M,i)
            else :
                if ppw[0]==1 :
                    str+=" t+=v{} as DPINT; ".format(i)
                else :
                    str+=" t+=(v{} as DPINT)*(p0 as DPINT); ".format(i)
                maxnum+=ppw[0]*maxdigit
        else :
            str+=" let v{}=((t as SPINT) & mask) as SPINT; ".format(i)
        if i==N-1 :
            if karatsuba :
                print("// Run using Comba in modmul for tighter overflow check ")
            else :
                print("// Overflow limit   =",2**(2*WL))
                print("// maximum possible =",maxnum)
                if maxnum >= 2**(2*WL) :
                    print("//Warning: Overflow possibility detected - change radix ")
        str+=" t>>={};\n".format(base)
        maxnum=2**(2*WL-base)
        if i<N-1 :
            str= getZMU(str,i+1)
     
    str=getZMD(str,N)
    if gone_neg :
        if PM :
            str+=" t+=({}*mask) as DPINT;".format(M)
        else :
            str+=" s=mask as SPINT;"
            mask_set=True
            #str+=" t+=mask;"

    if E :  
        k=0
        while k<N :
            str,gone_neg=mul_process(N-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" t+=s as DPINT;"
            mask_set=False        

        if fullmonty :
            str+=" let v{}=(((t as SPINT)*ndash) & mask) as SPINT; ".format(N)
            if PM :
                str+=" t-=({}*v{}) as DPINT; ".format(M,N)
            else :
                if ppw[0]==1 :
                    str+=" t+=v{} as DPINT; ".format(N)
                else :
                    str+=" t+=(v{} as DPINT)*(p0 as DPINT); ".format(N)
        else :
            str+=" let v{}=((t as SPINT) & mask) as SPINT; ".format(N)
        str+=" t>>={};\n".format(base)

        str=getZMD(str,N+1)

        for i in range(N+1,2*N) :
            if gone_neg :
                if PM :
                    str+=" t+=({}*mask) as DPINT;".format(M)
                else :
                    str+=" s=mask as SPINT;"
                    mask_set=True
            k=i-N 
            while k<=N : 
                str,gone_neg=mul_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" t+=s as DPINT;"
                mask_set=False        

            str+=" c[{}]=((t as SPINT) & mask) as SPINT; ".format(i-N-1)
            str+=" t>>={};\n".format(base)

            #if i<=2*N-1 :
            if i<2*N-2 :
                str=getZMD(str,i+1)
            else :
                str+="\t"
        
        if gone_neg :
            if PM :
                str+=" t+=(v{}-{}) as DPINT;".format(N,M)
            else :
                str+=" t+=(v{}-1) as DPINT;".format(N)
        else :
            str+=" t+=v{} as DPINT;".format(N)

        str+=" c[{}] = t as SPINT;\n".format(N-1)
    else :
        for i in range(N,2*N-1) :
            k=i-(N-1) 
            while k<=N-1 :
                str,gone_neg=mul_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" t+=s as DPINT;"
                mask_set=False   

            if i==2*N-1 :
                break
            str+=" c[{}]=((t as SPINT) & mask) as SPINT; ".format(i-N)
            str+=" t>>={};\n".format(base)
            if i<=2*N-3 :
                str=getZMD(str,i+1)
                if gone_neg :
                    if PM :
                        str+=" t+=({}*mask) as DPINT;".format(M)
                    else :
                        str+=" s=mask as SPINT;"
                        mask_set=True
        if gone_neg :
            if PM :
                str+="\tt-={};".format(M)
            else :
                str+="\tt-=1;"
        str+="\tc[{}] = t as SPINT;\n".format(N-1)
    str+="\treturn;\n}\n"
    return str

# modular squaring
# Note that allowing inlining gives significant speed-up
def modsqr(n,base) :
    mask=(1<<base)-1

    str="// Modular squaring, c=a*a  mod 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modsqr(c: &mut [SPINT]) {\n"

    str+="\tlet mut t:UDPINT;\n"
    str+="\tlet mut tot:UDPINT;\n"
    str+="\tlet mut s:SPINT;\n"
    for i in range(0,N) :
        if i==0 and PM :
            continue
        n=ppw[i]
        if abs(n)>1 :
            if i==0 or ispowerof2(n)<0 :
                str+="\tlet p{}={} as SPINT;\n".format(i,hex(ppw[i]))
    str+="\tlet q=(1 as SPINT)<<{}; // q is unsaturated radix \n".format(base)
    str+="\tlet mask=(q-1) as SPINT;\n"
    if fullmonty :
        str+="\tlet ndash=0x{:x} as SPINT;\n".format(ndash)
    str=getZSU(str,0)
    gone_neg=False

    if fullmonty :
        str+=" let v0=(((t as SPINT)*ndash)& mask) as SPINT;"
        if PM :
            gone_neg=True
            str+=" t+=({}*(q-v0)) as UDPINT; ".format(M)
        else :
            if ppw[0]==1 :
                str+=" t+=v0 as UDPINT;"
            else :
                str+=" t+=(v0 as UDPINT) * (p0 as UDPINT);"

    else :
        str+=" let v0=((t as SPINT) & mask) as SPINT;"
    str+=" t>>={};\n".format(base)
    
    str=getZSU(str,1)
    mask_set=False

    for i in range(1,N) :
        if gone_neg :
            if PM :
                str+=" t+=({}*mask) as UDPINT;".format(M)
            else :
                str+=" s=mask as SPINT;"
                mask_set=True

        k=1
        str,gone_neg=sqr_process(i,0,str,gone_neg,mask_set)
        while k<i :
            str,gone_neg=sqr_process(i-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" t+=s as UDPINT;"
            mask_set=False   

        if fullmonty :
            str+=" let v{}=(((t as SPINT)*ndash) & mask) as SPINT; ".format(i)
            if PM :
                str+=" t-=({}*v{}) as UDPINT; ".format(M,i)
            else :
                if ppw[0]==1 :
                    str+=" t+=v{} as UDPINT; ".format(i)
                else :
                    str+=" t+=(v{} as UDPINT)*(p0 as UDPINT); ".format(i)
        else :
            str+=" let v{}=((t as SPINT) & mask) as SPINT;".format(i)
        str+=" t>>={};\n".format(base)
        if i<N-1 :
            str= getZSU(str,i+1)
     
    str=getZSD(str,N)
    if gone_neg :
        if PM :
            str+=" t+=({}*mask) as UDPINT;".format(M)
        else :
            str+=" s=mask as SPINT;"
            mask_set=True

    if E :  
        k=0
        while k<N :
            str,gone_neg=sqr_process(N-k,k,str,gone_neg,mask_set)
            k+=1
        if mask_set :
            str+=" t+=s as UDPINT;"
            mask_set=False   
        if fullmonty :
            str+=" let v{}=(((t as SPINT)*ndash) & mask) as SPINT; ".format(N)
            if PM :
                str+=" t-=({}*v{}) as UDPINT; ".format(M,N)
            else :
                if ppw[0]==1 :
                    str+=" t+=v{} as UDPINT; ".format(N)
                else :
                    str+=" t+=(v{} as UDPINT)*(p0 as UDPINT); ".format(N)
        else :
            str+=" let v{}=((t as SPINT) & mask) as SPINT; ".format(N)
        str+=" t>>={};\n".format(base)

        str=getZSD(str,N+1)
        for i in range(N+1,2*N) :
            if gone_neg :
                if PM :
                    str+=" t+=({}*mask) as UDPINT;".format(M)
                else :
                    str+= "s=mask as SPINT;"
                    mask_set=True
                    #str+=" t+=mask;"
            k=i-N 
            while k<=N :
                str,gone_neg=sqr_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" t+=s as UDPINT;"
                mask_set=False   
            str+=" c[{}]=((t as SPINT) & mask) as SPINT; ".format(i-N-1)
            str+=" t>>={};\n".format(base)

            #if i<=2*N-1 :
            if i<2*N-2 :
                str=getZSD(str,i+1)
            else :
                str+="\t"

        if gone_neg :
            if PM :
                str+=" t+=(v{}-{}) as UDPINT;".format(N,M)
            else :
                str+=" t+=(v{}-1) as UDPINT;".format(N)
        else :
            str+=" t+=v{} as UDPINT;".format(N)
        str+=" c[{}] = t as SPINT;\n".format(N-1)
    else :
        for i in range(N,2*N-1) :
            k=i-(N-1) 
            while k<=N-1 :
                str,gone_neg=sqr_process(i-k,k,str,gone_neg,mask_set)
                k+=1
            if mask_set :
                str+=" t+=s as UDPINT;"
                mask_set=False   
            if i==2*N-1 :
                break
            str+=" c[{}]=((t as SPINT) & mask) as SPINT; ".format(i-N)
            str+=" t>>={};\n".format(base)
            if i<=2*N-3 :
                str=getZSD(str,i+1)
                if gone_neg :
                    if PM :
                        str+=" t+=({}*mask) as UDPINT;".format(M)
                    else :
                        str+= "s=mask as SPINT;"
                        mask_set=True
        if gone_neg :
            if PM :
                str+="\tt-={};".format(M)
            else :
                str+="\tt-=1;"
        str+="\tc[{}] = t as SPINT;\n".format(N-1)
    str+="\treturn;\n}\n"
    return str

def modcpy() :
    str="//copy\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modcpy(a: &[SPINT],c: &mut [SPINT]) {\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\tc[i]=a[i];\n"
    str+="\t}\n"
    str+="\treturn;\n}\n"
    return str

def modnsqr() :
    str="//square n times\n"
    if makepublic :
        str+="pub "
    str+="fn modnsqr(a:&mut [SPINT],n: isize) {\n"
    str+="\tfor _i in 0..n {\n"
    str+="\t\tmodsqr(a);\n"
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
    if makepublic :
        str+="pub "
    str+="fn modpro(w: &[SPINT],r: &mut [SPINT]) {\n"
    str+="\tlet mut x: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut z: [SPINT; {}] = [0; {}];\n".format(N,N)
    for i in range(0,ntmps) :
        str+="\tlet mut t{}: [SPINT; {}] = [0; {}];\n".format(i,N,N)   
    str+="\tmodcpy(w,&mut x);\n"
    for i in range(1,len(lines)) :
        info=lines[i].split()
        if info[0]=="double" :
            if info[1]!=info[2] :
                str+="\tmodcpy(&{},&mut {});\n".format(info[2],info[1])
            str+="\tmodsqr(&mut {});\n".format(info[1])
        if info[0]=="add" :
            if info[2]==info[1] :
                str+="\tmodmul(&{},&mut {});\n".format(info[3],info[1])
            if info[3]==info[1] :
                str+="\tmodmul(&{},&mut {});\n".format(info[2],info[1])
            if info[2]!=info[1] and info[3]!=info[1] :
                str+="\tmodcpy(&{},&mut {});\n".format(info[2],info[1])
                str+="\tmodmul(&{},&mut {});\n".format(info[3],info[1])
        if info[0]=="shift" :
            if info[2]!=info[1] :
                str+="\tmodcpy(&{},&mut {});\n".format(info[2],info[1])
            str+="\tmodnsqr(&mut {},{});\n".format(info[1],int(info[3]))
    str+="\tmodcpy(&z,r);\n"
    str+="\treturn;\n}\n"
    subprocess.call("rm ac.txt",shell=True)    
    return str

def modinv() :
    str="//calculate inverse, provide progenitor h if available\n"
    if makepublic :
        str+="pub "
    str+="fn modinv(h: Option<&[SPINT]>,z: &mut[SPINT]) {\n"
    str+="\tlet mut s: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut t: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tif let Some(hint) = h {\n"
    str+="\t\tmodcpy(&hint,&mut t);\n"
    str+="\t} else {\n"
    str+="\t\tmodpro(&z,&mut t);\n"
    str+="\t}\n"
    str+="\tmodcpy(&z,&mut s);\n"
    if PM1D2>1 :
        str+="\tfor _i in 0..{}-1 {{\n".format(PM1D2)
        str+="\t\tmodsqr(&mut s);\n"
        str+="\t\tmodmul(&z,&mut s);\n"
        str+="\t}\n"
    str+="\tmodnsqr(&mut t,{});\n".format(PM1D2+1)
    str+="\tmodcpy(&t,z);\n"
    str+="\tmodmul(&s,z);\n"
    str+="\treturn;\n}\n"
    return str

def modqr() :
    str="//Test for quadratic residue \n"
    if makepublic :
        str+="pub "
    str+="fn modqr(h: Option<&[SPINT]>,x: &[SPINT]) -> bool {\n"
    str+="\tlet mut r: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tif let Some(hint) = h {\n"
    str+="\t\tmodcpy(&hint,&mut r);\n"
    str+="\t} else {\n"
    str+="\t\tmodpro(&x,&mut r);\n"
    str+="\t}\n" 
    str+="\tmodsqr(&mut r);\n"
    str+="\tmodmul(&x,&mut r);\n"
    if PM1D2>1 :
        str+="\tmodnsqr(&mut r,{});\n".format(PM1D2-1)
    str+="\treturn modis1(&r) || modis0(x) ;\n}\n"
    return str

def modsqrt() :
    str="//Modular square root, provide progenitor h if available, NULL if not\n"
    if makepublic :
        str+="pub "
    str+="fn modsqrt(x: &[SPINT],h: Option<&[SPINT]>,r: &mut[SPINT]) {\n"
    if PM1D2>1 :
        str+="\tlet mut t: [SPINT; {}] = [0; {}];\n".format(N,N)
        str+="\tlet mut b: [SPINT; {}] = [0; {}];\n".format(N,N)
        str+="\tlet mut v: [SPINT; {}] = [0; {}];\n".format(N,N)
        str+="\tlet mut z: [SPINT; {}] = [".format(N)
        for i in range(0,N-1) :
            str+=hex(ROI[i])
            str+=","
        str+=hex(ROI[N-1])
        str+="];\n"
    str+="\tlet mut y: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut s: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tif let Some(hint) = h {\n"
    str+="\t\tmodcpy(&hint,&mut y);\n"
    str+="\t} else {\n"
    str+="\t\tmodpro(&x,&mut y);\n"
    str+="\t}\n"
    str+="\tmodcpy(&x,&mut s);\n"
    str+="\tmodmul(&y,&mut s);\n"
    if PM1D2>1 :
        str+="\tmodcpy(&s,&mut t);\n"
        str+="\tmodmul(&y,&mut t);\n"
        str+="\tnres(&mut z);\n"
        str+="\tfor k in (2..={}).rev() {{\n".format(PM1D2)
        str+="\t\tmodcpy(&t,&mut b);\n"

        str+="\t\tmodnsqr(&mut b,k-2);\n"
        str+="\t\tlet d= 1 - (modis1(&b) as usize);\n"
        str+="\t\tmodcpy(&z,&mut v);\n"
        str+="\t\tmodmul(&s,&mut v);\n"
        str+="\t\tmodcmv(d,&v,&mut s);\n"
        str+="\t\tmodsqr(&mut z);\n"
        str+="\t\tmodcpy(&z,&mut v);\n"
        str+="\t\tmodmul(&t,&mut v);\n"
        str+="\t\tmodcmv(d,&v,&mut t);\n"
        str+="\t}\n" 
    str+="\tmodcpy(&s,r);\n"
    str+="\treturn;\n}\n"
    return str

def modis1(n,base) :
    str="//is unity?\n"
    if makepublic :
        str+="pub "
    str+="fn modis1(a: &[SPINT]) -> bool {\n"
    str+="\tlet mut c: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut d=0 as SSPINT;\n"
    str+="\tmodcpy(a,&mut c);\n"
    str+="\tredc(&mut c);\n"
    str+="\tfor i in 1..{} {{\n".format(N)
    str+="\t\td|=c[i] as SSPINT;\n\t}\n"
    str+="\tlet c0=c[0] as SSPINT;\n"
    str+="\treturn (1 & ((d-1)>>{}) & (((c0^1)-1)>>{})) != 0;\n}}\n".format(base,base)
    return str

def modis0(n,base) :
    str="//is zero?\n"
    if makepublic :
        str+="pub "
    str+="fn modis0(a: &[SPINT]) -> bool {\n"
    str+="\tlet mut c: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut d=0 as SSPINT;\n"
    str+="\tmodcpy(a,&mut c);\n"
    str+="\tredc(&mut c);\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\td|=c[i] as SSPINT;\n\t}\n"
    str+="\treturn (1 & ((d-1)>>{})) != 0;\n}}\n".format(base)
    return str    

def modzer() :
    str="//set to zero\n"
    if makepublic :
        str+="pub "
    str+="fn modzer(a: &mut[SPINT]) {\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\ta[i]=0;\n\t}\n"
    str+="\treturn;\n}\n"
    return str

def modone() :
    str="//set to one\n"
    if makepublic :
        str+="pub "
    str+="fn modone(a: &mut[SPINT]) {\n"
    str+="\ta[0]=1;\n"
    str+="\tfor i in 1..{} {{\n".format(N)
    str+="\t\ta[i]=0;\n\t}\n"
    str+="\tnres(a);\n"
    str+="\treturn;\n}\n"
    return str

def modint() :
    str="//set to integer\n"
    if makepublic :
        str+="pub "
    str+="fn modint(x: usize,a: &mut[SPINT]) {\n"
    str+="\ta[0]=x as SPINT;\n"
    str+="\tfor i in 1..{} {{\n".format(N)
    str+="\t\ta[i]=0;\n\t}\n"
    str+="\tnres(a);\n"
    str+="\treturn;\n}\n"
    return str

# convert to Montgomery nresidue form
def nres(n,base) :
    str="//Convert n to n-residue form, n=nres(m) \n"
    if makepublic :
        str+="pub "
    str+="fn nres(n: &mut [SPINT]) {\n"
    str+="\tlet c: [SPINT; {}] = [".format(N)
    for i in range(0,N-1) :
        str+=hex(cw[i])
        str+=","
    str+=hex(cw[N-1])
    str+="];\n"
    str+="\tmodmul(&c,n);\n"
    str+="\treturn;\n}\n"
    return str 

#convert back from Montgomery to integer form
def redc(n,base) :
    str="//Convert m back to normal form, m=redc(n) \n"
    if makepublic :
        str+="pub "
    str+="fn redc(m: &mut [SPINT]) {\n"
    str+="\tlet mut c:[SPINT;{}]=[0;{}];\n".format(N,N)
    str+="\tc[0]=1;\n"
    str+="\tmodmul(&c,m);\n"
    str+="\tmodfsb(m);\n"
    str+="\treturn;\n}\n"
    return str 

#conditional swap - see Loiseau et al. 2021
def modcsw() :
    str="//conditional swap g and f if d=1\n"
    str+="#[inline(never)]\n"
    if makepublic :
        str+="pub "
    str+="fn modcsw(b: usize,g: &mut [SPINT],f: &mut [SPINT]) {\n"
    str+="\tlet bb = b as SPINT;\n"
    str+="\tstatic mut R:SPINT=0;\n"
    str+="\tlet w:SPINT;\n"
    str+="\tunsafe {\n"
    if WL==16 :
        str+="\t\tR+=0xa55a;\n"
    if WL==32 :
        str+="\t\tR+=0x5aa5a55a;\n"
    if WL==64 :
        str+="\t\tR+=0x3cc3c33c5aa5a55a;\n"
    str+="\t\tw=R;\n"
    str+="}\n"
    str+="\t\tlet c0=(!bb)&(w+1);\n"
    str+="\t\tlet c1=bb+w;\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\tlet s=g[i];\n"
    str+="\t\tlet t=f[i];\n"
    str+="\t\tlet v=w*(t+s);\n"
    str+="\t\tunsafe{core::ptr::write_volatile(&mut f[i],c0*t+c1*s);}\n"
    str+="\t\tf[i]-=v;\n"
    str+="\t\tunsafe{core::ptr::write_volatile(&mut g[i],c0*s+c1*t);}\n"
    str+="\t\tg[i]-=v;\n\t}\n"
    str+="\treturn;\n}\n"
    return str

#conditional move
def modcmv() :
    str="//conditional move g to f if d=1\n"
    str+="#[inline(never)]\n"
    if makepublic :
        str+="pub "
    str+="fn modcmv(b: usize,g: &[SPINT],f: &mut [SPINT]) {\n"
    str+="\tstatic mut R:SPINT=0;\n"
    str+="\tlet w:SPINT;\n"
    str+="\tunsafe {\n"
    if WL==16 :
        str+="\t\tR+=0xa55a;\n"
    if WL==32 :
        str+="\t\tR+=0x5aa5a55a;\n"
    if WL==64 :
        str+="\t\tR+=0x3cc3c33c5aa5a55a;\n"
    str+="\t\tw=R;\n"
    str+="}\n"
    str+="\tlet bb = b as SPINT;\n"
    str+="\t\tlet c0=(!bb)&(w+1);\n"
    str+="\t\tlet c1=bb+w;\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\tlet s=g[i];\n"
    str+="\t\tlet t=f[i];\n"
    str+="\t\tunsafe{core::ptr::write_volatile(&mut f[i],c0*t+c1*s);}\n"
    str+="\t\tf[i]-=w*(t+s);\n\t}\n"
    str+="\treturn;\n}\n"
    return str

#shift left
def modshl(n,base) :
    N=getN(n,base)
    mask=(1<<base)-1
    str="//shift left by less than a word\n"
    if makepublic :
        str+="pub "
    str+="fn modshl(n: isize,a: &mut [SPINT]) {\n"
    str+="\ta[{}-1]=((a[{}-1]<<n)) | (a[{}-2]>>({}-n));\n".format(N,N,N,base)
    str+="\tfor i in (1..{}-1).rev() {{\n".format(N)
    str+="\t\ta[i]=((a[i]<<n)&0x{:x}) | (a[i-1]>>({}-n));\n\t}}\n".format(mask,base)
    str+="\ta[0]=(a[0]<<n)&0x{:x};\n".format(mask)
    str+="\treturn;\n}\n"
    return str 

#shift right
def modshr(n,base) :
    N=getN(n,base)
    mask=(1<<base)-1
    str="//shift right by less than a word. Return shifted out part\n"
    if makepublic :
        str+="pub "
    str+="fn modshr(n: isize,a: &mut [SPINT]) -> isize {\n"
    str+="\tlet r=a[0]&((1<<n)-1);\n"
    str+="\tfor i in 0..{}-1 {{\n".format(N)
    str+="\t\ta[i]=(a[i]>>n) | ((a[i+1]<<({}-n))&0x{:x});\n\t}}\n".format(base,mask)
    str+="\ta[{}-1]=a[{}-1]>>n;\n".format(N,N)
    str+="\treturn r as isize;\n}\n"
    return str

#set a=2^r
def mod2r() :
    str="//set a= 2^r\n"
    if makepublic :
        str+="pub "
    str+="fn mod2r(r:usize, a: &mut [SPINT]) {\n"
    str+="\tlet n=r/{};\n".format(base)
    str+="\tlet m=r%{};\n".format(base)    
    str+="\tmodzer(a);\n"
    str+="\tif r>={}*8 {{\n".format(Nbytes)
    str+="\t\treturn;\n\t}\n"
    str+="\ta[n]=1;\n"
    str+="\ta[n] <<= m;\n"
    str+="\tnres(a);\n}\n"
    return str

#export to byte array
def modexp() :
    str="//export to byte array\n"
    if makepublic :
        str+="pub "
    str+="fn modexp(a: &[SPINT],b: &mut [u8]) {\n"
    str+="\tlet mut c: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tmodcpy(a,&mut c);\n"
    str+="\tredc(&mut c);\n"
    str+="\tfor i in (0..{}).rev() {{\n".format(Nbytes)
    str+="\t\tb[i]=(c[0]&0xff) as u8;\n"
    str+="\t\tmodshr(8,&mut c);\n\t}\n"
    str+="\treturn;\n}\n"
    return str 

#import from byte array
def modimp() :
    str="//import from byte array\n"
    if makepublic :
        str+="pub "
    str+="fn modimp(b: &[u8], a: &mut [SPINT]) -> bool {\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\ta[i]=0;\n\t}\n"
    str+="\tfor i in 0..{} {{\n".format(Nbytes)
    str+="\t\tmodshl(8,a);\n"
    str+="\t\ta[0]+=b[i] as SPINT;\n\t}\n"
    str+="\tlet res=modfsb(a);\n"
    str+="\tnres(a);\n"
    str+="\treturn res;\n}\n"
    return str 

#get sign (parity)
def modsign() :
    str="//determine sign\n"
    if makepublic :
        str+="pub "
    str+="fn modsign(a: &[SPINT]) -> usize {\n"
    str+="\tlet mut c:[SPINT;{}]=[0;{}];\n".format(N,N)
    str+="\tmodcpy(a,&mut c);\n"
    str+="\tredc(&mut c);\n"
    str+="\treturn (c[0]%2) as usize;\n}\n"
    return str 

#compare for equality
def modcmp() :
    str="//return true if equal\n"
    if makepublic :
        str+="pub "
    str+="fn modcmp(a: &[SPINT],b: &[SPINT]) -> bool {\n"
    str+="\tlet mut c:[SPINT;{}]=[0;{}];\n".format(N,N)
    str+="\tlet mut d:[SPINT;{}]=[0;{}];\n".format(N,N)
    str+="\tlet mut eq=1;\n"
    str+="\tmodcpy(a,&mut c); redc(&mut c);\n"
    str+="\tmodcpy(b,&mut d); redc(&mut d);\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\teq&=(((c[i]^d[i])-1)>>{})&1;\n\t}}\n".format(base)
    str+="\treturn eq==1;\n"
    str+="}\n"
    return str 

# for timings
def time_modmul(n,base,ra,rb) :
    N=getN(n,base)
    mask=(1<<base)-1
    bit=n%base
    smask=(1<<bit)-1
    rap=makebig(ra,base,N)
    rbp=makebig(rb,base,N)
    str="pub fn time_modmul() {\n"
    str+="\tlet mut x: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut y: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut z: [SPINT; {}] = [0; {}];\n".format(N,N)

    str+="\t"
    for i in range(0,N) :
        str+="x[{}]={}; ".format(i,hex(rap[i]))
    str+="\n\t"
    for i in range(0,N) :
        str+="y[{}]={}; ".format(i,hex(rbp[i]))
    str+="\n"

    str+="\tnres(&mut x);\n"
    str+="\tnres(&mut y);\n"
    if cyclesorsecs:
        str+="\tunsafe {\n"
        str+="\tlet pre = std::arch::x86_64::_rdtsc();\n"
    str+="\tlet begin = Instant::now();\n"
    str+="\tfor _i in 0..{} {{\n".format(100000//scale)
    str+="\t\tfor _j in 0..200 {\n"
    str+="\t\t\tmodcpy(&x,&mut z);\n"
    str+="\t\t\tmodmul(&y,&mut z);\n"
    str+="\t\t\tmodcpy(&z,&mut y);\n"
    str+="\t\t\tmodmul(&x,&mut y);\n"    
    str+="\t\t\tmodcpy(&y,&mut x);\n"
    str+="\t\t\tmodmul(&z,&mut x);\n"    
    str+="\t\t\tmodcpy(&x,&mut z);\n"
    str+="\t\t\tmodmul(&y,&mut z);\n"    
    str+="\t\t\tmodcpy(&z,&mut y);\n"
    str+="\t\t\tmodmul(&x,&mut y);\n"    
    str+="\t\t}\n"
    str+="\t}\n"
    if cyclesorsecs:
        str+="\tlet post = std::arch::x86_64::_rdtsc();\n"
    str+="\tlet elapsed = begin.elapsed();\n"
    str+="\tlet dur = (elapsed.as_secs() * 1_000_000_000) + (elapsed.subsec_nanos()) as u64;\n"
    str+="\tredc(&mut z);\n"
    if cyclesorsecs:
        str+='\tprintln!("modmul check 0x{{:x}} Clock cycles= {{}} Nanosecs= {{}}",(z[0] as i32)&0xFFFFFF,((post-pre)/{}),dur/{});\n'.format(100000000//scale,100000000//scale)
        str+="\t}\n"
    else :
        str+='\tprintln!("modmul check 0x{{:x}} Nanosecs= {{}}",(z[0] as i32)&0xFFFFFF,dur/{});\n'.format(100000000//scale)

    str+="\treturn;\n}\n"
    return str

def time_modsqr(n,base,r) :
    N=getN(n,base)
    mask=(1<<base)-1
    bit=n%base
    smask=(1<<bit)-1
    rp=makebig(r,base,N)
    str="pub fn time_modsqr() {\n"
    str+="\tlet mut x: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut z: [SPINT; {}] = [0; {}];\n".format(N,N)

    str+="\t"
    for i in range(0,N) :
        str+="x[{}]={}; ".format(i,hex(rp[i]))
    str+="\n"

    str+="\tnres(&mut x);\n"
    if cyclesorsecs:
        str+="\tunsafe {\n"
        str+="\tlet pre = std::arch::x86_64::_rdtsc();\n"
    str+="\tlet begin = Instant::now();\n"
    str+="\tfor _i in 0..{} {{\n".format(100000//scale)
    str+="\t\tfor _j in 0..500 {\n"
    str+="\t\t\tmodcpy(&x,&mut z);\n"
    str+="\t\t\tmodsqr(&mut z);\n"
    str+="\t\t\tmodcpy(&z,&mut x);\n"
    str+="\t\t\tmodsqr(&mut x);\n"
    str+="\t\t}\n"
    str+="\t}\n"
    if cyclesorsecs:
        str+="\tlet post = std::arch::x86_64::_rdtsc();\n"
    str+="\tlet elapsed = begin.elapsed();\n"
    str+="\tlet dur = (elapsed.as_secs() * 1_000_000_000) + (elapsed.subsec_nanos()) as u64;\n"
    str+="\tredc(&mut z);\n"
    if cyclesorsecs:
        str+='\tprintln!("modsqr check 0x{{:x}} Clock cycles= {{}} Nanosecs= {{}}",(z[0] as i32)&0xFFFFFF,((post-pre)/{}),dur/{});\n'.format(100000000//scale,100000000//scale)
        str+="\t}\n"
    else :
        str+='\tprintln!("modsqr check 0x{{:x}} Nanosecs= {{}}",(z[0] as i32)&0xFFFFFF,dur/{});\n'.format(100000000//scale)
    str+="\treturn;\n}\n"
    return str

def time_modinv(n,base,r) :
    N=getN(n,base)
    mask=(1<<base)-1
    bit=n%base
    smask=(1<<bit)-1
    rp=makebig(r,base,N)
    str="pub fn time_modinv() {\n"
    str+="\tlet mut x: [SPINT; {}] = [0; {}];\n".format(N,N)
    str+="\tlet mut z: [SPINT; {}] = [0; {}];\n".format(N,N)

    str+="\t"
    for i in range(0,N) :
        str+="x[{}]={}; ".format(i,hex(rp[i]))
    str+="\n"

    str+="\tnres(&mut x);\n"
    if cyclesorsecs:
        str+="\tunsafe {\n"
        str+="\tlet pre = std::arch::x86_64::_rdtsc();\n"
    str+="\tlet begin = Instant::now();\n"
    str+="\tfor _i in 0..{} {{\n".format(50000//scale)
    str+="\t\tmodcpy(&x,&mut z);\n"
    str+="\t\tmodinv(None,&mut z);\n"
    str+="\t\tmodcpy(&z,&mut x);\n"
    str+="\t\tmodinv(None,&mut x);\n"
    str+="\t}\n"
    if cyclesorsecs:
        str+="\tlet post = std::arch::x86_64::_rdtsc();\n"
    str+="\tlet elapsed = begin.elapsed();\n"
    str+="\tlet dur = (elapsed.as_secs() * 1_000_000_000) + (elapsed.subsec_nanos()) as u64;\n"
    str+="\tredc(&mut z);\n"
    if cyclesorsecs:
        str+='\tprintln!("modinv check 0x{{:x}} Clock cycles= {{}} Nanosecs= {{}}",(z[0] as i32)&0xFFFFFF,((post-pre)/{}),dur/{});\n'.format(100000//scale,100000//scale)
        str+="\t}\n"
    else :
        str+='\tprintln!("modinv check 0x{{:x}} Microsecs= {{}}",(z[0] as i32)&0xFFFFFF,dur/{});\n'.format(100000//scale)

    str+="\treturn;\n}\n"
    return str

def main() :
    str="fn main() {\n"
    str+='\tprintln!("Rust code - timing some functions - please wait");\n'
    str+="\ttime_modmul();\n"
    str+="\ttime_modsqr();\n"
    str+="\ttime_modinv();\n"
    str+="}\n"
    return str

if len(sys.argv)!=3 :
    print("Syntax error")
    print("Valid syntax - python monty_rust.py <word length> <prime> OR <prime> OR <prime expression>")
    print("For example - python monty_rust.py 64 NIST256")
    print("For example - python monty_rust.py 64 0x01fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffa51868783bf2f966b7fcc0148f709a5d03bb5c9b8899c47aebb6fb71e91386409")
    exit(2)

WL=int(sys.argv[1])
if WL !=16 and WL != 32 and WL !=64 :
    print("Only 16, 32 and 64-bit word lengths supported")
    exit(2)

prime=sys.argv[2]
p=0
base=0

algorithm=False
mp=0

# Note that number base MUST be less than wordlength WL (unsaturated).
# The radix for an n-bit prime is 2^base where base is usually in the range 51-60 for a 64-bit processor, 26-30 for 
# 32-bit and 12-13 for 16 bit
# The radix must be low enough so that the sum of L 2^(2*base) double length integers does not encroach on the sign bit, 
# that is bits(L)+2*base<2*WL, where L=n/base (rounded up) is the number of limbs in the number representation. But 
# at the same time radix should be high enough to minimize L. An underused top digit allows for lazy reduction.

# More named primes can be added here


### Start of user editable area


if WL!=64 :
    inline=False

if prime=="NIST256" :
    p=2**256-2**224+2**192+2**96-1
    #if WL==64 :                          # manual override of default
    #    base=48

if prime=="NIST384" :
    p=2**384-2**128-2**96+2**32-1
#    if WL==64:
#        base=60
#    if WL==32 :
#        base=29

if prime=="X25519" :
    p=2**255-19

if prime=="X448" :
    p=2**448-2**224-1
    if not generic :
        algorithm=True  # if algorithm is known, fix multiple of prime (for modular subtractions) as described in https://eprint.iacr.org/2017/437
        mp=2            # Make sure there is sufficient excess - otherwise change default base. Here assuming Montgomery ladder algorithm. Now no reduction required after modular additions/subtractions.

if prime=="ED448" :
    p=2**448-2**224-1
    if not generic :
        algorithm=True
        mp=4          # Assuming Edwards curve - see https://eprint.iacr.org/2017/437

if prime=="NIST521" :
    p=2**521-1
    #if WL==64 :
    #    base=60
    #if WL==32:
    #    karatsuba=True

if prime=="GM270" :
    p=2**270-2**162-1
    if WL==32 :
        base=29

if prime=="C414" :
    p=2**414-17

if prime=="PM512" :
    p=2**512-569
    if WL==64 :
        base=58

if prime=="GM270" :
    p=2**270-2**162-1
    if WL==32 :
        base=29

if prime=="GM266" :
    p=2**266-2**168-1
    if WL==32 :
        base=28
    if WL==64 :
        base=56

if prime=="GM240" :
    p=2**240-2**183-1
    if WL==32 :
        base=29
    if WL==64 :
        base=61

if prime=="GM360" :
    p=2**360-2**171-1
    if WL==32 :
        base=29
    if WL==64 :
        base=57

if prime=="GM480" :
    p=2**480-2**240-1
    if WL==64 :
        base=60
    if WL==32 :
        base=29

if prime=="GM384" :
    p=2**384-2**186-1
    if WL==32:
        base=29
    if WL==64:
        base=62

if prime=="GM512" :
    p=2**512-2**127-1
    if WL==32:
        base=29
        #karatsuba=True
    if WL==64:
        base=58

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
    if WL==64:
        base=52

if prime=="MFP7" :
    p=2**145*(3**9)*(59**3)*(311**3)*(317**3)*(503**3)-1
    if WL==64:
        base=52

if prime=="MFP1973" :
    p=0x34e29e286b95d98c33a6a86587407437252c9e49355147ffffffffffffffffff
    if WL==64:
        base=52

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
if p==0 :
    if not prime[0].isdigit() :
        print("This named prime not supported")
        exit(2)
    if prime[0]=='0' and prime[1]=='0' :
        field=False
        p=int(prime)
    else :
        p=eval(prime)

n=p.bit_length() 
if n<120 or pow(3,p-1,p)!=1 :
    print("Not a sensible modulus, too small or not a prime")
    exit(2)

PM=False

if base==0 :
    base=getbase(n) #use default

trin=trinomial(p,base)
if trin>0 :
    print("Prime is a lucky trinomial")

b=2**base
N=getN(n,base)  # Number of digits in modulus (to given radix)
xcess=N*base-n  # excess bits in representation - should be 2 or greater, or 0 in which case extra "virtual" limb will be added

if N>9 :
    inline=False

# Is it an exploitable Pseudo-Mersenne?
m=2**n-p
if m>1 and bits(m)+base<WL :
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

#get non-trivial root of unity
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
if karatsuba :
    print("Using Karatsuba for modmul")
else : 
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

if E :
    print("Extra virtual limb added")
    R=getR(n+base,base) #add extra digit
else :
    R=getR(n,base)

# should only be an issue for bad user choice of radix
if xcess<2 and (not E) :
    print("Error - Excess is only one bit, consider change of radix")
    exit(1)


# detect if Montgomery-friendly, or not
b=2**base
ndash=inverse((R-p)%b,b)

# check if Montgomery friendly
fullmonty=True
if ndash==1 :
    fullmonty=False

# represent Montgomery constant to given radix
c=(R*R)%p              # Montgomery constant
cw=makebig(c,base,N)

maxdigit=2**base-1
maxnum=0

# generate fused WL-bit modular multiplication and squaring code for shaped prime
# where the unsaturated radix is 2^base and base<WL


from contextlib import redirect_stdout

print("For correctness compare check value against that generated by C timing code")

# Note that the accumulated partial products must not exceed the double precision limit. 
# If NOT using karatsuba this limit can be 2**(2*WL), otherwise 2**(2*WL-1)
# If NOT using karatsuba use unsigned integer types to store limbs, otherwise use signed types

subprocess.call("rm time.rs", shell=True)

import random
modulus=p
random.seed(42)
ra=random.randint(0,modulus-1)
rb=random.randint(0,modulus-1)
rs=random.randint(0,modulus-1)
ri=random.randint(0,modulus-1)
with open('time.rs', 'w') as f:
    with redirect_stdout(f):
        print("use std::time::Instant;")
        print("pub type SPINT = u{};".format(WL))
        print("pub type SSPINT = i{};".format(WL))
        print("pub type UDPINT = u{};".format(2*WL))
        if karatsuba :
            print("pub type DPINT = i{};".format(2*WL))
        else :
            print("pub type DPINT = u{};".format(2*WL))
        print(prop(n,base))
        print(flat(n,base))
        print(modfsb(n,base))
        #if trin==0 :
        #    print(modint())
        #if trin>0 :
        print(modmli(n,base))
        #else :
        #    print(dummymodmli(n,base))
        print(modmul(n,base))
        print(modsqr(n,base))
        print(modcpy())
        print(modnsqr())
        print(modpro())
        print(modinv())
        print(nres(n,base))
        print(redc(n,base))
        print(time_modmul(n,base,ra,rb))
        print(time_modsqr(n,base,rs))
        print(time_modinv(n,base,ri))
        print(main())

f.close()

subprocess.call("rustc -C opt-level=3 -C target-cpu=native time.rs", shell=True)
print("For timings run ./time")

fname='field.rs'
if not field :
    fname='group.rs'

# output code in final form
with open(fname, 'w') as f:
    with redirect_stdout(f):
        print("\n//Automatically generated modular arithmetic Rust code")
        print("//Command line : python {} {} {}".format(sys.argv[0], sys.argv[1], sys.argv[2]))
        print("//Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)\n")
        if makepublic :
            print("pub type SPINT = u{};".format(WL))
            print("pub type SSPINT = i{};".format(WL))
            print("pub type UDPINT = u{};".format(2*WL))
            if karatsuba :
                print("pub type DPINT = i{};".format(2*WL))
            else :
                print("pub type DPINT = u{};".format(2*WL))
        else :
            print("type SPINT = u{};".format(WL))
            print("type SSPINT = i{};".format(WL))
            print("type UDPINT = u{};".format(2*WL))
            if karatsuba :
                print("type DPINT = i{};".format(2*WL))
            else :
                print("type DPINT = u{};".format(2*WL))

        print(prop(n,base))
        print(flat(n,base))
        print(modfsb(n,base))
        print(modadd(n,base))
        print(modsub(n,base))
        print(modneg(n,base))
        #if trin>0 :
        print(modmli(n,base))
        #else :
        #    print(dummymodmli(n,base))
        print(modmul(n,base))
        print(modsqr(n,base))
        print(modcpy())
        print(modnsqr())
        print(modpro())
        print(modinv())
        print(nres(n,base))
        print(redc(n,base))
        print(modis1(n,base))
        print(modis0(n,base))
        print(modzer())
        print(modone())
        print(modint())
        print(modqr())
        print(modcmv())
        print(modcsw())
        print(modsqrt())
        print(modshl(n,base))
        print(modshr(n,base))
        print(mod2r())
        print(modexp())
        print(modimp())
        print(modsign())
        print(modcmp())
        if makepublic :
            print("pub const NLIMBS: usize = {};".format(N))
            print("pub const RADIX: usize = {};".format(base))
            print("pub const NBITS: usize = {};".format(n))
            print("pub const NBYTES: usize = {};\n".format(Nbytes))
            print("pub const MERSENNE: bool = false;")
            print("pub const MONTGOMERY: bool = true;\n")
            if trin>0 :
                print("pub const MULBYINT: bool = true;\n")
            else :
                print("pub const MULBYINT: bool = false;\n")
        else :
            print("const NLIMBS: usize = {};".format(N))
            print("const RADIX: usize = {};".format(base))
            print("const NBITS: usize = {};".format(n))
            print("const NBYTES: usize = {};\n".format(Nbytes))
            print("const MERSENNE: bool = false;")
            print("const MONTGOMERY: bool = true;\n")
            if trin>0 :
                print("const MULBYINT: bool = true;\n")
            else :
               print("const MULBYINT: bool = false;\n")
f.close()

if formatted :
    if field :
        subprocess.call("rustfmt field.rs", shell=True)
    else :
        subprocess.call("rustfmt group.rs", shell=True)
if field :
    print("field code is in field.rs")
else :
    print("group code is in group.rs")

sys.exit(base)
