# Python program to generate reasonably efficient C/C++ modular arithmetic code for pseudo-mersenne primes on a 16, 32 or 64-bit processor
# Modulus should be a pseudo-mersenne of the form 2^n-m
# uses unsaturated radix
#
# In particular this script generates code for primes like
#
# X25519 = 2^255-19
# PM266 = 2^266-3
# NIST521 = 2^521-1
#
# requires addchain utility in the path - see https://github.com/mmcloughlin/addchain 
#
# Execute this program as: python pseudo_rust.py 64 X25519
# Code is output to file field.rs
#
# Mike Scott 22nd April 2024
# TII
#

# Some default settings

cyclesorsecs=True # use a cpu cycle counter, otherwise just provide timings
karatsuba=False  # default setting
formatted=True # pretty up the final output
inline=True # consider encouraging inlining
generic=True # set to False if algorithm is known in advance, in which case modadd and modsub can be faster - see https://eprint.iacr.org/2017/437
scale=1 # set to 10 or 100 for faster timing loops. Default to 1
makepublic=False # Make functions public

import sys
import subprocess

# Determine optimal radix
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
        while limbs*base<n :
            base=base+1
        if base>WL-3 :
            continue
        if limbs*(2**base-1)**2 < 2**limit :
            break
    return base


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
def caddp(x) :
    str=""
    str+="\tn[0]-=({})&(carry as SPINT);\n".format(m*x)
    str+="\tn[{}]+=(0x{:x} as SPINT)&(carry as SPINT);\n ".format(N-1,x*TW)
    return str

#add p
def addp(x) :
    str=""
    str+="\tn[0]-={} as SPINT;\n".format(m*x)
    str+="\tn[{}]+=0x{:x} as SPINT;\n".format(N-1,x*TW)
    return str

#subtract p
def subp(x) :
    str=""
    str+="\tn[0]+={} as SPINT;\n".format(m*x)
    str+="\tn[{}]-=0x{:x} as SPINT;\n".format(N-1,x*TW)
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
    str+="fn modfsb(n: &mut[SPINT]) -> bool{\n"
    #str+="\tlet q=(1 as SPINT)<<{};\n".format(base)
    str+=subp(1)
    str+="\treturn flatten(n);\n}\n"
    return str

#modular addition
def modadd(n,m,base) :
    str="//Modular addition - reduce less than 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modadd(b: &[SPINT],n: &mut [SPINT]) {\n"
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
def modsub(n,m,base) :
    str="//Modular subtraction - reduce less than 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modsub(b: &[SPINT],n: &mut [SPINT]) {\n"
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
def modneg(n,m,base) :
    str="//Modular negation\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modneg(n: &mut [SPINT]) {\n"
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

#multiplication macro
def getZM(str,row,n,m,base) :
    N=getN(n,base)
    xcess=N*base-n
    mm=m*2**xcess
    k=row+1
    L=N-1
    if karatsuba :
        i=row+N
        if row<N-1 :
            str+="\ttt=d{}-d{}; ".format(N-1,row)
            for m in range(N-1,int(i/2),-1) :
                str+="tt+=(((c[{}] as SSPINT)-(c[{}] as SSPINT)) as DPINT) * (((b[{}] as SSPINT)-(b[{}] as SSPINT)) as DPINT); ".format(m,i-m, i-m, m)
            if overflow :
                str+=" lo=(tt as SPINT)&mask;"
                if row==0 :
                    str+=" t+=d{}+(lo as DPINT)*0x{:x};".format(row,mm)
                else :
                    if bad_overflow :
                        str+=" t+=d{}+((lo as DPINT)+hi)*0x{:x};".format(row,mm)
                    else :
                        str+=" t+=d{}+((lo+hi) as DPINT)*0x{:x};".format(row,mm)
                if bad_overflow :
                    str+=" hi=tt>>{};".format(base)
                else :
                    str+=" hi=(tt>>{}) as SPINT;".format(base)
            else :
                str+="tt*=0x{:x};".format(mm)
                str+=" t+=d{}+tt;".format(row)
        else :
            str+="\tt+=d{};".format(N-1)

        i=row
        for m in range(i,int(i/2),-1) :
           str+=" t+=(((c[{}] as SSPINT)-(c[{}] as SSPINT)) as DPINT) * (((b[{}] as SSPINT)-(b[{}] as SSPINT)) as DPINT); ".format(m,i - m, i - m, m) 
        if row==N-1 and overflow :
            str+=" t+=(hi as DPINT)*0x{:x};".format(mm)
        str+=" let v{}=(t as SPINT) & mask; t=t>>{};\n".format(row,base)
        return str

    first=True
    while k<N :
        if EPM :
            if first :
                str+="\t"
            else :
                str+=" "
            str+="t+=(mc{} as DPINT)*(b[{}] as DPINT);".format(k,L)
        else :
            if first :
                str+="\ttt=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,L)
                first=False
            else :
                str+=" tt+=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,L)
        L-=1
        k+=1
    if row<N-1:
        if overflow :
            str+=" lo=(tt as SPINT)&mask;"
            if row==0 :
                str+=" t+=(lo as DPINT)*0x{:x};".format(mm)
            else :
                if bad_overflow :
                    str+=" t+=((lo as DPINT)+hi)*0x{:x};".format(mm)
                else :
                    str+=" t+=((lo+hi) as DPINT)*0x{:x};".format(mm)
            if bad_overflow :
                str+=" hi=tt>>{};".format(base)
            else :
                str+=" hi=(tt>>{}) as SPINT;".format(base)   
        else :
            if not EPM :
                str+=" tt*=0x{:x};".format(mm)
                str+=" t+=tt;"
    else :
        first=True
    
    k=0
    while k<=row :
        if first :
            str+="\tt+=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,row-k)
            first=False
        else :
            str+=" t+=(c[{}] as DPINT)*(b[{}] as DPINT);".format(k,row-k)
        k+=1
    if row==N-1 and overflow :
        str+=" t+=(hi as DPINT)*0x{:x};".format(mm)
    str+=" let v{}=(t as SPINT) & mask; t=t>>{};\n".format(row,base)
    return str

#squaring macro
def getZS(str,row,n,m,base) :
    N=getN(n,base)
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
                str+="t+=(mc{} as UDPINT)*(tc{} as UDPINT);".format(k,L)
            else :
                if first :
                    str+="\tt+=(mc{} as UDPINT)*(c[{}] as UDPINT);".format(k,L)
                    first=False
                else :
                    str+=" t+=(mc{} as UDPINT)*(c[{}] as UDPINT);".format(k,L)
        else :
            if first :
                str+="\ttt=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,L)
                first=False
            else :
                str+=" tt+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,L)
        L-=1
        k+=1
    if dble :
        if not EPM :
            str+=" tt*=2;"
    if k==L :
        if EPM :
            if first :
                str+="\t"
            else :
                str+=" "  
            str+="t+=(mc{} as UDPINT)*(c[{}] as UDPINT);".format(k,k)
        else :
            if first :
                str+="\ttt=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,k)
                first=False
            else :
                str+=" tt+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,k)
    first=True
    if row<N-1:
        if overflow :
            str+=" lo=(tt as SPINT)&mask;"
        else :
            if not EPM :
                str+=" tt*=0x{:x};".format(mm)
                str+=" t+=tt;"
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
            str+="t+=(c[{}] as UDPINT)*(tc{} as UDPINT);".format(k,L)
        else :
            if first :
                str+="t2=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,L)
                first=False
            else :
                str+=" t2+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,L)
        k+=1
        L-=1

    if dble :
        if not EPM :
            str+=" t2*=2;"
    if k==L :
        if EPM :
            str+=" t+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,k)
        else :
            if first :
                str+="t2=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,k)
                first=False
            else :
                str+=" t2+=(c[{}] as UDPINT)*(c[{}] as UDPINT);".format(k,k)
 
    if overflow :
        if row==N-1 : 
            str+=" t+=(hi as UDPINT)*0x{:x};".format(mm)
        else :
            if row==0 :
                str+=" t2+=(lo as UDPINT)*0x{:x};".format(mm)
            else :
                if bad_overflow :
                    str+=" t2+=((lo as UDPINT)+hi)*0x{:x};".format(mm)
                else :
                    str+=" t2+=((lo+hi) as UDPINT)*0x{:x};".format(mm)
            if bad_overflow :
                str+=" hi=tt>>{};".format(base)
            else :
                str+=" hi=(tt>>{}) as SPINT;".format(base)

    if not EPM :
        str+=" t+=t2;"
    str+=" let v{}=(t as SPINT) & mask; t=t>>{};\n".format(row,base)
    return str

# second reduction pass
def second_pass(str,n,m,base) :
    N=getN(n,base)
    mask=(1<<base)-1
    xcess=N*base-n
    smask=(1<<(base-xcess))-1
    str+="// second reduction pass\n"    

    if xcess>0 :
        str+="\tlet mut nv=v{};\n".format(N-1)
    else :
        str+="\tlet nv=v{};\n".format(N-1)
    if fred :
        str+="\tlet mut ut=t as SPINT;\n"
        if xcess>0 :
            str+= "\tut=(ut<<{})+(nv>>{}); nv&=0x{:x};\n".format(xcess,base-xcess,smask)
    else :
        str+="\tlet mut ut=t as UDPINT;\n"
        if xcess>0 :
            str+= "\tut=(ut<<{})+((nv>>{}) as UDPINT); nv&=0x{:x};\n".format(xcess,base-xcess,smask)

    k=0
    if m>1 :
        str+= "\tut*={};\n".format(m)

    if carry_on :
        str+= "\tlet mut s=v0+((ut as SPINT)&mask);\n"
        str+= "\tc[0]=(s&mask) as SPINT;\n"
        str+="\tut=((s>>{}) as UDPINT)+(ut>>{});\n".format(base,base)
        str+="\ts=v1+((ut as SPINT) & mask);\n"
        str+= "\tc[1]=(s&mask) as SPINT;\n"
        k+=1
    else :
        str+= "\tlet s=v0+((ut as SPINT)&mask);\n"
        str+= "\tc[0]=(s&mask) as SPINT;\n"

    str+= "\tlet carry=(s>>{})+((ut>>{}) as SPINT);\n".format(base,base)
    k=k+1

    str+= "\tc[{}]=v{}+carry;\n".format(k,k)

    for i in range(k+1,N-1) :
        str+= "\tc[{}]=v{};\n".format(i,i)
    str+="\tc[{}]=nv;\n".format(N-1)

    return str

# modular multiplication
# Note that allowing inlining gives significant speed-up
def modmul(n,m,base) :
    N=getN(n,base)
    mask=(1<<base)-1
    xcess=N*base-n
    str="// Modular multiplication, c=c*b mod 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modmul(b: &[SPINT],c: &mut [SPINT]) {\n"
    str+="\tlet mut t=0 as DPINT;\n"

    if karatsuba :
        str+="\tlet mut tt:DPINT;\n"
        str+="\tlet d0=(c[0] as DPINT)*(b[0] as DPINT);\n"
        for i in range(1,N) :
            str+="\tlet d{}=d{}+(c[{}] as DPINT)*(b[{}] as DPINT);\n".format(i, i-1, i, i)
    else :
        if EPM :
            for i in range(1,N) :
                str+="\tlet mc{}=c[{}]*0x{:x};\n".format(i,i,mm)
        else :
            str+="\tlet mut tt:DPINT;\n"
    if overflow :
        str+="\tlet mut lo:SPINT;\n"
        if bad_overflow :
            str+="\tlet mut hi:DPINT;\n"
        else :
            str+="\tlet mut hi:SPINT;\n"
    str+="\tlet mask=((1 as SPINT)<<{})-1;\n".format(base)

    for row in range(0,N) :
        str=getZM(str,row,n,m,base)

    str=second_pass(str,n,m,base)

    str+="\treturn;\n}\n"
    return str

# modular squaring
# Note that allowing inlining gives significant speed-up
def modsqr(n,m,base) :
    N=getN(n,base)
    mask=(1<<base)-1
    xcess=N*base-n
    str="// Modular squaring, c=c*c mod 2p\n"
    str+="#[allow(unused_variables)]\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modsqr(c: &mut [SPINT]) {\n"
    str+="\tlet mut t=0 as UDPINT;\n"
    if EPM :
       for i in range(1,N) :
           str+="\tlet tc{}=c[{}]*2;\n".format(i,i)
       for i in range(1,N) :
           str+="\tlet mc{}=c[{}]*0x{:x};\n".format(i,i,mm)     
    else :
        str+="\tlet mut tt:UDPINT;\n"
        str+="\tlet mut t2:UDPINT;\n"
    str+="\tlet mask=((1 as SPINT)<<{})-1;\n".format(base)

    if overflow :
        str+="\tlet mut lo:SPINT;\n"
        if bad_overflow :
            str+="\tlet mut hi:UDPINT;\n"
        else :
            str+="\tlet mut hi:SPINT;\n"

    for row in range(0,N) :
        str=getZS(str,row,n,m,base)
  
    str=second_pass(str,n,m,base)

    str+="\treturn;\n}\n"
    return str

def modmli(n,m,base) :
    N=getN(n,base)
    mask=(1<<base)-1
    xcess=N*base-n
    str="// Modular multiplication by an integer, c=c*b mod 2p\n"
    if inline :
        str+="#[inline]\n"
    if makepublic :
        str+="pub "
    str+="fn modmli(b: usize,c: &mut [SPINT]) {\n"
    str+="\tlet mut t=0 as UDPINT;\n"
    str+="\tlet mask=((1 as SPINT)<<{})-1;\n".format(base)

    for i in range(0,N) :
        str+="\tt+=(c[{}] as UDPINT)*(b as UDPINT); ".format(i)
        str+="let v{}=(t as SPINT) & mask; t=t>>{};\n".format(i,base)

    str=second_pass(str,n,m,base)

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
    str="//Calculate progenitor - use optimal addition chain\n"
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
            #str+="\tfor _i in 0..{} {{\n".format(int(info[3]))
            #str+="\t\tmodsqr(&mut {});\n".format(info[1])
            #str+="\t}\n"
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
    #str+="\tfor _i in 0..={} {{\n".format(PM1D2)
    #str+="\t\tmodsqr(&mut t);\n"
    #str+="\t}\n"
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
        #str+="\tfor _i in 0..{}-1 {{\n".format(PM1D2)
        #str+="\t\tmodsqr(&mut r);\n\t}\n"
    str+="\treturn modis1(&r) || modis0(x);\n}\n"
    return str

#modular square root
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
        #str+="\t\tfor _ in 1..k-1 { \n"
        #str+="\t\t\tmodsqr(&mut b);\n\t\t}\n"

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

# convert to internal form
def nres(n,base) :
    str="//Convert n to n-residue form, n=nres(m) \n"
    if makepublic :
        str+="pub "
    str+="fn nres(_n: &mut [SPINT]) {\n"
    str+="\treturn;\n}\n"
    return str 

#convert back to integer form
def redc(n,base) :
    str="//Convert m back to normal form, m=redc(n) \n"
    if makepublic :
        str+="pub "
    str+="fn redc(m: &mut [SPINT]) {\n"
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
    if WL==16 :
        str+="let r=0xa55a;\n"
    if WL==32 :
        str+="let r=0x5aa5a55a;\n"
    if WL==64 :
        str+="let r=0x3cc3c33c5aa5a55a;\n"
    str+="\t\tlet c0=(!bb)&(r+1);\n"
    str+="\t\tlet c1=bb|r;\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\tlet s=g[i];\n"
    str+="\t\tlet t=f[i];\n"
    str+="\t\tlet w=r*(t+s);\n"
    str+="\t\tunsafe{core::ptr::write_volatile(&mut f[i],c0*t+c1*s);}\n"
    str+="\t\tf[i]-=w;\n"
    str+="\t\tunsafe{core::ptr::write_volatile(&mut g[i],c0*s+c1*t);}\n"
    str+="\t\tg[i]-=w;\n\t}\n"
    str+="\treturn;\n}\n"
    return str

#conditional move
def modcmv() :
    str="//conditional move g to f if d=1\n"
    str+="#[inline(never)]\n"
    if makepublic :
        str+="pub "
    str+="fn modcmv(b: usize,g: &[SPINT],f: &mut [SPINT]) {\n"
    if WL==16 :
        str+="let r=0xa55a;\n"
    if WL==32 :
        str+="let r=0x5aa5a55a;\n"
    if WL==64 :
        str+="let r=0x3cc3c33c5aa5a55a;\n"
    str+="\tlet bb = b as SPINT;\n"
    str+="\t\tlet c0=(!bb)&(r+1);\n"
    str+="\t\tlet c1=bb|r;\n"
    str+="\tfor i in 0..{} {{\n".format(N)
    str+="\t\tlet s=g[i];\n"
    str+="\t\tlet t=f[i];\n"
    str+="\t\tunsafe{core::ptr::write_volatile(&mut f[i],c0*t+c1*s);}\n"
    str+="\t\tf[i]-=r*(t+s);\n\t}\n"
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
    str+="\ta[n] <<= m;\n}\n"
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
    print("Valid syntax - python pseudo_rust.py <word length> <prime> OR <prime> OR <prime expression>")
    print("For example - python pseudo_rust.py 64 X25519")
    print("For example - python pseudo_rust.py 64 2**255-19")
    exit(2)

WL=int(sys.argv[1])
if WL !=16 and WL != 32 and WL !=64 :
    print("Only 16, 32 and 64-bit word lengths supported")
    exit(2)

prime=sys.argv[2]

n=0
base=0
m=0

algorithm=False
mp=0
# Note that number base MUST be less than wordlength WL (unsaturated).
# The radix for an n-bit prime is 2^base where base is usually in the range 51-60 for a 64-bit processor, 26-30 for 
# 32-bit and 12-13 for 16-bit.
# The radix must be low enough so that the sum of L 2^(2*base) double length integers does not encroach on the sign bit, 
# that is bits(L)+2*base<2*WL, where L=n/base (rounded up) is the number of limbs in the number representation. But 
# at the same time radix should be high enough to minimize L. An underused top digit allows for lazy reduction.


### Start of user editable area

if WL!=64 :
    inline=False


# More named primes can be added here
p=0
if prime=="PM266" :
    p=2**266-3
    #if WL==64:      # manual override possible
    #    base=54
    #if WL==32 :
    #    base=29

if prime=="NUMS256W" :
    p=2**256-189

if prime=="NUMS256E" :
    p=2**256-189

if prime=="NIST521" :
    p=2**521-1
    #if WL==32 :
    #    karatsuba=True

if prime=="X25519" :
    p=2**255-19
    if not generic :
        algorithm=True  # if algorithm is known, fix multiple of prime mp (for modular subtractions) as described in https://eprint.iacr.org/2017/437
        mp=2            # Make sure there is sufficient excess - otherwise change default radix. Here assuming Montgomery ladder algorithm. Now no reduction required after modular additions/subtractions.

if prime=="ED25519" :
    p=2**255-19
    if not generic :
        algorithm=True  # if algorithm is known, fix multiple of prime mp (for modular subtractions) as described in https://eprint.iacr.org/2017/437
        mp=4            # Make sure there is sufficient excess - otherwise change default radix. Here assuming Montgomery ladder algorithm. Now no reduction required after modular additions/subtractions.

if prime=="C2065" :
    p=2**206-5

if prime=="PM336" :
    p=2**336-3
    #if WL==32 :
    #    karatsuba=True

if prime=="PM383" :
    p=2**383-187

if prime=="C414" :
    p=2**414-17

if prime=="PM512" :
    p=2**512-569
    if WL==16 :
        base=12

if prime=="SECP256K1" :
    p=2**256-2**32-977

if p==0 :
    if prime[0].isdigit() :
        p=eval(prime)
    else :
        print("This named prime not supported (not a pseudo-Mersenne?)")
        exit(2)

### End of user editable area

n=p.bit_length() 
if n<120 or pow(3,p-1,p)!=1 :
    print("Not a sensible modulus, too small or not a prime")
    exit(2)
#if n>360 and WL==16 :
#	print("Modulus probably too big for 16-bit processor")
#	exit(1)

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
PE=(p-1-e)//(2*e)

if m>=b : 
    print("Not an exploitable pseudo-Mersenne - Use Montgomery method instead")
    exit(1)

# get number of limbs
N=getN(n,base)
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

# get number of bytes for export
Nbytes=n//8
if (n%8)!=0 :
    Nbytes+=1

# get non-trivial root of unity
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

overflow=False
bad_overflow=False
if (b-1)*(b-1)*mm*N >= 2**(2*WL) :
    overflow=True
    print("Possibility of overflow... using alternate method")
    if karatsuba :
        if (N-1)*(b-1)**2 >= 2**(2*WL-4) :
            bad_overflow=True
    else :
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
    print("Must propagate carries one limb further in second pass")
    carry_on=True

modulus=p

# generate WL bit modular multiplication and squaring code for Mersenne prime 2^n-m
# where the unsaturated radix is 2^base and base<WL 

from contextlib import redirect_stdout

print("For correctness compare check value against that generated by C timing code")

# Note that the accumulated partial products must not exceed the double precision limit. 
# If NOT using karatsuba this limit can be 2**(2*WL), otherwise 2**(2*WL-1)
# If NOT using karatsuba use unsigned integer types to store limbs, otherwise use signed types

# output code in a form suitable for timing

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
        print(modmul(n,m,base))
        print(modsqr(n,m,base))
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

# output code in final form
with open('field.rs', 'w') as f:
    with redirect_stdout(f):
        print("\n//Automatically generated modular arithmetic Rust code for pseudo-Mersenne primes")
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
        print(modmli(n,m,base))
        print(modadd(n,m,base))
        print(modsub(n,m,base))
        print(modneg(n,m,base))
        print(modmul(n,m,base))
        print(modsqr(n,m,base))
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
        print(modexp())
        print(modimp())
        print(modsign())
        print(modcmp())
        if makepublic :
            print("pub const NLIMBS: usize = {};".format(N))
            print("pub const RADIX: usize = {};".format(base))
            print("pub const NBITS: usize = {};".format(n))
            print("pub const NBYTES: usize = {};\n".format(Nbytes))
            print("pub const MERSENNE: bool = true;")
            print("pub const MONTGOMERY: bool = false;\n")
            print("pub const MULBYINT: bool = true;\n")
        else :
            print("const NLIMBS: usize = {};".format(N))
            print("const RADIX: usize = {};".format(base))
            print("const NBITS: usize = {};".format(n))
            print("const NBYTES: usize = {};\n".format(Nbytes))
            print("const MERSENNE: bool = true;")
            print("const MONTGOMERY: bool = false;\n")
            print("const MULBYINT: bool = true;\n")
f.close()

if formatted :
    subprocess.call("rustfmt field.rs", shell=True)
print("Field code is in field.rs")

sys.exit(0)
