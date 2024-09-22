# Script for automatic generation of code for elliptic curve cryptography
# rust version
# See below for curve definitions for some popular curves. Insert your own!
# Mike Scott 3rd September 2024
# TII
#
import sys
import subprocess

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

def replacefromfile(namefile,oldtext,newfile):
    f = open(namefile,'r')
    filedata = f.read()
    f.close()

    nf = open(newfile,'r')
    newdata = nf.read()
    nf.close

    newdata = filedata.replace(oldtext,newdata)

    f = open(namefile,'w')
    f.write(newdata)
    f.close()


if len(sys.argv)!=3 :
    print("Syntax error")
    print("Valid syntax - python curve_rust.py <word length> <curve>")
    print("For example - python curve_rust.py 64 ED25519")
    exit(0);

WL=int(sys.argv[1])
if WL !=16 and WL != 32 and WL !=64 :
    print("Only 16, 32 and 64-bit word lengths supported")
    exit(0)

curve=sys.argv[2]

inline=True # consider encouraging inlining
if WL!=64 :
    inline=False

PSEUDO=1
MONTY=2

WEIERSTRASS=1
EDWARDS=2
MONTGOMERY=3

# More curves can be added here
p=0
q=0
cof=0
prime_type=0
curve_type=0
B=0
X=0
Y=0
A=0

if curve=="ED25519" :
    p=2**255-19
    q=0x1000000000000000000000000000000014DEF9DEA2F79CD65812631A5CF5D3ED
    cof=3
    A=-1
    prime_type=PSEUDO
    curve_type=EDWARDS
    B=0x52036CEE2B6FFE738CC740797779E89800700A4D4141D8AB75EB4DCA135978A3
    X=0x216936D3CD6E53FEC0A4E231FDD6DC5C692CC7609525A7B2C9562D608F25D51A
    Y=0x6666666666666666666666666666666666666666666666666666666666666658

if curve=="ED448" :
    p=2**448-2**224-1
    q=(p + 1 - 28312320572429821613362531907042076847709625476988141958474579766324) // 4;
    cof=2
    A=1
    prime_type=MONTY
    curve_type=EDWARDS
    B=-39081
    X=0x4f1970c66bed0ded221d15a622bf36da9e146570470f1767ea6de324a3d3a46412ae1af72ab66511433b80e18b00938e2626a82bc70cc05e
    Y=0x693f46716eb6bc248876203756c9c7624bea73736ca3984087789c1e05a0c2d73ad3ff1ce67c39c4fdbd132c4ed7c8ad9808795bf230fa14

if curve=="NUMS256E" :
    p=2**256-189
    q=0x4000000000000000000000000000000041955AA52F59439B1A47B190EEDD4AF5
    cof=2
    A=1
    prime_type=PSEUDO
    curve_type=EDWARDS
    B=-15342
    X=34

if curve=="NUMS256W" :
    p=2**256-189
    q=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE43C8275EA265C6020AB20294751A825
    cof=1
    A=-3
    prime_type=PSEUDO
    curve_type=WEIERSTRASS
    B=152961
    X=2

if curve=="NIST256" :
    p=115792089210356248762697446949407573530086143415290314195533631308867097853951
    q=115792089210356248762697446949407573529996955224135760342422259061068512044369
    cof=1
    A=-3
    prime_type=MONTY
    curve_type=WEIERSTRASS
    B=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
    X=0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
    Y=0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5

if curve=="SECP256K1" :
    p=2**256 - 2**32 - 977
    q=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
    cof=1
    A=0
    if WL==16 or WL==32 :
        prime_type=MONTY
    else :
        prime_type=PSEUDO
    curve_type=WEIERSTRASS
    B=7
    X=0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
    Y=0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8

if p==0:
    print("This curve not supported")

radix=0
if prime_type==PSEUDO :
    radix=subprocess.run("python pseudo_rust.py "+str(WL)+" "+curve, shell=True).returncode
if prime_type==MONTY :
    radix=subprocess.run("python monty_rust.py "+str(WL)+" "+curve, shell=True).returncode

if radix==0 :
    print("Bad curve");
    exit(0)

base=2**radix

mulbyint=True
if prime_type==MONTY and not trinomial(p,radix) :
    mulbyint=False

bts=p.bit_length()

limbs=int(bts/radix)
if bts%radix != 0 :
    limbs+=1

lim=2**28
small_b=False
small_x=False

if abs(B)<lim :
    small_b=True
if abs(X)<lim : 
    small_x=True

R=1
if bts%radix==0 :
    R=2**((limbs+1)*radix)
else :
    R=2**(limbs*radix)

if prime_type==MONTY :
    if not small_b :
        B=(B*R)%p
    if not small_x :
        X=(X*R)%p
        Y=(Y*R)%p

B3=3*B

strng="\n"

strng+="const COF:usize = {};\n".format(cof)
strng+="const CONSTANT_A: isize = {};\n".format(A)

if not small_b or not mulbyint:
    B3=B3%p
    strng+="const CONSTANT_B: isize = 0;\n"
    strng+="const constant_b: [SPINT;{}]=[".format(limbs)
    for i in range(0,limbs-1) :
        strng+="{},".format(hex(B%base))
        B//=base
    strng+="{}".format(hex(B))
    strng+=("];\n")
    if curve_type==WEIERSTRASS :
        strng+="const constant_b3: [SPINT;{}]=[".format(limbs)
        for i in range(0,limbs-1) :
            strng+="{},".format(hex(B3%base))
            B3//=base
        strng+="{}".format(hex(B3))
        strng+=("];\n")
else :
    strng+="const CONSTANT_B: isize = {};\n".format(B)

if not small_x or not mulbyint:
    strng+="const CONSTANT_X: usize = 0;\n"
    strng+="const constant_x: [SPINT;{}]=[".format(limbs)
    for i in range(0,limbs-1) :
        strng+="{},".format(hex(X%base))
        X//=base
    strng+="{}".format(hex(X))
    strng+="];\n"
else :
    strng+="const CONSTANT_X: usize = {};\n".format(X)

if not small_x or not mulbyint :
    strng+="const constant_y: [SPINT;{}]=[".format(limbs)
    for i in range(0,limbs-1) :
        strng+="{},".format(hex(Y%base))
        Y//=base
    strng+="{}".format(hex(Y))
    strng+="];\n"

from contextlib import redirect_stdout

with open('curve.rs', 'w') as f:
    with redirect_stdout(f):
        print(strng)
f.close()

with open('point.rs', 'w') as f:
    with redirect_stdout(f):
        print("// elliptic curve point in projective coordinates")
        print("const WORDLENGTH: usize = {};".format(WL))
        print("#[derive(Clone)]")
        print("pub struct ECP {")
        print("\tx: [u{};{}],".format(WL,limbs))
        print("\ty: [u{};{}],".format(WL,limbs))
        print("\tz: [u{};{}],".format(WL,limbs))
        print("}\n")
        print("#[allow(non_snake_case)]")
        print("impl ECP {")
        print("\tpub fn new() -> ECP {")
        print("\t\tECP {")
        print("\t\t\tx: [0;{}],".format(limbs))
        print("\t\t\ty: [0;{}],".format(limbs))
        print("\t\t\tz: [0;{}],".format(limbs))
        print("\t\t}")
        print("\t}")
        print("}")
        f.close()

subprocess.run("python monty_rust.py "+str(WL)+" "+curve.lower(), shell=True)

print("field code is in field.rs")
print("curve definition is in curve.rs")
print("point definition is in point.rs")

if curve_type==WEIERSTRASS :
    replacefromfile("weierstrass.rs","@field@","field.rs")
    replacefromfile("weierstrass.rs","@curve@","curve.rs");
    replacefromfile("weierstrass.rs","@point@","point.rs")
if curve_type==EDWARDS:
    replacefromfile("edwards.rs","@field@","field.rs")
    replacefromfile("edwards.rs","@curve@","curve.rs");
    replacefromfile("edwards.rs","@point@","point.rs")