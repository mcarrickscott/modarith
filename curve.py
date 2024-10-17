# Script for automatic generation of code for elliptic curve cryptography
# C version
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
    print("Valid syntax - python curve.py <word length> <curve>")
    print("For example - python curve.py 64 ED25519")
    exit(0);

WL=int(sys.argv[1])
if WL !=16 and WL != 32 and WL !=64 :
    print("Only 16, 32 and 64-bit word lengths supported")
    exit(0)

curve=sys.argv[2]

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

if curve=="NIST384" :
    p=39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319
    q=p + 1 - 1388124618062372383606759648309780106643088307173319169677
    cof=1
    A=-3
    prime_type=MONTY
    curve_type=WEIERSTRASS
    B=27580193559959705877849011840389048093056905856361568521428707301988689241309860865136260764883745107765439761230575
    X=0xaa87ca22be8b05378eb1c71ef320ad746e1d3b628ba79b9859f741e082542a385502f25dbf55296c3a545e3872760ab7
    Y=0x3617de4a96262c6f5d9e98bf9292dc29f8f41dbd289a147ce9da3113b5f0b8c00a60b1ce1d7e819d7a431d7c90ea0e5f

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
    radix=subprocess.run("python pseudo.py "+str(WL)+" "+curve, shell=True).returncode
if prime_type==MONTY :
    radix=subprocess.run("python monty.py "+str(WL)+" "+curve, shell=True).returncode

if radix==0 :
    print("Bad curve");
    exit(0)



base=2**radix

mulbyint=True
if prime_type==MONTY and trinomial(p,radix)==0 :
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

strng+="#define COF {}\n".format(cof)
strng+="#define CONSTANT_A {}\n".format(A)

if not small_b or not mulbyint:
    B3=B3%p
    strng+="const spint constant_b[{}]={{".format(limbs)
    for i in range(0,limbs-1) :
        strng+="{},".format(hex(B%base))
        B//=base
    strng+="{}".format(hex(B))
    strng+=("};\n")
    if curve_type==WEIERSTRASS :
        strng+="const spint constant_b3[{}]={{".format(limbs)
        for i in range(0,limbs-1) :
            strng+="{},".format(hex(B3%base))
            B3//=base
        strng+="{}".format(hex(B3))
        strng+=("};\n")
else :
    strng+="#define CONSTANT_B {}\n".format(B)

if not small_x or not mulbyint:
    strng+="const spint constant_x[{}]={{".format(limbs)
    for i in range(0,limbs-1) :
        strng+="{},".format(hex(X%base))
        X//=base
    strng+="{}".format(hex(X))
    strng+="};\n"
else :
    strng+="#define CONSTANT_X {}\n".format(X)

if not small_x or not mulbyint :
    strng+="const spint constant_y[{}]={{".format(limbs)
    for i in range(0,limbs-1) :
        strng+="{},".format(hex(Y%base))
        Y//=base
    strng+="{}".format(hex(Y))
    strng+="};\n"

from contextlib import redirect_stdout

with open('curve.c', 'w') as f:
    with redirect_stdout(f):
        print(strng)
f.close()

with open('point.h', 'w') as f:
    with redirect_stdout(f):
        print("// elliptic curve point in projective coordinates")
        print("#include <stdint.h>\n")
        print("#ifndef",curve)
        print("#define",curve)
        print("#endif")
        print("#define WORDLENGTH {}".format(WL))
        print("struct xyz {")
        print("\tuint{}_t x[{}];".format(WL,limbs))
        print("\tuint{}_t y[{}];".format(WL,limbs))
        print("\tuint{}_t z[{}];".format(WL,limbs))
        print("};")
        print("typedef struct xyz point;")
        f.close()

subprocess.run("python monty.py "+str(WL)+" "+curve.lower(), shell=True)

print("field code is in field.c")
print("curve definition is in curve.c")
print("point definition is in point.h")

if curve_type==WEIERSTRASS :
    replacefromfile("weierstrass.c","@field@","field.c")
    replacefromfile("weierstrass.c","@curve@","curve.c");
    replacefromfile("curve.h","@point@","point.h")
if curve_type==EDWARDS:
    replacefromfile("edwards.c","@field@","field.c")
    replacefromfile("edwards.c","@curve@","curve.c");
    replacefromfile("curve.h","@point@","point.h")
