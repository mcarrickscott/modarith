# generate corner case testvectors for finite field arithmetic, given elliptic curve modulus
# Test Vectors generated in Python, then tested against C code generated from scripts
# Completely automatic
#
# Download and install modarith code from github
# For example in modarith directory, to test the 32-bit code, run
# python3 edge.py 32 X25519
# or
# python3 edge.py 32 2**255-19
# 1. generates corner case test vectors in Python, writes them to edge.txt
# 2. uses scripts to generate C code for field arithmetic, automatically inserted into edge.c
# 3. Compiles and runs program edge.c, which reads from edge.txt and compares outputs against testvectors
#
# NOTE : Ensure generic=True in pseudo.py or monty.py (otherwise some test failures may occur)
#
# See main program below for corner case generation. Add your own!
#
import random
from contextlib import redirect_stdout
import sys
import subprocess

compiler="gcc"  # choose C compiler

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

def replaceback(namefile,newtext,oldfile):
    f = open(namefile,'r')
    filedata = f.read()
    f.close()

    nf = open(oldfile,'r')
    olddata = nf.read()
    nf.close

    newdata = filedata.replace(olddata,newtext)

    f = open(namefile,'w')
    f.write(newdata)
    f.close()

def bits(n) :
    b=0
    m=n
    while m!=0 :
        b+=1
        m>>=1
    return b

# return number of hex digits in n
def hexdigits(n) :
    b=bits(n)
    d=b//4
    if b%4 != 0 :
        d+=1
    if d%2==1 :   # make it even, for easy conversion to bytes
        d+=1
    return d

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

def tostring(a) :
    x='{:x}'.format(a)
    d=len(x)
    w=""
    while d<digits :
        w+='0'
        d+=1
    w+=x
    return w

# generate test vectors
def corner(a,b) :

    print(tostring(a))
    print(tostring(b))

    x=a%p
    print(tostring(x))

    x=(a+b)%p
    print(tostring(x))    # test modular addition

    x=(a-b)%p
    if x<0 :
        x+=p;
    print(tostring(x))    # test modular subtraction

    x=(b-a)%p              
    if x<0 :
        x+=p;
    print(tostring(x))    # test reverse modular subtraction

    x=(a*b)%p
    print(tostring(x))    # test modular multiplication

    x=(a*a)%p
    print(tostring(x))    # test modular squaring

    x=(2*a)%p
    print(tostring(x))

    x=(2*a+2*b)%p
    print(tostring(x))

    x=(2*a-2*b)%p
    if x<0 :
        x+=p;
    print(tostring(x))

    x=(2*b-2*a)%p
    if x<0 :
        x+=p;
    print(tostring(x))

    x=(2*a*2*b)%p
    print(tostring(x))

    x=(2*a*2*a)%p
    print(tostring(x))
   

if len(sys.argv)!=3 :
    print("Syntax error")
    print("Valid syntax - python edge.py <word length> <prime> OR <prime name>")
    print("For example - python edge.py 64 X25519")
    print("For example - python edge.py 64 2**255-19")
    exit(2)

WL=int(sys.argv[1])
if WL != 32 and WL !=64 :
    print("Only 32 and 64-bit word lengths supported")
    exit(2)

prime=sys.argv[2]

### Start of user editable area

# More named primes can be added here
p=0

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

if prime=="X25519" :
    p=2**255-19
    
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

if prime=="NIST256" :
    p=2**256-2**224+2**192+2**96-1

if prime=="NIST384" :
    p=2**384-2**128-2**96+2**32-1

if prime=="ED448" :
    p=2**448-2**224-1

if prime=="X448" :
    p=2**448-2**224-1

if prime=="GM270" :
    p=2**270-2**162-1

if prime=="GM240" :
    p=2**240-2**183-1

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

if prime=="NIST224" :
    p=2**224-2**96+1

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

if p==0 :
    if prime[0].isdigit() :
        p=eval(prime)
    else :
        print("This named modulus not supported (not a pseudo-Mersenne?)")
        exit(2)

n=p.bit_length() 
if n<120 or pow(3,p-1,p)!=1 :
    print("Not a sensible modulus, too small or not a prime")
    exit(2)

#Is it a pseudo-Mersenne?
c=2**n-p
PM=True
if c>=2**(WL-5) :
    PM=False  

# Generate finite field code for chosen modulus, and insert into edge.c

if PM:
    if WL==32: 
        subprocess.call("python3 pseudo.py 32 "+prime,shell=True)
    else :
        subprocess.call("python3 pseudo.py 64 "+prime,shell=True)
else :
    if WL==32: 
        subprocess.call("python3 monty.py 32 "+prime,shell=True)
    else :
        subprocess.call("python3 monty.py 64 "+prime,shell=True)

replacefromfile("edge.c","@field@","field.c")


digits=hexdigits(p)
r=random.randint(0,p-1)
i=inverse(r,p)

#Generate Test Vectors

#Add your own new corner cases here!
#n is number of bits in the prime modulus p. r is a random value, i is its inverse. c=2^n-p
#all entries *must* be positive and less than 2^n

#output testvectors to file edge.txt
with open('edge.txt', 'w') as f:
    with redirect_stdout(f):
        corner(p-1,p-1)
        corner(0,p-1)
        corner(0,0)
        corner(r,r)
        corner(r,i)
        corner(p,1)
        corner(p-1,1)
        corner(p-2,2)
        corner(p-1,r)
        corner(p-2,r)
        corner(2**64,2**64)
        corner(2**(n-1),2**(n-1)-1)
        corner(c,1)
        corner(c,2**n-1)
        corner(2**n-1,0)
        corner(2**n-1,1)
        corner(2**n-1,2**n-1)
    f.close()

# Compile and run C code to check outputs against generated test vectors
subprocess.call(compiler + " -O2 edge.c -o edge", shell=True)
subprocess.call("./edge", shell=True)

replaceback("edge.c","@field@","field.c") # put it back the way you found it!
