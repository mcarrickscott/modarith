def bits(n) :
    b=0
    m=n
    while m!=0 :
        b+=1
        m>>=1
    return b

q = 0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409
print(hex(q))

bts=bits(q)
bytes=bts//8;
if bytes*8 != bts :
    bytes+=1
size=8*(bytes-1)

print("Split occurs at ",size)

#WL=32, from monty.py

limbs=19
radix=28

excess=limbs*radix-bts
if excess==0 :
    print("Special case")
    exit(0)

R=2**(limbs*radix)

c=(2**size)%q
c =(c*R)%q # convert to nresidue form
#print(hex(c))
base=2**radix

print("static const spint constant_c[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c%base)+",",end="")
    c//=base
print(hex(c),end="")
print("};")

#WL=64, from monty.py

limbs=9
radix=59

excess=limbs*radix-bts
if excess==0 :
    print("Special case")
    exit(0)

R=2**(limbs*radix)

c=(2**size)%q
c =(c*R)%q # convert to nresidue form
#print(hex(c))
base=2**radix

print("static const spint constant_c[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c%base)+",",end="")
    c//=base
print(hex(c),end="")
print("};")

