def bits(n) :
    b=0
    m=n
    while m!=0 :
        b+=1
        m>>=1
    return b

q = 39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319
q = q + 1 - 1388124618062372383606759648309780106643088307173319169677
print(hex(q))

bts=bits(q)
bytes=bts//8;
if bytes*8 != bts :
    bytes+=1
size=8*(bytes-1)

print("Split occurs at ",size)

#WL=32, from monty.py

limbs=14
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

limbs=7
radix=56

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
