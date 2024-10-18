q = 0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409

print(hex(q))

c=(2**512)%q

WL=32

limbs=19
radix=28
R=2**532
if WL==64 :
	limbs=9
	radix=59
	R=2**531


c =(c*R)%q # convert to nresidue form
print(hex(c))
base=2**radix

print("const spint constant_c[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c%base)+",",end="")
    c//=base
print(hex(c),end="")
print("};")

