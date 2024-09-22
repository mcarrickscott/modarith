q=(2**448-2**224-28312320572429821613362531907042076847709625476988141958474579766324)//4
print(hex(q))

c1=(2**472)%q
c2=2**440

R=2**448  # Montgomery reduction R value

c1 =(c1*R)%q # convert to nresidue form
c2= (c2*R)%q

print(hex(c1))
print(hex(c2))

limbs=16    # 32 bit target
radix=28

limbs=8     # 64 bit target
radix=56

base=2**radix

print("const spint constant_c1[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c1%base)+",",end="")
    c1//=base
print(hex(c1),end="")
print("};")

print("const spint constant_c2[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c2%base)+",",end="")
    c2//=base
print(hex(c2),end="")
print("};")
