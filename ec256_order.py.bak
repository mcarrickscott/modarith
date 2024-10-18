q=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
print(hex(q))

c=(2**248)%q

R=2**256  # Montgomery reduction R value

c =(c*R)%q # convert to nresidue form

print(hex(c))

limbs=9    # 32 bit target
radix=29

#limbs=5     # 64 bit target
#radix=52

base=2**radix

print("const spint constant_c[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c%base)+",",end="")
    c//=base
print(hex(c),end="")
print("};")

