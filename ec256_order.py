q=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
print(hex(q))

c=(2**248)%q

WL=32

limbs=9
radix=29
R=2**261
if WL==64 :
	limbs=5
	radix=52
	R=2**260


c =(c*R)%q # convert to nresidue form
print(hex(c))
base=2**radix

print("const spint constant_c[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c%base)+",",end="")
    c//=base
print(hex(c),end="")
print("};")

