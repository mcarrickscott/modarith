q = 39402006196394479212279040100143613805079739270465446667948293404245721771496870329047266088258938001861606973112319
q = q + 1 - 1388124618062372383606759648309780106643088307173319169677

print(hex(q))

c=(2**376)%q

WL=32

limbs=14
radix=28
R=2**392
if WL==64 :
	limbs=7
	radix=56
	R=2**392


c =(c*R)%q # convert to nresidue form
print(hex(c))
base=2**radix

print("const spint constant_c[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(c%base)+",",end="")
    c//=base
print(hex(c),end="")
print("};")

