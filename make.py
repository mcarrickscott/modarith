# Tool for generating large constants for legacy curves
# python make.py monty 2**448-2**224-1 56 0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa955555555555555555555555555555555555555555555555555555555

import sys

if len(sys.argv)!=5 :
    print("Syntax error")
    print("Valid syntax - python make.py <type> <prime> <radix> <constant>")
    print("For example - python make.py monty 2**255-19 51 0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa955555555555555555555555555555555555555555555555555555555")
    exit(0);

type=sys.argv[1]
prime=eval(sys.argv[2])
radix=int(sys.argv[3])
constant=eval(sys.argv[4])

bts=prime.bit_length()

limbs=int(bts/radix)
if bts%radix != 0 :
    limbs+=1

R=1
if type=="monty" :
    if bts%radix==0 :
        R=2**((limbs+1)*radix)
    else :
        R=2**(limbs*radix)

constant=(constant*R)%prime

base=2**radix
print("#if Radix == ",radix)
print("const spint constant[{}]={{".format(limbs),end="")
for i in range(0,limbs-1) :
    print(hex(constant%base)+",",end="")
    constant//=base
print(hex(constant),end="")
print("};")
print("#endif")

