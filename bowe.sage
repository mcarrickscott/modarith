# Adapted by M. Scott from script by
# Authors: Sean Bowe, Alessandro Chiesa, Matthew Green, Ian Miers, Pratyush Mishra, Howard Wu
# Finds a good strong Montgomery curve

from sage.all import *

# return True for good curve, else False
def IsGoodMontgomeryCurve(F,A,B):
    a=(3 - A^2)/(3 * B^2)		# convert to Weierstrass form for Sage
    b=(2 * A^3 - 9 * A)/(27 * B^3)
    ec = EllipticCurve(F,[0,0,0,a,b])   # set curve
    order = ec.order() # order of curve ...
    orderm8 = order % 8

    if orderm8==4 :
    	order//=4
    if orderm8==0 :
    	order//=8
   
    print("order: %s" % hex(order))
    if not (is_prime(order)): return False # only consider curves with co-factors 4 or 8
    print ("Found secure Montgomery curve with A=%s and B=%s" % (A, B))
    if orderm8 == 0:
        print ("cofactor = 8")
    elif orderm8 == 4:
        print ("cofactor = 4")
    return True

# set field, initial value for A >= 3. Fix B=1
# looking for curve y^2=x^3+Ax^2+x mod q
# incrementing A in search for good curve of order 4 or 8 times a prime
# NOT twist secure (no need for this condition IMHO)

q=2**255-19
A=6   # SHOULD be 2 mod 4

#q=2**266-3
#A=6

#q=2**266-2**168-1
#A=6

B=1
Fq = GF(q)
while true:
	print("Montgomery A= ",A);
	if IsGoodMontgomeryCurve(Fq, Fq(A), Fq(B)) :
		break;
	A=A+4  # increment A and try again. We want A+2 to be divisible by 4
