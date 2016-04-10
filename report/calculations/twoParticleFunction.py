from sympy.functions import exp
import sympy as s

alpha,beta,a,x1,x2,y1,y2,r12,r1,r2 = s.symbols('alpha beta a x1 x2 y1 y2 r12 r1 r2')

#r12 = s.sqrt((x1-x2)**2)
#r1  = x1**2
#r2  = x2**2
r12 = s.sqrt((x1-x2)**2 + (y1-y2)**2)
r1  = x1**2 + y1**2
r2  = x2**2 + y2**2
sep = x1 -x2
obr12 = 1+beta*r12
Exp1 = exp(-alpha*(r1+r2)/2)
Exp2 = exp(a*r12/(1+r12*beta))
s.pretty_print(s.diff(Exp1,x1))
print("\n\n\n")
s.pretty_print(s.diff(Exp2,x1))
#s.pretty_print(1/r12)
#print("\n")
#s.pretty_print(sep/obr12**2*s.diff(1/r12,x1))
#print("\n")
#s.pretty_print(1/obr12**2)
#print("\n")
#s.pretty_print(sep/r12*s.diff(1/obr12**2,x1))
#ExpDiffX1 = s.diff(Exp1*Exp2,x1)
#ExpDiffY1 = s.diff(Exp1*Exp2,y1)
#
#secxx = s.diff(ExpDiffX1,x1)
#secxy = s.diff(ExpDiffX1,y1)
#secyx = s.diff(ExpDiffY1,x1)
#secyy = s.diff(ExpDiffY1,y1)
#print("\nX1X1")
#e1e2 = s.simplify(Exp1*Exp2)
#print(s.simplify(s.simplify(secxx)/e1e2))
#print("\n\n")
#print(secxy)
#print("\n\n")
#print(secyx)
#print("\n\n")
#print(secyy)

#s.pretty_print(s.simplify(s.diff(s.simplify(s.diff(s.simplify(Exp2),x1)),x1)))
