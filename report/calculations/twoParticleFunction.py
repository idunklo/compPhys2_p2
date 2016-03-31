from sympy.functions import exp
import sympy as s

alpha,beta,a,x1,x2,y1,y2 = s.symbols('alpha beta a x1 x2 y1 y2')

posDiff = s.sqrt((x1-x2)**2 + (y1-y2)**2)
Exp1 = exp(-alpha*(x1**2 + x2**2 + y1**2 + y2**2)/2)
Exp2 = exp(a*posDiff/(1+posDiff*beta))
ExpDiffX1 = s.diff(Exp1*Exp2,x1)
ExpDiffX2 = s.diff(Exp1*Exp2,x2)
ExpDiffY1 = s.diff(Exp1*Exp2,y1)
ExpDiffY2 = s.diff(Exp1*Exp2,y2)
#print("X1")
#print(s.simplify(ExpDiffX1))
#print("\nX2")
#print(s.simplify(ExpDiffX2))
#print("\nY1")
#print(s.simplify(ExpDiffY1))
#print("\nY2")
#print(s.simplify(ExpDiffY2))

firstDerivative = -alpha*(x1+x2+y1+y2)*Exp1*Exp2

print("\nX1X1")
print(s.diff(firstDerivative,x1))
print("\n\nX2X2")
print(s.diff(firstDerivative,x2))
print("\n\nY1Y1")
print(s.diff(firstDerivative,y1))
print("\n\nY2Y2")
print(s.diff(firstDerivative,y2))
