from sympy.functions import exp
import sympy as s

Hy0,Hy1,Hy2,Hy3,Hx0,Hx1,Hx2,Hx3,phi, omega,x,y =\
        s.symbols('Hy0 Hy1 Hy2 Hy3 Hx0 Hx1 Hx2 Hx3 phi omega x y')
Hx0 = 1
Hx1 = 2*x*s.sqrt(omega)
Hx2 = 4*x*s.sqrt(omega)*x*s.sqrt(omega) - 2
Hx3 = 8*x*s.sqrt(omega)*x*s.sqrt(omega)*x*s.sqrt(omega) - 12*x*s.sqrt(omega)
Hy0 = 1
Hy1 = 2*y*s.sqrt(omega)
Hy2 = 4*y*s.sqrt(omega)*y*s.sqrt(omega) - 2
Hy3 = 8*y*s.sqrt(omega)*y*s.sqrt(omega)*y*s.sqrt(omega) - 12*y*s.sqrt(omega)

phi  = exp(-omega*(x*x+y*y)/2)
phi0 = phi*Hx0*Hy0
phi1 = phi*Hx1*Hy1
phi2 = phi*Hx2*Hy2
phi3 = phi*Hx3*Hy3


print("Grad H0")
s.pretty_print(((s.diff(phi0,x)/phi)).simplify().simplify())
print("\nGrad H1")
s.pretty_print(((s.diff(phi1,x)/phi)).simplify().simplify())
print("\nGrad H2")
s.pretty_print(((s.diff(phi2,x)/phi)).simplify().simplify())
print("\nGrad H2")
s.pretty_print(((s.diff(phi3,x)/phi)).simplify().simplify())
print("-----------------------------------------------------------------------")
print("Lap H0")
s.pretty_print(((s.diff(s.diff(phi0,x),x)/phi)+s.diff(s.diff(phi0,y),y)/phi).simplify().simplify())
print("\nLap H1")
s.pretty_print(((s.diff(s.diff(phi1,x),x)/phi)+s.diff(s.diff(phi1,y),y)/phi).simplify().simplify()) 
print("\nLap H2")
s.pretty_print(((s.diff(s.diff(phi2,x),x)/phi)+s.diff(s.diff(phi2,y),y)/phi).simplify().simplify())
print("\nLap H3")
s.pretty_print(((s.diff(s.diff(phi3,x),x)/phi)+s.diff(s.diff(phi3,y),y)/phi).simplify().simplify())
