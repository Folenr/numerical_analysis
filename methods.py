import numpy as np
from sympy import Symbol, diff, lambdify

class Methods:
    def biSection(f,p1,p2,tol):
        if np.sign(f(p1)) == np.sign(f(p2)):
            raise Exception("The scalars a and b do not bound a root")
        m = (p1 + p2)/2   
        if np.abs(f(m)) < tol:
            return m
        elif np.sign(f(p1)) == np.sign(f(m)):
            return Methods.biSection(f, m, p2, tol)
        elif np.sign(f(p2)) == np.sign(f(m)):
            return Methods.biSection(f, p1, m, tol)
    
    def fixedPoint(f,g,p,tol):
        if np.abs(f(p)) < tol:
            return p
        p = g(p)
        return Methods.fixedPoint(f, g, p, tol)
    
    def newton(f,g,p,tol):
        if g == None:
            x = Symbol('x')
            f_sympy = f(x)
            g = diff(f_sympy, x)
            g = lambdify(x, g)
        if np.abs(f(p)) < tol:
            return p
        p = p - f(p)/g(p)
        return Methods.newton(f, g, p, tol)
    
    def secant(f,p1,p2,tol):
        if np.abs(f(p1)) < tol:
            return p1
        temp = p1
        p1 = p1 - (f(p1)*(p1-p2))/(f(p1)-f(p2))
        p2 = temp
        return Methods.secant(f, p1, p2, tol)

f = lambda x: x**3 + 4*x**2 - 10
g = lambda x: 0.5*(10-x**3)**0.5

print(Methods.biSection(f,0,2,0.00000000001))
print(Methods.fixedPoint(f,g,1.5,0.00000000001))
print(Methods.newton(f,None,1.5,0.00000000001))
print(Methods.secant(f,2,1,0.00000000001))
