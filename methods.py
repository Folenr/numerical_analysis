import numpy as np
from sympy import Symbol, diff, lambdify

class Methods:
    def biSection(f,p1,p2,tol,n):
        if np.sign(f(p1)) == np.sign(f(p2)):
            raise Exception("The scalars a and b do not bound a root")
        m = (p1 + p2)/2   
        if np.abs(f(m)) < tol:
            return m, n
        elif np.sign(f(p1)) == np.sign(f(m)):
            n+=1
            return Methods.biSection(f, m, p2, tol,n)
        elif np.sign(f(p2)) == np.sign(f(m)):
            n+=1
            return Methods.biSection(f, p1, m, tol,n)
    
    def fixedPoint(f,g,p,tol,n):
        if np.abs(f(p)) < tol:
            return p, n
        n+=1
        p = g(p)
        return Methods.fixedPoint(f, g, p, tol,n)
    
    def newton(f,g,p,tol,n):
        if g == None:
            x = Symbol('x')
            f_sympy = f(x)
            g = diff(f_sympy, x)
            g = lambdify(x, g)
        if np.abs(f(p)) < tol:
            return p, n
        n+=1
        p = p - f(p)/g(p)
        return Methods.newton(f, g, p, tol,n)
    
    def secant(f,p1,p2,tol,n):
        if np.abs(f(p1)) < tol:
            return p1, n
        n+=1
        temp = p1
        p1 = p1 - (f(p1)*(p1-p2))/(f(p1)-f(p2))
        p2 = temp
        return Methods.secant(f, p1, p2, tol,n)
    
    def falsePosition (f,p1,p2,tol,n):
        if np.abs(f(p1)) < tol:
            return p1, n
        n+=1
        temp = p1
        p1 = p1 - (f(p1)*(p1-p2))/(f(p1)-f(p2))
        if np.sign(f(p1)) != np.sign(f(temp)):
            p2 = temp
        return Methods.secant(f, p1, p2, tol,n)
    
    def modifiedNewton(f,f1,f2,p,tol,n):
        if f1 == None:
            x = Symbol('x')
            f_sympy = f(x)
            f1 = diff(f_sympy, x)
            f1 = lambdify(x, f1)
        if f2 == None:
            x = Symbol('x')
            f_sympy = f1(x)
            f2 = diff(f_sympy, x)
            f2 = lambdify(x, f2)
        if np.abs(f(p)) < tol:
            return p, n
        n+=1
        p = p - (f(p)*f1(p))/(f1(p)*f1(p) - f(p)*f2(p))
        return Methods.modifiedNewton(f, f1,f2, p, tol,n)

f = lambda x: x**3 + 4*x**2 - 10
g = lambda x: 0.5*(10-x**3)**0.5

print("biSection      : ",Methods.biSection(f,0,2,0.00000000001,0))
print("fixedPoint     : ",Methods.fixedPoint(f,g,1.5,0.00000000001,0))
print("newton         : ",Methods.newton(f,None,1.5,0.00000000001,0))
print("secant         : ",Methods.secant(f,2,1,0.00000000001,0))
print("falsePosition  : ",Methods.falsePosition(f,2,1,0.00000000001,0))
print("modifiedNewton : ",Methods.modifiedNewton(f,None,None,1,0.00000000001,0))
