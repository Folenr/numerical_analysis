import numpy as np

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
