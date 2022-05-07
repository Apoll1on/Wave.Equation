import numpy as np
import math

if __name__ == '__main__':
    a = 0
    b = 1
    n = 20
    h = 1 / (n - 1)

    Omega = np.linspace(a, b, n)
    i_count = np.linspace(0, n - 1, n)

    
        def DUDX(u):
        diff = np.zeros_like(u)
        diff[1:-1] = (u[2:] - u[:-2]) / (2 *h)
        return diff

    def D2UDX2(u):
        diff = np.zeros_like(u)
        diff[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (h * h)
        return diff
    
    
    def f1(i):
        return (i * h - 1 / 2) ** 2 + i * h


    def f2(i):
        (i * h - 1 / 2) ** 3 + (i * h - 1 / 2) ** 2 + i * h


    def f3(i):
        return np.sqrt(i * h)


    def s1(i):
        return np.sin(12 * np.pi * i * h)


    def s2(i):
        return np.sin(12 * np.pi * i * h) ** 4


    def g(i, a):
        return np.exp(-a * (i * h) ** 2)
    
    
    def Df1DX(i):
        return 2 * i * h
    
    def Df2DX(i):
        return 3 * (i * h)**2 - i * h + 3 / 4
    
    def Df3DX(i):
        return 0.5 / np.sqrt(i * h)
    
    def Ds1DX(i):
        return 12 * np.pi * np.cos(12 * np.pi * i * h)
    
    def Ds2DX(i):
        return 48 * np.pi * np.sin(12 * np.pi * i * h)**3 * np.cos(12 * np.pi * i * h)
    
    def DgDX(i, a):
        return -2 * a * i * h * np.exp(-a * (i * h)**2)
    
    afasfjlasf
    
