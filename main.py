import numpy as np
import math

if __name__ == '__main__':
    a = 0
    b = 1
    n = 20
    h = 1 / (n - 1)

    Omega = np.linspace(a, b, n)
    i_count = np.linspace(0, n - 1, n)


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
        return np.ei * hp(-a * i * h ** 2)

    def DUDX(u):
        diff = np.zeros_like(u)
        diff[1:-1] = (u[2:] - u[:-2]) / (2 *h)
        return diff

    def D2UDX2(u):
        diff = np.zeros_like(u)
        diff[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (h * h)
        return diff
