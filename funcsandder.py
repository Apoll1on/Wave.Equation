import numpy as np

a = 0
b = 1
n = 20
h = 1 / (n - 1)


def DUDX(u, h):
    diff = np.zeros_like(u)
    diff[1:-1] = (u[2:] - u[:-2]) / (2 * h)
    diff[0] = diff[1]
    diff[-1] = diff[-2]
    return diff


def D2UDX2(u, h):
    diff = np.zeros_like(u)
    diff[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (h * h)
    return diff


def f1(x):
    return (x - 1 / 2) ** 2 + x


def f2(x):
    return (x - 1 / 2) ** 3 + (x - 1 / 2) ** 2 + x


def f3(x):
    return np.sqrt(x)


def s1(x):
    return np.sin(12 * np.pi * x)


def s2(x):
    return np.sin(12 * np.pi * x) ** 4


def g(x, a):
    return np.exp(-a * x ** 2)


def Df1DX(x):
    return 2 * x


def Df2DX(x):
    return 3 * (x) ** 2 - x + 3 / 4


def Df3DX(x):
    return 0.5 / np.sqrt(x)


def Ds1DX(x):
    return 12 * np.pi * np.cos(12 * np.pi * x)


def Ds2DX(x):
    return 48 * np.pi * np.sin(12 * np.pi * x) ** 3 * np.cos(12 * np.pi * x)


def DgDX(x, a):
    return -2 * a * x * np.exp(-a * (x) ** 2)