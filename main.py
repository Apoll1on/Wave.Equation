import numpy as np
import math
from matplotlib import pyplot as plt


def plotting(xvalues, mathfunc, derivative, numder, funcname):
    values = mathfunc(xvalues)
    fig, ax = plt.subplots(3)
    ax[0].plot(xvalues, derivative(xvalues))
    ax[1].plot(xvalues, numder(values))
    ax[2].plot(xvalues, numder(values) - derivative(xvalues))
    plt.title(funcname)
    plt.show()


def DUDX(u):
    diff = np.zeros_like(u)
    diff[1:-1] = (u[2:] - u[:-2]) / (2 * h)
    return diff


def D2UDX2(u):
    diff = np.zeros_like(u)
    diff[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (h * h)
    return diff


def f1(x):
    return (x - 1 / 2) ** 2 + x


def f2(x):
    (x - 1 / 2) ** 3 + (x - 1 / 2) ** 2 + x


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


if __name__ == '__main__':
    a = 0
    b = 1
    n = 20
    h = 1 / (n - 1)

    Omega = np.linspace(a, b, n)
    i_count = np.linspace(0, n - 1, n)

    plotting(Omega, f1, Df1DX, DUDX, "f1")
    # print(Omega)
    # print(f1(Omega))
    # print(DUDX(f1(Omega)))
    # print(Df1DX(Omega))
