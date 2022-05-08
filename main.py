import numpy as np
from matplotlib import pyplot as plt
from funcsandder import *
import rungekutta as rk


def plotting(xvalues, mathfunc, derivative, numder, funcname):
    values = mathfunc(xvalues)
    fig, ax = plt.subplots(3)
    ax[0].plot(xvalues, derivative(xvalues))
    ax[1].plot(xvalues, numder(values))
    ax[2].plot(xvalues, numder(values) - derivative(xvalues))
    plt.title(funcname)
    plt.show()


if __name__ == '__main__':

# Aufgabe 5
# taking n,a,b,h from funcsandder
# ugly but works...

# Omega = np.linspace(a, b, n)
# i_count = np.linspace(0, n - 1, n)

# plotting(Omega, f1, Df1DX, DUDX, "f1")
# plotting(Omega, f2, Df2DX, DUDX, "f2")
# print(Omega)
# print(f1(Omega))
# print(DUDX(f1(Omega)))
# print(Df1DX(Omega))
