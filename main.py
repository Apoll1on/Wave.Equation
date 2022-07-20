import numpy as np
from matplotlib import pyplot as plt

import a8
from funcsandder import *
import rungekutta as rk
import a7
import a5


# hh


def plotting(xvalues, mathfunc, derivative, numder, funcname):
    values = mathfunc(xvalues)
    fig, ax = plt.subplots(3)
    ax[0].plot(xvalues, derivative(xvalues))
    ax[1].plot(xvalues, numder(values))
    ax[2].plot(xvalues, numder(values) - derivative(xvalues))
    plt.title(funcname)
    plt.show()


# For Harmonic Oscillator we get:
# dx/dt=p firsteq
# dp/dt=-x by Newton Law
# Therefore runge kutta for both:



if __name__ == '__main__':
    # Aufgabe 6
    # rk.a6()
    # Aufgabe
    # a5.plot()
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
    #
    # aufgabe 7
    # a7.solving(1001, 1000, "calculateddata", "periodic")

    # a8

    a8.stabtest()