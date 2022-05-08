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


# For Harmonic Oscillator we get:
# dx/dt=p firsteq
# dp/dt=-x by Newton Law
# Therefore runge kutta for both:
def firsteq(t, p):
    return p


def secondeq(t, x):
    return -x


if __name__ == '__main__':

    # Aufgabe 6
    fileName = "calculateddata.txt"

    x0 = 1
    p0 = 0
    tsteps = 100
    delt = 2 * np.pi / tsteps
    tmax = 2 * np.pi
    t0 = 0

    f = open(fileName, "a")
    f.write("\n\nNew Run")

    x = x0
    p = p0
    t = t0
    iter = 0
    tarray = np.zeros(tsteps + 1)
    xarray = np.zeros(tsteps + 1)
    parray = np.zeros(tsteps + 1)

    while t < tmax:
        k1x = firsteq(t, p)
        k1p = secondeq(t, x)
        k2x = firsteq(t + delt / 2, p + delt * k1p / 2)
        k2p = secondeq(t + delt / 2, x + delt * k1x / 2)
        k3x = firsteq(t + delt / 2, p + delt * k2p / 2)
        k3p = secondeq(t + delt / 2, x + delt * k2x / 2)
        k4x = firsteq(t + delt, p + k3p)
        k4p = secondeq(t + delt, x + k3x)

        x = x + delt * (k1x / 6 + k2x / 3 + k3x / 3 + k4x / 6)
        p = p + delt * (k1p / 6 + k2p / 3 + k3p / 3 + k4p / 6)

        xarray[iter] = x
        parray[iter] = p
        tarray[iter] = t

        f.write("\n")
        f.write("Time: " + str(t) + "  p: " + str(p) + "  x: " + str(x))

        t += delt
        iter += 1

    f.close()

    fig, ax = plt.subplots(2)
    ax[0].plot(tarray, xarray)
    ax[1].plot(tarray, parray)

    fig1, ax1 = plt.subplots()
    ax1.plot(xarray, parray)
    plt.show()

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
