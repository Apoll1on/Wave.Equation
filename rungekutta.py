from typing import Callable
import numpy as np
from matplotlib import pyplot as plt


def RK4(f: Callable, u0, delt, t0, tmax, fileName: str):
    f = open(fileName, "a")
    f.write("New Run")
    iter = (tmax - t0) / delt

    u = u0
    t = t0

    for it in range(iter):
        k1 = f(t, u)
        k2 = f(t + delt / 2, u0 + delt * k1 / 2)
        k3 = f(t + delt / 2, u0 + delt * k2 / 2)
        k4 = f(t + delt, u0 + delt * k3)
        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
        t += delt

        f.write("Iteration: " + it + "   Time: " + t + "  u: " + u)

    f.close()


def firsteq(t, p):
    return p


def secondeq(t, x):
    return -x


def a6():
    fileName = "calculateddata.txt"

    x0 = 1
    p0 = 0
    tsteps = 7500  # 100 per 2pi is good
    delt = 10 * np.pi / tsteps
    tmax = 10 * np.pi
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
        k4x = firsteq(t + delt, p + delt * k3p)
        k4p = secondeq(t + delt, x + delt * k3x)

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

    cosarray = np.cos(tarray)
    sinarray = np.sin(tarray)

    fig1, ax1 = plt.subplots()
    ax1.plot(xarray, parray)

    fig2, ax2 = plt.subplots()
    ax2.plot(tarray, parray ** 2 + xarray ** 2)
    ax2.plot(tarray, xarray)
    ax2.plot(tarray, parray)

    fig3, ax3 = plt.subplots()
    ax3.plot(tarray, cosarray - xarray)
    ax3.plot(tarray, sinarray + parray)

    plt.show()
