# stability tests:

from a7 import solving
from matplotlib import pyplot as plt
import numpy as np


def stabtest():
    a = 101
    c = 1000
    linestoread = [0, int(0.125 * c), int(0.25 * c), int(0.375 * c), int(0.5 *c), int(0.625 * c), int(0.75 * c), int(0.875 * c), c -1]
    xarray, times, phiarray, piarray = solving(a, c + 1, linestoread=linestoread, alpha = .1)

    fig, ax = plt.subplots(1)
    for i in range(len(linestoread)):
        ax.plot(xarray, phiarray[i, :], label = i)

    ax.legend()
    plt.show()



def convergence():
    a = 101
    c = 3500
    linestoread = [0, 999, 1999, 2999, c - 1]
    xarray, times, phiarray, piarray = solving(a, c + 1, linestoread=linestoread, alpha = .1)

    fig, ax = plt.subplots(1)
    for i in range(len(linestoread)):
        ax.plot(xarray, phiarray[i, :], label = i)

    ax.legend()
    plt.show()


def selfconvergence():
    a = 101
    c = 1000
    linestoread = [c - 1] #0, int(c * 0.125), int(c * 0.25), int(c * 0.375), int(c * 0.5), int(c * 0.625), int(c * 0.75), int(c * 0.875),
    xarray, times, phiarray, piarray = solving(a, c + 1, linestoread=linestoread, alpha = .1)
    xarray2, times, phiarray2, piarray2 = solving(2 * a - 1, c + 1, linestoread=linestoread)
    xarray4, times, phiarray4, piarray4 = solving(4 * (a - 1) + 1, c + 1, linestoread=linestoread)


    fig3, ax3 = plt.subplots(1)
    ax3.plot(xarray, 8 * (phiarray[0] - phiarray2[0][::2])) #
    ax3.plot(xarray, (phiarray2[0][::2] - phiarray4[0][::4])) #
    # ax3.plot(xarray, phiarray4[0][::4])


    plt.show()