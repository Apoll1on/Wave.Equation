# stability tests:

import solver
import a7
from matplotlib import pyplot as plt
import numpy as np


def stabtest():
    a = 501
    c = 13721
    linestoread = [int(0 * c), int(c * 0.125), int(c * 0.25), int(c * 0.375), int(c * 0.5), int(c * 0.625),
                   int(c * 0.75), int(c * 0.875), c]
    xarray, times, phiarray, piarray = solver.solving(a, c + 2, linestoread=linestoread,
                                                      boundaryCondition="extrapolation")

    print(times)
    fig, ax = plt.subplots(2, 1)
    for i in range(len(linestoread)):
        # print("Phi " + str(i) + " :")
        # print(np.mean(phiarray[i]))
        # print("Pi " + str(i) + " :")
        # print(np.mean(piarray[i]))
        ax[0].plot(xarray, phiarray[i, :], label=format(float(times[i]), '.4f'))
        ax[1].plot(xarray, piarray[i, :], label=format(float(times[i]), '.4f'))

    ax[0].set_title('phi')
    ax[1].set_title('pi')
    plt.show()


def convergence():
    periods = 10
    a = 201
    c = periods * (a - 1)
    linestoread = [0]
    for i in range(1, periods):
        linestoread.append(int(i * (a - 1)))
    xarray, times, phiarray, piarray = solver.solving(a, c + 2, linestoread=linestoread, alpha=1)
    print(times)
    fig, ax = plt.subplots(1)
    for i in range(len(linestoread)):
        ax.plot(xarray, phiarray[i, :], label=format(float(times[i]), '.4f'))

    ax.legend()
    plt.show()


def selfconvergence():
    a = 101
    c = 1000
    linestoread = [
        c - 1]  # 0, int(c * 0.125), int(c * 0.25), int(c * 0.375), int(c * 0.5), int(c * 0.625), int(c * 0.75), int(c * 0.875),
    xarray, times, phiarray, piarray = solver.solving(a, c + 2, linestoread=linestoread, alpha=.1)
    xarray2, times, phiarray2, piarray2 = solver.solving(2 * a - 1, c + 2, linestoread=linestoread)
    xarray4, times, phiarray4, piarray4 = solver.solving(4 * (a - 1) + 1, c + 2, linestoread=linestoread)

    fig, ax = plt.subplots(1)
    ax.plot(xarray, (phiarray[0] - phiarray2[0][::2]), label='h - h/2')  #
    ax.plot(xarray, 16 * (phiarray2[0][::2] - phiarray4[0][::4]), label='h/2 - h/4')  #
    # ax3.plot(xarray, phiarray4[0][::4])
    ax.legend()
    plt.show()
