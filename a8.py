# stability tests:

from a7 import solving
from matplotlib import pyplot as plt
import numpy as np


def stabtest():
    a = 101
    c = 1000
    linestoread = [int(0 * c), int(c * 0.125), int(c * 0.25), int(c * 0.375), int(c * 0.5), int(c * 0.625),
                   int(c * 0.75), int(c * 0.875), c]
    xarray, times, phiarray, piarray = solving(a, c + 1, linestoread=linestoread)

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
    linestoread = [
        c - 1]  # 0, int(c * 0.125), int(c * 0.25), int(c * 0.375), int(c * 0.5), int(c * 0.625), int(c * 0.75), int(c * 0.875),
    xarray, times, phiarray, piarray = solving(a, c + 1, linestoread=linestoread, alpha=.1)
    xarray2, times, phiarray2, piarray2 = solving(2 * a - 1, c + 1, linestoread=linestoread)
    xarray4, times, phiarray4, piarray4 = solving(4 * (a - 1) + 1, c + 1, linestoread=linestoread)

    fig3, ax3 = plt.subplots(1)
    ax3.plot(xarray, 8 * (phiarray[0] - phiarray2[0][::2]))  #
    ax3.plot(xarray, (phiarray2[0][::2] - phiarray4[0][::4]))  #
    # ax3.plot(xarray, phiarray4[0][::4])
    fig0, ax0 = plt.subplots()
    fig1, ax1 = plt.subplots()
    for i in range(len(linestoread)):
        print("Phi " + str(i) + " :")
        print(np.mean(phiarray[i]))
        print("Pi " + str(i) + " :")
        print(np.mean(piarray[i]))
        ax0.plot(xarray, phiarray[i, :])
        ax1.plot(xarray, piarray[i, :])

    print(times)

    # mylist=[]
    # for i in range(2000,c+1):
    #     if np.max(phiarray[i,:]) >= 0.999:
    #         mylist.append(i)
    #
    # print(mylist)
    # ma=10
    # for i in mylist:
    #     if np.sum(phiarray[i,1:-1]-phiarray[0,1:-1]) >= np.sum(phiarray[ma,1:-1]-phiarray[0,1:-1]):
    #         ma=i
    #
    # print(ma)

    plt.show()
