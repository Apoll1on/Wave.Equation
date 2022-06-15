# stability tests:

from a7 import solving
from matplotlib import pyplot as plt
import numpy as np


def stabtest():
    a = 4
    b = 0.001
    c = 3
    xarray, times, phiarray, piarray = solving(a, b, c + 1)

    # fig1, ax1 = plt.subplots(3, 3)
    # ax1[0, 0].plot(xarray, piarray[0, 1:-1])
    # ax1[0, 1].plot(xarray, piarray[12, 1:-1])
    # ax1[0, 2].plot(xarray, piarray[25, 1:-1])
    # ax1[1, 0].plot(xarray, piarray[37, 1:-1])
    # ax1[1, 1].plot(xarray, piarray[50, 1:-1])
    # ax1[1, 2].plot(xarray, piarray[62, 1:-1])
    # ax1[2, 0].plot(xarray, piarray[75, 1:-1])
    # ax1[2, 1].plot(xarray, piarray[87, 1:-1])
    # ax1[2, 2].plot(xarray, piarray[99, 1:-1])
    #
    # fig2, ax2 = plt.subplots(3, 3)
    # ax2[0, 0].plot(xarray, phiarray[0, 1:-1])
    # ax2[0, 1].plot(xarray, phiarray[12, 1:-1])
    # ax2[0, 2].plot(xarray, phiarray[25, 1:-1])
    # ax2[1, 0].plot(xarray, phiarray[37, 1:-1])
    # ax2[1, 1].plot(xarray, phiarray[50, 1:-1])
    # ax2[1, 2].plot(xarray, phiarray[62, 1:-1])
    # ax2[2, 0].plot(xarray, phiarray[75, 1:-1])
    # ax2[2, 1].plot(xarray, phiarray[87, 1:-1])
    # ax2[2, 2].plot(xarray, phiarray[99, 1:-1])

    fig3, ax3 = plt.subplots(1)
    ax3.plot(xarray, phiarray[int(0 * c), :])
    ax3.plot(xarray, phiarray[int(c * 0.125), :])
    ax3.plot(xarray, phiarray[int(c * 0.25), :])
    ax3.plot(xarray, phiarray[int(c * 0.375), :])
    ax3.plot(xarray, phiarray[int(c * 0.5), :])
    ax3.plot(xarray, phiarray[int(c * 0.625), :])
    ax3.plot(xarray, phiarray[int(c * 0.75), :])
    ax3.plot(xarray, phiarray[int(c * 0.875), :])
    ax3.plot(xarray, phiarray[c, :])

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
