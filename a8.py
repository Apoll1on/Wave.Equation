# stability tests:

from a7 import solving
from matplotlib import pyplot as plt
import numpy as np


def stabtest():
    a = 201
    c = 100
    linestoread = [c] # int(0 * c), int(c * 0.125), int(c * 0.25), int(c * 0.375), int(c * 0.5), int(c * 0.625), int(c * 0.75), int(c * 0.875),
    xarray, times, phiarray, piarray = solving(a, c + 1, linestoread=linestoread, alpha=0.1)
    xarray2, times2, phiarray2, piarray2 = solving(2 * a - 1, c + 1, linestoread=linestoread, alpha=0.2)
    #xarray4, times4, phiarray4, piarray4 = solving(4 * (a - 1) + 1, c + 1, linestoread=linestoread, alpha=0.4)

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

    # fig3, ax3 = plt.subplots(1)
    # for i in range(len(linestoread)):
    #     ax3.plot(xarray, phiarray[i, :])

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

    fig3, ax3 = plt.subplots(1)
    ax3.plot(xarray, phiarray[0] ) #- phiarray2[0][::2]
    #ax3.plot(xarray, phiarray2[0][::2] - phiarray4[0][::4]) #
    #ax3.plot(xarray, phiarray4[0][::4])


    plt.show()


#def selfconvergence():