import misc
import solver
import a10old
from matplotlib import pyplot as plt
import numpy as np


def calcplot():
    a = 4001
    periods=100
    c=periods*1000
    linestoread = [0]
    for i in range(1, periods):
        linestoread.append(i * 1000)
    xarray, times, phiarray, piarray = a10.solving(a, c + 2, alpha=1, boundaryCondition="extrapolation",
                                                   linestoread=linestoread)
    print(times)

    fig, ax = plt.subplots(2, 1)
    c = plt.get_cmap('gist_rainbow')
    n = len(linestoread)
    for i in range(n):
        # print("Phi " + str(i) + " :")
        # print(np.mean(phiarray[i]))
        # print("Pi " + str(i) + " :")
        # print(np.mean(piarray[i]))
        ax[0].plot(xarray, phiarray[i, :], label=format(float(times[i]), '.4f'), color = c(i/(n-1)))
        ax[1].plot(xarray, piarray[i, :], label=format(float(times[i]), '.4f'), color = c(i/(n-1)))

    ax[0].set_title('phi')
    ax[1].set_title('pi')
    ax[0].legend()
    ax[1].legend()
    plt.show()

def plotten(xpoints,xarray,linestoread,fileName):

    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)

    print(times)

    fig, ax = plt.subplots(2, 1)
    c = plt.get_cmap('gist_rainbow')
    n = len(linestoread)
    for i in range(n):
        # print("Phi " + str(i) + " :")
        # print(np.mean(phiarray[i]))
        # print("Pi " + str(i) + " :")
        # print(np.mean(piarray[i]))
        ax[0].plot(xarray, phiarray[i, :], label=format(float(times[i]), '.4f'), color=c(i / (n - 1)))
        ax[1].plot(xarray, piarray[i, :], label=format(float(times[i]), '.4f'), color=c(i / (n - 1)))

    ax[0].set_title('phi')
    ax[1].set_title('pi')
    ax[0].legend()
    ax[1].legend()
    plt.show()