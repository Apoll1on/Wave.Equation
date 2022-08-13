import misc
import a10ram
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)
import time
import a10ram
import numpy as np

divider = 100

def calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread):
    start_time = time.time()
    xarray, times, phiarray, piarray = a10ram.solving(x0,xmax,xpoints,t0,timesteps+2,alpha,
                                                   phiinit,piinit,boundaryCondition,fileName,linestoread)
    print("--- %s seconds ---" % (time.time() - start_time))
    print(times)
    # fig,ax=plt.subplots()
    # ax.plot(xarray,0.15 / ((np.cosh(0.18 * (xarray) + 0.43)) * (np.cosh(0.18 * (xarray) + 0.43))))
    # plt.show()
    # fig, ax = plt.subplots(2, 1)
    # c = plt.get_cmap('gist_rainbow')
    # n = len(linestoread)
    # ax[0].plot(xarray, 0.15 / ((np.cosh(0.18 * (xarray) + 0.43)) * (np.cosh(0.18 * (xarray) + 0.43))),color="black")
    # for i in range(int(n/divider)):
    #     # print("Phi " + str(i) + " :")
    #     # print(np.mean(phiarray[i]))
    #     # print("Pi " + str(i) + " :")
    #     # print(np.mean(piarray[i]))
    #     ax[0].plot(xarray, phiarray[divider*i, :]+divider*i*a10ram.PTpotential(xarray), label=format(float(times[divider*i]), '.4f'), color = c(divider*i/(n-1)))
    #     ax[1].plot(xarray, piarray[divider*i, :], label=format(float(times[divider*i]), '.4f'), color = c(divider*i/(n-1)))
    # ax[0].set_title('phi')
    # ax[1].set_title('pi')
    # ax[0].legend()
    # ax[1].legend()
    # plt.show()

    fig, ax = plt.subplots()
    c = plt.get_cmap('gist_rainbow')
    n = len(linestoread)
    ax.plot(xarray, 0.15 / ((np.cosh(0.18 * (xarray) + 0.43)) * (np.cosh(0.18 * (xarray) + 0.43))), color="black")
    for i in range(int(n / divider)):
        ax.plot(xarray, phiarray[divider * i, :],
                label=format(float(times[divider * i]), '.4f'), color=c(divider * i / (n - 1)))
    ax.legend()
    ax.set_title('phi')
    plt.show()

    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax2.set_ylim(pow(10, -7), pow(10, 0))
    ax.plot(times[::], phiarray[::, int(xpoints * 0.75)], color="blue")
    ax2.semilogy(times[::], np.absolute(phiarray[::, int(xpoints * 0.75)]), color="red")
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.0f'))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.xaxis.set_minor_formatter(FormatStrFormatter('% 1.0f'))
    ax.set_xlim(180, 550)
    ax.set_xlabel("Time")
    ax.set_ylabel(r'$\Phi$', color="blue")
    ax2.set_ylabel(r'|$\Phi$|', color="red")
    plt.show()
    print(times)

def plotten(xpoints,xarray,linestoread,fileName):
    times, phiarray, piarray = misc.readfromshortfile(fileName, xpoints, lines=linestoread)
    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    ax2.set_ylim(pow(10, -7), pow(10, 0))
    ax.plot(times[::], phiarray[::, int(xpoints * 0.75)], color="blue")
    ax2.semilogy(times[::], np.absolute(phiarray[::, int(xpoints * 0.75)]), color="red")
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.0f'))
    ax.xaxis.set_minor_locator(MultipleLocator(50))
    ax.xaxis.set_minor_formatter(FormatStrFormatter('% 1.0f'))
    ax.set_xlim(180, 550)
    ax.set_xlabel("Time")
    ax.set_ylabel(r'$\Phi$', color="blue")
    ax2.set_ylabel(r'|$\Phi$|', color="red")
    plt.show()
    print(times)

    fig, ax = plt.subplots()
    c = plt.get_cmap('gist_rainbow')
    n = len(linestoread)
    ax.plot(xarray, 0.15 / ((np.cosh(0.18 * (xarray) + 0.43)) * (np.cosh(0.18 * (xarray) + 0.43))), color="black")
    for i in range(int(n / divider)):
        # print("Phi " + str(i) + " :")
        # print(np.mean(phiarray[i]))
        # print("Pi " + str(i) + " :")
        # print(np.mean(piarray[i]))
        ax.plot(xarray, phiarray[divider * i, :],
                label=format(float(times[divider * i]), '.4f'), color=c(divider * i / (n - 1)))

    ax.set_title('phi')
    ax.legend()

    plt.show()