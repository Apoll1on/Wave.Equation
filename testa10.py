import misc
import a10
from matplotlib import pyplot as plt
import time

def calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread):
    start_time = time.time()
    xarray, times, phiarray, piarray = a10.solving(x0,xmax,xpoints,t0,timesteps+2,alpha,
                                                   phiinit,piinit,boundaryCondition,fileName,linestoread)
    print("--- %s seconds ---" % (time.time() - start_time))
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