import solver
from matplotlib import pyplot as plt


def stabtest(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread):

    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, alpha,
                                                   phiinit, piinit, boundaryCondition, fileName, linestoread)
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


# plot q over time for every timestep
def convergence(x0, xmax, xpoints, t0, timesteps, alpha,
                                                   phiinit, piinit, boundaryCondition, fileName,linestoread,periods):
    timesteps = periods * (xpoints - 1)/alpha
    linestoread = [0]
    for i in range(1, periods):
        linestoread.append(int(i * (xpoints - 1)))
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps+2, alpha,
                                                   phiinit, piinit, boundaryCondition, fileName, linestoread)
    print("alpha",alpha,"\ntimes",times)
    fig, ax = plt.subplots(1)
    for i in range(len(linestoread)):
        ax.plot(xarray, phiarray[i, :], label=format(float(times[i]), '.4f'))

    ax.legend()
    ax.set_title("alpha"+format(float(alpha), '.2f'))
    plt.show()


def selfconvergence(x0, xmax, xpoints, t0, timesteps, alpha,
                                                   phiinit, piinit, boundaryCondition, fileName, linestoread):

    linestoread = [timesteps - 1]
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.25*alpha,
                                                   phiinit, piinit, boundaryCondition, fileName, linestoread)
    xarray2, times, phiarray2, piarray2 = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.5*alpha,
                                                   phiinit, piinit, boundaryCondition, fileName, linestoread)
    xarray4, times, phiarray4, piarray4 = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, alpha,
                                                   phiinit, piinit, boundaryCondition, fileName, linestoread)

    fig, ax = plt.subplots(1)
    ax.plot(xarray, (phiarray[0] - phiarray2[0][::2]), label='h - h/2')  #
    ax.plot(xarray, 16 * (phiarray2[0][::2] - phiarray4[0][::4]), label='h/2 - h/4')  #
    # ax3.plot(xarray, phiarray4[0][::4])

    ax.legend()
    plt.show()
