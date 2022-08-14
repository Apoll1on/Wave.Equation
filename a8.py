import numpy as np
import funcsandder
import solver
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)
import misc
import solvernew


def selfconv_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    func_phiint = funcsandder.gausswave
    func_piint = funcsandder.dergaus

    linestoread = [timesteps - 1]
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                      func_phiint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                      func_piint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                      boundaryCondition, fileName,
                                                      linestoread)
    xarray2, times, phiarray2, piarray2 = solver.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                         func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5, 0.05),
                                                         func_piint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5, 0.05),
                                                         boundaryCondition, fileName, linestoread)
    xarray4, times, phiarray4, piarray4 = solver.solving(x0, xmax, 4 * (xpoints - 1) + 1, t0, timesteps + 2, 0.4,
                                                         func_phiint(np.linspace(x0, xmax, 4 * (xpoints - 1) + 1), 0.5, 0.05),
                                                         func_piint(np.linspace(x0, xmax, 4 * (xpoints - 1) + 1), 0.5, 0.05),
                                                         boundaryCondition, fileName, linestoread)

    fig, ax = plt.subplots(1)
    ax.plot(xarray, (phiarray[0] - phiarray2[0][::2]), label='$\phi\'^{(h)} - \phi\'^{(h/2)}$')
    ax.plot(xarray, 4 * (phiarray2[0][::2] - phiarray4[0][::4]), label='$4 (\phi\'^{(h/2)} - \phi\'^{(h/4)})$',
            linestyle='dashed')
    # ax.plot(xarray, (phiarray[0]), label='$\phi\'^{(h)} - \phi\'^{(h/2)}$')
    # ax.plot(xarray, func_phiint(np.linspace(x0, xmax, xpoints), 0.5, 0.05), label='$4 (\phi\'^{(h/2)} - \phi\'^{(h/4)})$', linestyle='dashed')

    ax.legend(loc='upper left')
    plt.show()


def plotconvfd(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    xarray = np.linspace(x0, xmax, xpoints)
    xarray2 = np.linspace(x0, xmax, 2 * xpoints - 1)
    times, phiarray, piarray = misc.readdata("calch.txt", xpoints, lines=linestoread)
    times2, phiarray2, piarray2 = misc.readdata("calch2.txt", 2 * xpoints - 1, lines=linestoread)
    qarray = np.zeros(len(times))
    fig, ax = plt.subplots(1)

    for i in range(len(times)):
        top = np.sqrt(np.sum((funcsandder.gausswave_verschoben(np.linspace(x0, xmax, xpoints), 0.5, 0.05,
                                                               float(times[i])) - phiarray[i]) *
                             (funcsandder.gausswave_verschoben(np.linspace(x0, xmax, xpoints), 0.5, 0.05,
                                                               float(times[i])) - phiarray[i])))
        bottom = np.sqrt(np.sum((funcsandder.gausswave_verschoben(np.linspace(x0, xmax, xpoints), 0.5, 0.05,
                                                                  float(times[i])) - phiarray2[i, ::2]) *
                                (funcsandder.gausswave_verschoben(np.linspace(x0, xmax, xpoints), 0.5, 0.05,
                                                                  float(times[i])) - phiarray2[i, ::2])))
        qarray[i] = top / bottom

    ax.plot(times[15:], qarray[15:], color="black")
    ax.xaxis.set_major_locator(MultipleLocator(200))
    ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.0f'))
    # ax.xaxis.set_minor_locator(MultipleLocator(200))
    # ax.xaxis.set_minor_formatter(FormatStrFormatter('% 1.0f'))
    # ax.set_xlim(180, 550)
    ax.set_xlabel("Time")
    ax.set_ylabel(r'$2^p$')

    plt.show()


def plotcompconv(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    xarray = np.linspace(x0, xmax, xpoints)
    xarray2 = np.linspace(x0, xmax, 2 * xpoints - 1)
    timesfd, phiarrayfd, piarrayfd = misc.readfromshortfile("calcfdh.txt", xpoints, lines=linestoread)
    times2fd, phiarray2fd, piarray2fd = misc.readfromshortfile("calcfdh2.txt", 2 * xpoints - 1, lines=linestoread)
    qarrayfd = np.zeros(len(timesfd))
    timesad, phiarrayad, piarrayad = misc.readfromshortfile("calcadh.txt", xpoints, lines=linestoread)
    times2ad, phiarray2ad, piarray2ad = misc.readfromshortfile("calcadh2.txt", 2 * xpoints - 1, lines=linestoread)
    qarrayad = np.zeros(len(timesad))
    timesex, phiarrayex, piarrayex = misc.readfromshortfile("calcexh.txt", xpoints, lines=linestoread)
    times2ex, phiarray2ex, piarray2ex = misc.readfromshortfile("calcexh2.txt", 2 * xpoints - 1, lines=linestoread)
    qarrayex = np.zeros(len(timesex))
    print(timesfd, times2fd, timesad, times2ad, timesex, times2ex)

    for i in range(len(timesfd)):
        topfd = np.sqrt(
            np.sum((funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarrayfd[i]) *
                   (funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarrayfd[i])))
        bottomfd = np.sqrt(
            np.sum((funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarray2fd[i, ::2]) *
                   (funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarray2fd[i, ::2])))
        qarrayfd[i] = topfd / bottomfd

        topad = np.sqrt(
            np.sum((funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarrayad[i]) *
                   (funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarrayad[i])))
        bottomad = np.sqrt(
            np.sum((funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarray2ad[i, ::2]) *
                   (funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarray2ad[i, ::2])))
        qarrayad[i] = topad / bottomad

        topex = np.sqrt(
            np.sum((funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarrayex[i]) *
                   (funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarrayex[i])))
        bottomex = np.sqrt(
            np.sum((funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarray2ex[i, ::2]) *
                   (funcsandder.gausswave_verschoben(xarray, 0.5, 0.05, float(timesfd[i])) - phiarray2ex[i, ::2])))
        qarrayex[i] = topex / bottomex

    fig, ax = plt.subplots()
    ax.plot(timesfd[15:], qarrayfd[15:], color="red", label="Characteristics analysis")
    ax.plot(timesad[15:], qarrayad[15:], color="blue", label="Advection based approach")
    ax.plot(timesex[15:], qarrayex[15:], color="black", label="Extrapolation")
    # ax.xaxis.set_major_locator(MultipleLocator(200))
    # ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.0f'))
    # ax.xaxis.set_minor_locator(MultipleLocator(200))
    # ax.xaxis.set_minor_formatter(FormatStrFormatter('% 1.0f'))
    # ax.set_xlim(180, 550)
    ax.set_xlabel("Time t")
    ax.set_ylabel(r'$2^p$')
    ax.legend()

    plt.show()


def compareconv(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    func_phiint = funcsandder.gausswave
    func_piint = funcsandder.dergaus

    # FD
    xarray, times, phiarray, piarray = solvernew.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                         func_phiint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                         func_piint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                         "FDstencil", "calcfdh.txt", linestoread)
    xarray2, times, phiarray2, piarray2 = solvernew.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                            func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5,
                                                                        0.05),
                                                            func_piint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5,
                                                                       0.05),
                                                            "FDstencil", "calcfdh2.txt", linestoread)
    # adv
    xarray, times, phiarray, piarray = solvernew.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                         func_phiint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                         func_piint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                         "advection", "calcadh.txt", linestoread)
    xarray2, times, phiarray2, piarray2 = solvernew.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                            func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5,
                                                                        0.05),
                                                            func_piint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5,
                                                                       0.05),
                                                            "advection", "calcadh2.txt", linestoread)
    # extra
    xarray, times, phiarray, piarray = solvernew.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                         func_phiint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                         func_piint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                         "extrapolation", "calcexh.txt", linestoread)
    xarray2, times, phiarray2, piarray2 = solvernew.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                            func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5,
                                                                        0.05),
                                                            func_piint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5,
                                                                       0.05),
                                                            "extrapolation", "calcexh2.txt", linestoread)


def conv_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    func_phiint = funcsandder.gausswave
    func_piint = funcsandder.dergaus

    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                      func_phiint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                      func_piint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                      boundaryCondition, "calch.txt",
                                                      linestoread)
    xarray2, times, phiarray2, piarray2 = solver.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                         func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5, 0.05),
                                                         func_piint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5, 0.05),
                                                         boundaryCondition, "calch2.txt", linestoread)



def sineconv(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    linestoread = [timesteps - 1]
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                      funcsandder.s1(np.linspace(x0, xmax, xpoints)),
                                                      funcsandder.Ds1DX(np.linspace(x0, xmax, xpoints)),
                                                      boundaryCondition, fileName,
                                                      linestoread)
    xarray2, times, phiarray2, piarray2 = solver.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                         funcsandder.s1(np.linspace(x0, xmax, 2 * xpoints - 1)),
                                                         funcsandder.Ds1DX(np.linspace(x0, xmax, 2 * xpoints - 1)),
                                                         boundaryCondition, fileName, linestoread)

    fig, ax = plt.subplots(1)
    ax.plot(xarray, phiarray[0, :], color="red")
    ax.plot(xarray, phiarray2[0, ::2])
    plt.show()

    top = np.sqrt(np.sum((funcsandder.s1(xarray) - phiarray[0]) *
                         (funcsandder.s1(xarray) - phiarray[0])))
    bottom = np.sqrt(np.sum((funcsandder.s1(xarray) - phiarray2[0][::2]) *
                            (funcsandder.s1(xarray) - phiarray2[0][::2])))
    print(top / bottom)

    ax.plot(xarray, np.abs((funcsandder.s1(xpoints) - phiarray[0]) / (
            funcsandder.Ds1DX(xpoints) - phiarray2[0][::2])))

    plt.show()

def QoverT_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    func_phiint = funcsandder.gausswave
    func_piint = funcsandder.dergaus

    linestoread = [timesteps - 1]
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                      func_phiint(np.linspace(x0, xmax, xpoints), 0.5, 0.05), func_piint(np.linspace(x0, xmax, xpoints), 0.5, 0.05),
                                                      boundaryCondition, fileName,
                                                      linestoread)
    xarray2, times, phiarray2, piarray2 = solver.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                         func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5, 0.05),
                                                         func_piint(np.linspace(x0, xmax, 2 * xpoints - 1), 0.5, 0.05),
                                                         boundaryCondition, fileName, linestoread)

    fig, ax = plt.subplots(1)
    c = plt.get_cmap('gist_rainbow')
    n = len(linestoread)
    print(linestoread)
    for i in range(n):
        ax[0].plot(xarray, phiarray[i, :], label=format(i * alpha * timesteps / ((len(linestoread) - 1) * (xpoints - 1)), '.2f'), color=c(i / (n - 1)))  #

    ax.legend(loc='upper left')
    plt.show()


def stabtest(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, alpha,
                                                      phiinit, piinit, boundaryCondition, fileName, linestoread)
    print(times)

    fig, ax = plt.subplots(2, 1)
    c = plt.get_cmap('gist_rainbow')
    n = len(linestoread)
    print(linestoread)
    for i in range(n):
        # print("Phi " + str(i) + " :")
        # print(np.mean(phiarray[i]))
        # print("Pi " + str(i) + " :")
        # print(np.mean(piarray[i]))
        ax[0].plot(xarray, phiarray[i, :], label=format(i * alpha * timesteps / ((len(linestoread) - 1) * (xpoints - 1)), '.2f'), color=c(i / (n - 1)))  #
        ax[1].plot(xarray, piarray[i, :], label=format(i * alpha * timesteps / ((len(linestoread) - 1) * (xpoints - 1)), '.2f'), color=c(i / (n - 1)))

    ax[0].set_title('$\phi$')
    ax[1].set_title('$\Pi$')
    ax[0].legend()
    ax[1].legend()
    fig.tight_layout()
    plt.show()


# plot Q over time for every timestep
def convergence(x0, xmax, xpoints, t0, timesteps, alpha,
                phiinit, piinit, boundaryCondition, fileName, linestoread, periods):
    timesteps = periods * (xpoints - 1) / alpha
    linestoread = [0]
    for i in range(1, int(periods / 10 + 1)):
        linestoread.append(int(10 * i * (xpoints - 1)))
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, alpha,
                                                      phiinit, piinit, boundaryCondition, fileName, linestoread)
    print("alpha", alpha, "\ntimes", times)
    fig, ax = plt.subplots(1)
    for i in range(len(linestoread)):
        ax.plot(xarray, phiarray[i, :], label=format(float(times[i]), '.0f'))

    ax.legend(loc='upper left')
    # ax.set_title("alpha" + format(float(alpha), '.2f'))
    plt.show()


def selfconvergence(x0, xmax, xpoints, t0, timesteps, alpha,
                    phiinit, piinit, boundaryCondition, fileName, linestoread):
    func_phiint = funcsandder.s1
    func_piint = np.zeros

    linestoread = [timesteps - 1]
    xarray, times, phiarray, piarray = solver.solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                      func_phiint(np.linspace(x0, xmax, xpoints)), func_piint(xpoints), boundaryCondition, fileName,
                                                      linestoread)
    xarray2, times, phiarray2, piarray2 = solver.solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                         func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1)), func_piint(2 * xpoints - 1),
                                                         boundaryCondition, fileName, linestoread)
    xarray4, times, phiarray4, piarray4 = solver.solving(x0, xmax, 4 * (xpoints - 1) + 1, t0, timesteps + 2, 0.4,
                                                         func_phiint(np.linspace(x0, xmax, 4 * (xpoints - 1) + 1)), func_piint(4 * (xpoints - 1) + 1),
                                                         boundaryCondition, fileName, linestoread)

    fig, ax = plt.subplots(1)
    ax.plot(xarray, (phiarray[0] - phiarray2[0][::2]), label='$\phi\'^{(h)} - \phi\'^{(h/2)}$')
    ax.plot(xarray, 4 * (phiarray2[0][::2] - phiarray4[0][::4]), label='$4 (\phi\'^{(h/2)} - \phi\'^{(h/4)})$', linestyle='dashed')
    # ax.plot(xarray, phiarray[0], label = 'h')
    # ax.plot(xarray, phiarray2[0][::2], label = 'h/2')
    # ax.plot(xarray, phiarray4[0][::4], label = 'h/4')

    ax.legend(loc='upper left')
    plt.show()
