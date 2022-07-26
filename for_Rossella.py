import numpy as np
import misc
import os
import funcsandder
from matplotlib import pyplot as plt


#### solver file:

def boundaryConditions(u, boundaryCondition, delx, alpha=1, phi_old_left=0, phi_old_right=0):
    if boundaryCondition == "periodic":
        u[:, -1] = u[:, 2]
        u[:, 0] = u[:, -3]
    elif boundaryCondition == "extrapolation":
        u[:, -1] = u[:, -2] + (u[:, -2] - u[:, -3])
        u[:, 0] = u[:, 1] + (u[:, 1] - u[:, 2])
    elif boundaryCondition == "advection":
        u[0, 0] = u[0, 2] - 2 * delx * u[1, 1]
        u[0, -1] = u[0, -3] - 2 * delx * u[1, -2]
        u[1, -1] = u[1, -2] + (u[1, -2] - u[1, -3])
        u[1, 0] = u[1, 1] + (u[1, 1] - u[1, 2])
    elif boundaryCondition == "FDstencil":
        u[0, 0] = (alpha * (2 * u[0, 1] - u[0, 2]) + u[0, 2] + (2 * phi_old_left - 2 * u[0, 1]) / alpha) / (alpha + 1)
        u[0, -1] = (alpha * (2 * u[0, -2] - u[0, -3]) + u[0, -3] + (2 * phi_old_right - 2 * u[0, -2]) / alpha) / (alpha + 1)
        u[1, -1] = u[1, -2] + (u[1, -2] - u[1, -3])
        u[1, 0] = u[1, 1] + (u[1, 1] - u[1, 2])


def calcRHS(u, delx, xpoints, boundaryCondition, alpha=1, phi_old_left=0, phi_old_right=0):
    result = np.zeros((2, xpoints + 2), dtype=np.double)
    result[0, 1:-1] = u[1, 1:-1]
    result[1, 1:-1] = (u[0, 2:] - 2 * u[0, 1:-1] + u[0, 0:-2]) / (delx * delx)
    boundaryConditions(result, boundaryCondition, delx, alpha, phi_old_left, phi_old_right)
    return result


# for phi in x direction: p[0], p[xpoints + 1] ghostpoint; p[1], p[xpoints] = x0, xmax; from x0 to xmax (xpoints - 1) xsteps

def solving(x0, xmax, xpoints, t0, timesteps, alpha,
            phiinit, piinit, boundaryCondition, fileName, linestoread):
    t = t0
    delt = alpha / (xpoints - 1)
    delx = (xmax - x0) / (xpoints - 1)
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)

    # File to write Data to
    if os.path.exists(fileName):
        os.remove(fileName)
    f = open(fileName, "a")

    # Set initial values to one of the function s,g. 0 so far.
    u = np.zeros((2, xpoints + 2), dtype=np.double)
    u[0, 1:-1] = phiinit
    u[1, 1:-1] = piinit

    if boundaryCondition != "FDstencil":
        # Ghost Points according to boundary conditions:
        boundaryConditions(u, boundaryCondition, delx)

        misc.savedata(f, t, u[0], u[1])
        tstep = 1
        while tstep < timesteps:
            k1 = calcRHS(u, delx, xpoints, boundaryCondition)
            k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints, boundaryCondition)
            k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, boundaryCondition)
            k4 = calcRHS(u + delt * k3, delx, xpoints, boundaryCondition)

            u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

            boundaryConditions(u, boundaryCondition, delx)

            # Advance time
            tstep = tstep + 1
            t = t + delt
            misc.savedata(f, t, u[0], u[1])

    else:
        # Ghost Points according to boundary conditions:
        boundaryConditions(u, "advection", delx)

        misc.savedata(f, t, u[0], u[1])
        phi_old = np.zeros((2, 2))
        phi_old[1, 0] = u[0, 1]
        phi_old[1, 1] = u[0, -2]

        k1 = calcRHS(u, delx, xpoints, "advection")
        k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints, "advection")
        k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, "advection")
        k4 = calcRHS(u + delt * k3, delx, xpoints, "advection")

        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
        phi_old[0, 0] = u[0, 1]
        phi_old[0, 1] = u[0, -2]

        boundaryConditions(u, "advection", delx)

        # Advance time
        t = t + delt
        misc.savedata(f, t, u[0], u[1])

        tstep = 2
        while tstep < timesteps:
            k1 = calcRHS(u, delx, xpoints, boundaryCondition, alpha, phi_old[1, 0], phi_old[1, 1])
            k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints, boundaryCondition, alpha, phi_old[1, 0], phi_old[1, 1])
            k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, boundaryCondition, alpha, phi_old[1, 0], phi_old[1, 1])
            k4 = calcRHS(u + delt * k3, delx, xpoints, boundaryCondition, alpha, phi_old[1, 0], phi_old[1, 1])

            u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

            boundaryConditions(u, boundaryCondition, delx, alpha, phi_old[1, 0], phi_old[1, 1])
            phi_old[1, :] = phi_old[0, :]
            phi_old[0, 0] = u[0, 1]
            phi_old[0, 1] = u[0, -2]

            # Advance time
            t = t + delt
            misc.savedata(f, t, u[0], u[1])
            tstep = tstep + 1


    f.close()
    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)
    return xarray, times, phiarray, piarray



#### a8 file:

def stabtest(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread):
    xarray, times, phiarray, piarray = solving(x0, xmax, xpoints, t0, timesteps + 2, alpha,
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
        ax[0].plot(xarray, phiarray[i, :], label='t: ' + format(float(times[i]), '.4f'), color=c(i / (n - 1)))
        ax[1].plot(xarray, piarray[i, :], label='t: ' + format(float(times[i]), '.4f'), color=c(i / (n - 1)))

    ax[0].set_title('Phi')
    ax[1].set_title('Pi')
    ax[0].legend()
    ax[1].legend()
    plt.show()


def convergence(x0, xmax, xpoints, t0, timesteps, alpha,
                phiinit, piinit, boundaryCondition, fileName, linestoread, periods):
    timesteps = periods * (xpoints - 1) / alpha
    linestoread = [0]
    for i in range(1, periods + 1):
        linestoread.append(int(i * (xpoints - 1)))
    xarray, times, phiarray, piarray = solving(x0, xmax, xpoints, t0, timesteps + 2, alpha,
                                                      phiinit, piinit, boundaryCondition, fileName, linestoread)
    print("alpha", alpha, "\ntimes", times)
    fig, ax = plt.subplots(1)
    for i in range(len(linestoread)):
        ax.plot(xarray, phiarray[i, :], label=format(float(times[i]), '.4f'))

    ax.legend()
    ax.set_title("alpha" + format(float(alpha), '.2f'))
    plt.show()


def selfconvergence(x0, xmax, xpoints, t0, timesteps, alpha,
                    phiinit, piinit, boundaryCondition, fileName, linestoread):
    func_phiint = funcsandder.s1
    func_piint = np.zeros

    linestoread = [timesteps - 1]
    xarray, times, phiarray, piarray = solving(x0, xmax, xpoints, t0, timesteps + 2, 0.1,
                                                      func_phiint(np.linspace(x0, xmax, xpoints)), func_piint(xpoints), boundaryCondition, fileName,
                                                      linestoread)
    xarray2, times, phiarray2, piarray2 = solving(x0, xmax, 2 * xpoints - 1, t0, timesteps + 2, 0.2,
                                                         func_phiint(np.linspace(x0, xmax, 2 * xpoints - 1)), func_piint(2 * xpoints - 1),
                                                         boundaryCondition, fileName, linestoread)
    xarray4, times, phiarray4, piarray4 = solving(x0, xmax, 4 * (xpoints - 1) + 1, t0, timesteps + 2, 0.4,
                                                         func_phiint(np.linspace(x0, xmax, 4 * (xpoints - 1) + 1)), func_piint(4 * (xpoints - 1) + 1),
                                                         boundaryCondition, fileName, linestoread)

    fig, ax = plt.subplots(1)
    ax.plot(xarray, (phiarray[0] - phiarray2[0][::2]), label='h - h/2')  #
    ax.plot(xarray, 4 * (phiarray2[0][::2] - phiarray4[0][::4]), label='h/2 - h/4')  #
    # ax.plot(xarray, phiarray[0], label = 'h')
    # ax.plot(xarray, phiarray2[0][::2], label = 'h/2')
    # ax.plot(xarray, phiarray4[0][::4], label = 'h/4')

    ax.legend()
    plt.show()




#### main file:


# xsteps
x0 = 0
xmax = 1
xpoints = 1001

# timesteps
timesteps = 1500  # number of timesteps
tcount = 6  # Number of lines to read out/plot
t0 = 0  # starting time

# alpha
alpha = 1

# Initial conditions
xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)
phiinit = funcsandder.gausswave(xarray, 0.5, 0.05)
piinit = funcsandder.dergaus(xarray, 0.5, 0.05)
# phiinit = funcsandder.s1(xarray)
# piinit = np.zeros(xpoints)


# phiinit = np.zeros(xpoints)
# phiinit[int(-x0 * (xpoints - 1) / (xmax - x0)): int((1-x0) * (xpoints - 1) / (xmax - x0) + 1)] = funcsandder.s1(np.linspace(0, 1, int((xpoints - 1) / (xmax - x0) + 1)))

# Boundary condition
boundaryCondition = "advection"

# read data
fileName = "calculateddata.txt"
linestoread = [0]
for i in range(1, tcount + 1):
    linestoread.append(i * timesteps / tcount)

stabtest(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)
# convergence(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread, 15) # Last, additional parameter is for number of periods
# selfconvergence(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)