import matplotlib.pyplot as plt
import numpy as np
import funcsandder

import misc
import os


def boundaryConditions(k, boundaryCondition, delx):
    if boundaryCondition == "periodic":
        k[:, -1] = k[:, 2]
        k[:, 0] = k[:, -3]
    elif boundaryCondition == "extrapolation":
        k[:, -1] = k[:, -2] + (k[:, -2] - k[:, -3])
        k[:, 0] = k[:, 1] + (k[:, 1] - k[:, 2])
    elif boundaryCondition == "advection":
        pass
    elif boundaryCondition == "FDstencil":
        pass


def PTpotential(xarray):
    return 0.15 / ((np.cosh(0.18 * (xarray) + 0.43)) * (np.cosh(0.18 * (xarray) + 0.43)))


def calcRHS(k, delx, xpoints, xarray, boundaryCondition):
    result = np.zeros((2, xpoints + 2), dtype=np.double)
    pot = PTpotential(xarray)
    result[0, 1:-1] = k[1, 1:-1]
    result[1, 1:-1] = (k[0, 2:] - 2 * k[0, 1:-1] + k[0, 0:-2]) / (delx * delx) - pot  # subtracting the potential
    boundaryConditions(result, boundaryCondition, delx)
    return result


# for phi in x direction: p[0], p[xpoints + 1] ghostpoint; p[1], p[xpoints] = x0, xmax; from x0 to xmax (xpoints - 1) xsteps

def solving(xpoints, tsteps, alpha, boundaryCondition, linestoread=[0], fileName="calculateddata.txt"):
    t = 0
    delt = alpha / (xpoints - 1)
    x0 = -25
    xmax = 15
    delx = (xmax - x0) / (xpoints - 1)
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)

    fig, ax = plt.subplots(1)
    ax.plot(xarray, PTpotential(xarray))
    plt.show()

    # File to write Data to
    if os.path.exists(fileName):
        os.remove(fileName)
    f = open(fileName, "a")

    # Set initial values to one of the function s,g. 0 so far.
    u = np.zeros((2, xpoints + 2), dtype=np.double)
    u[0, 1:-1] = funcsandder.gausswave(xarray, -20., 0.05)  # funcsandder.s1(xarray)
    u[1, 1:-1] = funcsandder.dergaus(xarray, -20., 0.05)

    # Ghost Points according to boundary conditions:
    boundaryConditions(u, boundaryCondition, delx)

    misc.savedata(f, (t, u[0], u[1]))
    tstep = 1
    while tstep < tsteps:
        k1 = calcRHS(u, delx, xpoints, xarray,boundaryCondition)
        k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints,xarray, boundaryCondition)
        k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, xarray,boundaryCondition)
        k4 = calcRHS(u + delt * k3, delx, xpoints, xarray,boundaryCondition)

        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

        boundaryConditions(u, boundaryCondition, delx)

        # Advance time
        tstep = tstep + 1
        t = t + delt
        misc.savedata(f, (t, u[0], u[1]))

    f.close()
    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)
    return (xarray, times, phiarray, piarray)
