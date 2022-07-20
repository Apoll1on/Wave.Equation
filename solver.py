import numpy as np
import funcsandder

import misc
import os


def boundaryConditions(k, boundaryCondition, delx):
    if boundaryCondition == "periodic":
        k[:, -1] = k[:, 2]
        k[:, 0] = k[:, -3]
    elif boundaryCondition == "extrapolation":
        k[:, -1] = k[:, -2] + (k[:, -2] - k[:, -3]) / delx
        k[:, 0] = k[:, 1] + (k[:, 1] - k[:, 2]) / delx
    elif boundaryCondition == "advection":
        pass
    elif boundaryCondition == "FDstencil":
        pass


def calcK(k, delx, xpoints, boundaryCondition):
    result = np.zeros((2, xpoints + 2), dtype=np.double)
    result[0, 1:-1] = k[1, 1:-1]
    result[1, 1:-1] = (k[0, 2:] - 2 * k[0, 1:-1] + k[0, 0:-2]) / (delx * delx)
    boundaryConditions(result, boundaryCondition, delx)
    return result


# for phi in x direction: p[0], p[xpoints + 1] ghostpoint; p[1], p[xpoints] = x0, xmax; from x0 to xmax (xpoints - 1) xsteps

def solving(xpoints, timesteps, linestoread=[0], fileName="calculateddata.txt", boundaryCondition="periodic", alpha=1):
    t = 0
    delt = alpha / (xpoints - 1)
    tsteps = timesteps
    x0 = 0
    xmax = 1
    delx = (xmax - x0) / (xpoints - 1)
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)

    # File to write Data to
    if os.path.exists(fileName):
        os.remove(fileName)
    f = open(fileName, "a")

    # Set initial values to one of the function s,g. 0 so far.
    u = np.zeros((2, xpoints + 2), dtype=np.double)
    phi = u[0]
    pi = u[1]
    phi[1:-1] = funcsandder.s1(xarray)
    u[0, 1:-1] = funcsandder.s1(xarray)
    # pi[1:-1] = funcsandder.s2(xarray)

    # Ghost Points according to boundary conditions:
    boundaryConditions(u, boundaryCondition, delx)

    tstep = 1
    while tstep < tsteps:
        rhs = calcK(u, delx, xpoints, boundaryCondition)

        k1 = rhs
        k2 = calcK(u + 0.5 * delt * k1, delx, xpoints, boundaryCondition)
        k3 = calcK(u + 0.5 * delt * k2, delx, xpoints, boundaryCondition)
        k4 = calcK(u + delt * k3, delx, xpoints, boundaryCondition)

        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

        boundaryConditions(u, boundaryCondition, delx)

        misc.savedata(f, (t, u[0], u[1]))

        # Advance time
        tstep = tstep + 1
        t = t + delt

    f.close()
    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)
    return (xarray, times, phiarray, piarray)
