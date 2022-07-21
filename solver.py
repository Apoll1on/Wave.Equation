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


def calcRHS(k, delx, xpoints, boundaryCondition):
    result = np.zeros((2, xpoints + 2), dtype=np.double)
    result[0, 1:-1] = k[1, 1:-1]
    result[1, 1:-1] = (k[0, 2:] - 2 * k[0, 1:-1] + k[0, 0:-2]) / (delx * delx)
    boundaryConditions(result, boundaryCondition, delx)
    return result


# for phi in x direction: p[0], p[xpoints + 1] ghostpoint; p[1], p[xpoints] = x0, xmax; from x0 to xmax (xpoints - 1) xsteps

def solving(xpoints, tsteps, alpha, boundaryCondition, linestoread=[0], fileName="calculateddata.txt"):
    t = 0
    delt = alpha / (xpoints - 1)
    x0 = -4
    xmax = 4
    delx = (xmax - x0) / (xpoints - 1)
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)

    # File to write Data to
    if os.path.exists(fileName):
        os.remove(fileName)
    f = open(fileName, "a")

    # Set initial values to one of the function s,g
    u = np.zeros((2, xpoints + 2), dtype=np.double)
    u[0, 1:-1] = funcsandder.gausswave(xarray, 0., 0.05)  # funcsandder.s1(xarray)
    u[1, 1:-1] = -funcsandder.dergaus(xarray, 0., 0.05)

    # Ghost Points according to boundary conditions:
    boundaryConditions(u, boundaryCondition, delx)

    tstep = 1
    while tstep < tsteps:
        rhs = calcRHS(u, delx, xpoints, boundaryCondition)

        k1 = rhs
        k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints, boundaryCondition)
        k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, boundaryCondition)
        k4 = calcRHS(u + delt * k3, delx, xpoints, boundaryCondition)

        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

        boundaryConditions(u, boundaryCondition, delx)

        misc.savedata(f, (t, u[0], u[1]))

        # Advance time
        tstep = tstep + 1
        t = t + delt

    f.close()
    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)
    return (xarray, times, phiarray, piarray)
