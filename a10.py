import matplotlib.pyplot as plt
import numpy as np
import funcsandder
from numba import jit
import misc
import os
import time


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


def calcRHS(k, delx, xpoints, xarray, boundaryCondition, pot):
    result = np.zeros((2, xpoints + 2), dtype=np.double)
    result[0, 1:-1] = k[1, 1:-1]
    result[1, 1:-1] = (k[0, 2:] - 2 * k[0, 1:-1] + k[0, 0:-2]) / (delx * delx) - pot  # subtracting the potential
    boundaryConditions(result, boundaryCondition, delx)
    return result


# for phi in x direction: p[0], p[xpoints + 1] ghostpoint; p[1], p[xpoints] = x0, xmax; from x0 to xmax (xpoints - 1) xsteps

def solving(x0,xmax,xpoints,t0,timesteps,alpha,
                                                   phiinit,piinit,boundaryCondition,fileName,linestoread):
    t=t0
    delt = alpha / (xpoints - 1)
    delx = (xmax - x0) / (xpoints - 1)
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)

    pot = PTpotential(xarray)

    # File to write Data to
    if os.path.exists(fileName):
        os.remove(fileName)
    f = open(fileName, "a")

    # Set initial values to one of the function s,g. 0 so far.
    u = np.zeros((2, xpoints + 2), dtype=np.double)
    u[0, 1:-1] = phiinit
    u[1, 1:-1] = piinit

    # Ghost Points according to boundary conditions:
    boundaryConditions(u, boundaryCondition, delx)

    misc.savedata(f, (t, u[0], u[1]))
    tstep = 1
    while tstep < timesteps:
        #start_time = time.time()
        k1 = calcRHS(u, delx, xpoints, xarray, boundaryCondition, pot)
        k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints, xarray, boundaryCondition, pot)
        k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, xarray, boundaryCondition, pot)
        k4 = calcRHS(u + delt * k3, delx, xpoints, xarray, boundaryCondition, pot)
        #print("--- %s seconds ---" % (time.time() - start_time))
        #start_time = time.time()
        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
        #print("--- %s seconds ---" % (time.time() - start_time))
        #start_time = time.time()
        boundaryConditions(u, boundaryCondition, delx)
        #print("--- %s seconds ---" % (time.time() - start_time))
        # Advance time
        tstep = tstep + 1
        t = t + delt

        start_time = time.time()
        misc.savedata(f, (t, u[0], u[1]))
        print("--- %s seconds ---" % (time.time() - start_time))


    f.close()
    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)
    return (xarray, times, phiarray, piarray)
