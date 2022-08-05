import time

import numpy as np
import misc
import os


def boundaryConditions(k, boundaryCondition, delx):
    """Applying different boundary conditions."""
    # We only used extrapolation here so far because the boundaries do not really matter anyway
    if boundaryCondition == "periodic":
        k[:, -1] = k[:, 2]
        k[:, 0] = k[:, -3]
    elif boundaryCondition == "extrapolation":
        k[:, -1] = k[:, -2] + (k[:, -2] - k[:, -3])
        k[:, 0] = k[:, 1] + (k[:, 1] - k[:, 2])
    elif boundaryCondition == "advection":
        pass
    elif boundaryCondition == "FDstencil":
        k[0, 0] = k[0, 2] - 2 * delx * k[1, 1]
        k[0, -1] = k[0, -3] - 2 * delx * k[1, -2]
        k[1, -1] = k[1, -2] + (k[1, -2] - k[1, -3])
        k[1, 0] = k[1, 1] + (k[1, 1] - k[1, 2])


def PTpotential(xarray):
    """Potential calculated for every point"""
    return 0.15 / ((np.cosh(0.18 * (xarray) + 0.43)) * (np.cosh(0.18 * (xarray) + 0.43)))


def calcRHS(k, delx, xpoints, xarray, boundaryCondition, pot):
    """calculating the RHS of equation du/dt=F(u) for RK4 step"""
    result = np.zeros((2, xpoints + 2), dtype=np.double)
    result[0, 1:-1] = k[1, 1:-1]

    # This is the part where the potential comes into play and I hope I implemented it correctly
    result[1, 1:-1] = (k[0, 2:] - 2 * k[0, 1:-1] + k[0, 0:-2]) / (delx * delx) - pot * k[0, 1:-1] # subtracting the potential
    boundaryConditions(result, boundaryCondition, delx)
    return result


def solving(x0, xmax, xpoints, t0, timesteps, alpha,
            phiinit, piinit, boundaryCondition, fileName, linestoread):
    """Here is where everything is put together. First initial data is set up and then the time stepping loop begins."""
    starttime=time.time()
    t = t0
    delt = alpha / (xpoints - 1)
    delx = (xmax - x0) / (xpoints - 1)
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)

    # calc. potential to reuse later
    pot = PTpotential(xarray)


    # save data in array to plot it later:
    phiarray = np.zeros((len(linestoread), xpoints + 2), dtype=np.double)
    piarray = np.zeros((len(linestoread), xpoints + 2), dtype=np.double)
    times = []
    index = 0

    # Set initial values to one of the function s,g. 0 so far.
    # u[0] is phi, where u[0,0] and u[0,xpoints+1] are the ghostpoints.
    # u[1] is pi, ghostpoints accordingly
    u = np.zeros((2, xpoints + 2), dtype=np.double)
    u[0, 1:-1] = phiinit
    u[1, 1:-1] = piinit

    # this is just for printing data later on.
    # not all timesteps needed for plotting
    if 0 in linestoread:
        times.append(t)
        phiarray[index, :] = u[0, :]
        piarray[index, :] = u[1, :]
        index = index + 1

    # Ghost Points according to boundary conditions:
    boundaryConditions(u, boundaryCondition, delx)

    # time stepping loop
    tstep = 1
    while tstep < timesteps:

        # calculate the k
        k1 = calcRHS(u, delx, xpoints, xarray, boundaryCondition, pot)
        k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints, xarray, boundaryCondition, pot)
        k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, xarray, boundaryCondition, pot)
        k4 = calcRHS(u + delt * k3, delx, xpoints, xarray, boundaryCondition, pot)

        # new values for phi and pi
        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

        boundaryConditions(u, boundaryCondition, delx)


        # Advance time
        t = t + delt

        # again for saving/plotting data
        if tstep in linestoread:
            times.append(t)
            phiarray[index, :] = u[0, :]
            piarray[index, :] = u[1, :]
            index = index + 1
            print(index, " of ", len(linestoread))

        tstep = tstep + 1

    print(time.time()-starttime)
    if os.path.exists(fileName):
        os.remove(fileName)
    f = open(fileName, "a")

    for i in range(len(linestoread)):
        misc.savedata(f, times[i], phiarray[i], piarray[i])

    #when returning we can leave out the ghost points
    return (xarray, times, phiarray[:, 1:-1], piarray[:, 1:-1])
