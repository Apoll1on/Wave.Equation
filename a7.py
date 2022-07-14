import numpy as np
import funcsandder
from matplotlib import pyplot as plt

import misc


def calcRHS(phi, pi, delx, xsteps):
    result = np.zeros((2, xsteps + 2), dtype=np.double)
    result[1, 1:-1] = (phi[2:] - 2 * phi[1:-1] + phi[:-2]) / (delx * delx)
    result[0, 1:-1] = pi[1:-1] + 0.

    # for i in range(1, xsteps + 1):
    #     result[1, i] = (phi[i + 1] - 2 * phi[i] + phi[i - 1]) / (delx * delx)
    #     result[0, i] = pi[i]
    return result

# for p: in x direction: p[0], p[xsteps + 1] ghostpoint; p[1], p[xsteps] = x0, xmax; from x0 to xmax (xsteps - 1) steps

def solving(xsteps, timesteps, linestoread=[0], fileName="calculateddata.txt", boundaryCondition="periodic",
            BCimpl=None, alpha=1):
    t = 0
    delt = 0.001 #alpha / (xsteps - 1)
    tsteps = timesteps
    x0 = 0
    xmax = 1
    delx = (xmax - x0) / (xsteps - 1)
    xarray = np.array(np.linspace(x0, xmax, xsteps), dtype=np.double)

    # piarray = np.zeros((tsteps, xsteps + 2), dtype=np.double)
    # phiarray = np.zeros((tsteps, xsteps + 2), dtype=np.double)

    # File to write Data to
    f = open(fileName, "a")

    # Set initial values to one of the function s,g. 0 so far.
    pi = np.zeros(xsteps + 2, dtype=np.double)
    phi = np.zeros(xsteps + 2, dtype=np.double)
    phi[1:-1] = funcsandder.f3(xarray)
    pi[1:-1] = np.zeros(xsteps)

    # Ghost Points according to boundary conditions:
    if boundaryCondition == "periodic":
        phi[-1] = phi[2]
        phi[0] = phi[-3]

        pi[-1] = pi[2]
        pi[0] = pi[-3]
    elif boundaryCondition == "open":
        if BCimpl == "extrapolation":
            pass
        elif BCimpl == "advection":
            pass
        elif BCimpl == "FDstencil":
            pass

    # For Saving Data
    # piarray[0, :] = pi[:] + 0.
    # phiarray[0, :] = phi[:] + 0.

    tstep = 1
    while tstep < tsteps:

        rhs = np.zeros((2, xsteps + 2), dtype=np.double)
        rhs = calcRHS(phi, pi, delx, xsteps)

        # debugging
        # fig, ax = plt.subplots(2)
        # ax[0].plot(xarray, rhs[1, 1:-1])
        # ax[0].plot(xarray, funcsandder.D2s1D2X(xarray))
        # ax[1].plot(xarray, funcsandder.D2s1D2X(xarray) - rhs[1, 1:-1])
        # plt.title("tstep= "+str(tstep))
        # plt.show()

        xstep = 1  # set to one because the zeroth entry is the ghost point
        while xstep <= xsteps:
            # Berechnung der Ki:
            kphi1 = rhs[0, xstep] + 0.
            kpi1 = rhs[1, xstep] + 0.
            kphi2 = rhs[0, xstep] + delt * kpi1 / 2
            kpi2 = rhs[1, xstep] + 0.  # TODO
            kphi3 = rhs[0, xstep] + delt * kpi2 / 2
            kpi3 = rhs[1, xstep] + 0.
            kphi4 = rhs[0, xstep] + delt * kpi3
            kpi4 = rhs[1, xstep] + 0.

            pi[xstep] = pi[xstep] + delt * (kpi1 / 6 + kpi2 / 3 + kpi3 / 3 + kpi4 / 6)
            phi[xstep] = phi[xstep] + delt * (kphi1 / 6 + kphi2 / 3 + kphi3 / 3 + kphi4 / 6)
            xstep = xstep + 1

        # Set new BC:
        if boundaryCondition == "periodic":
            phi[1] = phi[-2] + 0.
            pi[1] = pi[-2] + 0.
            pass
        elif boundaryCondition == "open":
            if BCimpl == "extrapolation":
                pass
            elif BCimpl == "advection":
                pass
            elif BCimpl == "FDstencil":
                pass

        # Ghost Points
        if boundaryCondition == "periodic":
            phi[-1] = phi[2]
            phi[0] = phi[-3]

            pi[-1] = pi[2]
            pi[0] = pi[-3]
        elif boundaryCondition == "open":
            if BCimpl == "extrapolation":
                phi[0] = phi[1] + (phi[1] - phi[2]) / delx
            elif BCimpl == "advection":
                pass
            elif BCimpl == "FDstencil":
                pass

        # Save data
        # piarray[tstep, :] = pi[:] + 0.
        # phiarray[tstep, :] = phi[:] + 0.
        misc.savedata(f, (t, phi, pi))

        # Advance time
        tstep = tstep + 1
        t = t + delt

    # plotting

    f.close()
    times, phiarray, piarray = misc.readdata("calculateddata.txt", xsteps, lines=linestoread)
    return (xarray, times, phiarray, piarray)
