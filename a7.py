import numpy as np
import funcsandder
from matplotlib import pyplot as plt
import misc
import os


def calcRHS(phi, pi, delx, xpoints):
    result = np.zeros((2, xpoints + 2))
    result[0, 1:-1] = pi[1:-1]
    result[1, 1:-1] = (phi[2:] - 2 * phi[1:-1] + phi[0:-2]) / (delx * delx)
    return result


# for phi in x direction: p[0], p[xpoints + 1] ghostpoint; p[1], p[xpoints] = x0, xmax; from x0 to xmax (xpoints - 1) xsteps

def solving(xpoints, timesteps, linestoread=[0], fileName="calculateddata.txt", boundaryCondition="periodic",
            BCimpl=None, alpha=0.1):
    t = 0
    delt = alpha / (xpoints - 1)
    tsteps = timesteps
    x0 = 0
    xmax = 1
    delx = (xmax - x0) / (xpoints - 1)
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)

    # File to write Data to
    os.remove(fileName)
    f = open(fileName, "a")

    # Set initial values to one of the function s,g. 0 so far.
    pi = np.zeros(xpoints + 2, dtype=np.double)
    phi = np.zeros(xpoints + 2, dtype=np.double)
    phi[1:-1] = funcsandder.s1(xarray)
    pi[1:-1] = np.zeros(xpoints)

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


    tstep = 1
    while tstep < tsteps:

        rhs = np.zeros((2, xpoints + 2), dtype=np.double)
        rhs = calcRHS(phi, pi, delx, xpoints)

        # debugging
        # fig, ax = plt.subplots(2)
        # ax[0].plot(xarray, rhs[1, 1:-1])
        # ax[0].plot(xarray, funcsandder.D2s1D2X(xarray))
        # ax[1].plot(xarray, funcsandder.D2s1D2X(xarray) - rhs[1, 1:-1])
        # plt.title("tstep= "+str(tstep))
        # plt.show()

        # xstep = 1  # set to one because the zeroth entry is the ghost point
        # while xstep <= xpoints:
        #     # Berechnung der Ki:
        #     kphi1 = rhs[0, xstep] + 0.
        #     kpi1 = rhs[1, xstep] + 0.
        #     kphi2 = rhs[0, xstep] + delt * kpi1 / 2
        #     kpi2 = rhs[1, xstep] + 0.  # TODO
        #     kphi3 = rhs[0, xstep] + delt * kpi2 / 2
        #     kpi3 = rhs[1, xstep] + 0.
        #     kphi4 = rhs[0, xstep] + delt * kpi3
        #     kpi4 = rhs[1, xstep] + 0.
        #
        #     pi[xstep] = pi[xstep] + delt * (kpi1 / 6 + kpi2 / 3 + kpi3 / 3 + kpi4 / 6)
        #     phi[xstep] = phi[xstep] + delt * (kphi1 / 6 + kphi2 / 3 + kphi3 / 3 + kphi4 / 6)
        #     xstep = xstep + 1

        kphi1 = rhs[0]
        kpi1 = rhs[1]


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

        misc.savedata(f, (t, phi, pi))

        # Advance time
        tstep = tstep + 1
        t = t + delt

    f.close()
    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)
    return (xarray, times, phiarray, piarray)
