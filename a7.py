import numpy as np
import funcsandder
from matplotlib import pyplot as plt


def calcRHS(phi, pi, delx, xsteps):
    result = np.zeros((2, xsteps + 2), dtype=np.double)
    result[1, 1:-1] = (phi[2:] - 2 * phi[1:-1] + phi[:-2]) / (delx * delx)
    result[0, 1:-1] = pi[1:-1] + 0.

    # Boundary Conditions by copying:
    result[:, 0] = result[:, 1] + 0.
    result[:, -1] = result[:, -2] + 0.

    return result


def solving(xsteps=21, fileName="calculateddata"):
    t = 0
    tsteps = 100
    tmax = 1
    delt = tmax / tsteps
    x0 = 0
    xmax = 1
    delx = (xmax - x0) / (xsteps - 1)
    print(delt, delx)

    piarray = np.zeros((tsteps, xsteps + 2), dtype=np.double)
    phiarray = np.zeros((tsteps, xsteps + 2), dtype=np.double)

    f = open(fileName, "a")
    f.write("\n\nNew Run")

    # Set initial values to one of the function s,g. 0 so far.
    pi = np.zeros(xsteps + 2, dtype=np.double)
    phi = np.zeros(xsteps + 2, dtype=np.double)
    print(np.shape(np.arange(0, 1 + delx, delx, dtype=np.double)))
    print(np.shape(funcsandder.s1(np.arange(0, 1 + delx, delx, dtype=np.double))))
    pi[1:-1] = funcsandder.s1(np.arange(0, 1 + delx, delx, dtype=np.double))
    phi[1:-1] = funcsandder.s1(np.arange(0, 1 + delx, delx, dtype=np.double))
    piarray[0, :] = pi[:] + 0.
    phiarray[0, :] = phi[:] + 0.

    while t < tmax:
        tstep = 0
        rhs = np.zeros((2, xsteps + 2), dtype=np.double)
        rhs = calcRHS(phi, pi, delx, xsteps)
        if tstep == 0:
            print(rhs)

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

            pi[xstep + 1] = pi[xstep + 1] + delt * (kpi1 / 6 + kpi2 / 3 + kpi3 / 3 + kpi4 / 6)
            phi[xstep + 1] = phi[xstep + 1] + delt * (kphi1 / 6 + kphi2 / 3 + kphi3 / 3 + kphi4 / 6)
            xstep = xstep + 1

        t = t + delt
        tstep = tstep + 1
        piarray[tstep, :] = pi[:] + 0.
        phiarray[tstep, :] = phi[:] + 0.

    xarray = np.array(np.arange(x0, xmax + delx, delx), dtype=np.double)
    tarray = np.array(np.arange(0, tmax + delt, delt), dtype=np.double)

    fig1, ax1 = plt.subplots(3, 3)
    ax1[0, 0].plot(xarray, piarray[0, 1:-1])
    ax1[0, 1].plot(xarray, piarray[12, 1:-1])
    ax1[0, 2].plot(xarray, piarray[25, 1:-1])
    ax1[1, 0].plot(xarray, piarray[37, 1:-1])
    ax1[1, 1].plot(xarray, piarray[50, 1:-1])
    ax1[1, 2].plot(xarray, piarray[62, 1:-1])
    ax1[2, 0].plot(xarray, piarray[75, 1:-1])
    ax1[2, 1].plot(xarray, piarray[87, 1:-1])
    ax1[2, 2].plot(xarray, piarray[99, 1:-1])

    fig2, ax2 = plt.subplots(3, 3)
    ax2[0, 0].plot(xarray, phiarray[0, 1:-1])
    ax2[0, 1].plot(xarray, phiarray[12, 1:-1])
    ax2[0, 2].plot(xarray, phiarray[25, 1:-1])
    ax2[1, 0].plot(xarray, phiarray[37, 1:-1])
    ax2[1, 1].plot(xarray, phiarray[50, 1:-1])
    ax2[1, 2].plot(xarray, phiarray[62, 1:-1])
    ax2[2, 0].plot(xarray, phiarray[75, 1:-1])
    ax2[2, 1].plot(xarray, phiarray[87, 1:-1])
    ax2[2, 2].plot(xarray, phiarray[99, 1:-1])

    plt.show()
