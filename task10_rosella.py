import numpy as np
from matplotlib import pyplot as plt
import time
import funcsandder


def boundaryConditions(k, boundaryCondition, delx):
    """Applying different boundary conditions"""
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
    """Potential calculated for every point"""
    return 0.15 / ((np.cosh(0.18 * (xarray) + 0.43)) * (np.cosh(0.18 * (xarray) + 0.43)))



def calcRHS(k, delx, xpoints, xarray, boundaryCondition, pot):
    """calculating the RHS of equation du/dt=F(u) for RK4 step"""
    result = np.zeros((2, xpoints + 2), dtype=np.double)
    result[0, 1:-1] = k[1, 1:-1]

    # This is the part where the potential comes into play and I hope I implented it correctly
    result[1, 1:-1] = (k[0, 2:] - 2 * k[0, 1:-1] + k[0, 0:-2]) / (delx * delx) - pot  # subtracting the potential
    boundaryConditions(result, boundaryCondition, delx)
    return result


def solving(x0, xmax, xpoints, t0, timesteps, alpha,
            phiinit, piinit, boundaryCondition, fileName, linestoread):
    """Here is where the calculation is put together. First initial data is set up and then the time stepping loop begins."""
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

    #when returning we can leave out the ghost points
    return (xarray, times, phiarray[:, 1:-1], piarray[:, 1:-1])





def calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread):
    """Here we calculate for some given initial data and then also plot it."""

    #initial data
    # xsteps
    x0 = -50
    xmax = 200
    xpoints = 5001

    # timesteps
    periods = 1500#just to simplify plotting
    timesteps = periods * 1000  # number of timesteps
    t0 = 0  # starting time

    # alpha
    alpha = 1

    # Initial conditions
    xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)
    phiinit = funcsandder.gausswave(xarray, -40., 0.1)  # phi. just a gauss wave pulse with mu=-40 and sigma=0.1
    piinit = funcsandder.dergaus(xarray, -40, 0.1)  # pi. timelike derivative of gauss pulse above

    # Boundary condition
    boundaryCondition = "extrapolation"



    #calculation
    xarray, times, phiarray, piarray = solving(x0,xmax,xpoints,t0,timesteps+2,alpha,
                                                   phiinit,piinit,boundaryCondition,fileName,linestoread)

    #plotting phi and pi in subplot 1 and 2.
    #Only printing every 100*1000=100000 timesteps, therefore in the loop it says 100*i
    fig, ax = plt.subplots(2, 1)
    c = plt.get_cmap('gist_rainbow')
    n = len(linestoread)
    for i in range(int(n/100)):
        # print("Phi " + str(i) + " :")
        # print(np.mean(phiarray[i]))
        # print("Pi " + str(i) + " :")
        # print(np.mean(piarray[i]))
        ax[0].plot(xarray, phiarray[100*i, :]+100*i*a10ram.PTpotential(xarray), label=format(float(times[100*i]), '.4f'), color = c(100*i/(n-1)))
        ax[1].plot(xarray, piarray[100*i, :], label=format(float(times[100*i]), '.4f'), color = c(100*i/(n-1)))

    ax[0].set_title('phi')
    ax[1].set_title('pi')
    ax[0].legend()
    ax[1].legend()
    plt.show()
