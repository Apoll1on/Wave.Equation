import numpy as np
import misc
import os


def boundaryConditions(u, boundaryCondition, delx, alpha=1, phi_old_left=0, phi_old_right=0):
    if boundaryCondition == "periodic":
        u[:, -1] = u[:, 2]
        u[:, 0] = u[:, -3]
    elif boundaryCondition == "extrapolation":
        u[:, -1] = 2 * u[:, -2] - u[:, -3]
        u[:, 0] = 2 * u[:, 1]- u[:, 2]
    elif boundaryCondition == "advection":
        u[:, -1] = 2 * u[:, -2] - u[:, -3]
        u[:, 0] = 2 * u[:, 1] - u[:, 2]
        u[1, 1] = (u[0, 2] - u[0,0]) / (2 * delx)
        u[1, -2] = (u[0, -3] - u[0, -1]) / (2 * delx)
    elif boundaryCondition == "FDstencil":
        u[0, 0] = u[0, 2] - 2 * delx * u[1, 1]
        u[0, -1] = u[0, -3] - 2 * delx * u[1, -2]
        u[1, -1] = 2 * u[1, -2] - u[1, -3]
        u[1, 0] = 2 * u[1, 1] - u[1, 2]
        # u[0, 0] = (alpha * (2 * u[0, 1] - u[0, 2]) + u[0, 2] + (2 * phi_old_left - 2 * u[0, 1]) / alpha) / (alpha + 1)
        # u[0, -1] = (alpha * (2 * u[0, -2] - u[0, -3]) + u[0, -3] + (2 * phi_old_right - 2 * u[0, -2]) / alpha) / (alpha + 1)
        # u[1, -1] = u[1, -2] + (u[1, -2] - u[1, -3])
        # u[1, 0] = u[1, 1] + (u[1, 1] - u[1, 2])


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

        misc.savedata(f, (t, u[0], u[1]))
        tstep = 1
        while tstep < timesteps:
            k1 = calcRHS(u, delx, xpoints, boundaryCondition)
            k2 = calcRHS(u + 0.5 * delt * k1, delx, xpoints, boundaryCondition)
            k3 = calcRHS(u + 0.5 * delt * k2, delx, xpoints, boundaryCondition)
            k4 = calcRHS(u + delt * k3, delx, xpoints, boundaryCondition)

            u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

            boundaryConditions(u, boundaryCondition, delx)
            # print(u[0, 1], u[0, -2], u[1, 1], u[1, -2])

            # Advance time
            tstep = tstep + 1
            t = t + delt
            misc.savedata(f, (t, u[0], u[1]))

    else:
        # Ghost Points according to boundary conditions:
        boundaryConditions(u, "advection", delx)

        misc.savedata(f, (t, u[0], u[1]))
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
        print(u[0, 1], u[0, -2], u[1, 1], u[1, -2])

        # Advance time
        t = t + delt
        misc.savedata(f, (t, u[0], u[1]))

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

            # print(u[0, 1], u[0, -2], u[1, 1], u[1, -2])

            # Advance time
            t = t + delt
            misc.savedata(f, (t, u[0], u[1]))
            tstep = tstep + 1


    f.close()
    times, phiarray, piarray = misc.readdata(fileName, xpoints, lines=linestoread)
    return xarray, times, phiarray, piarray
