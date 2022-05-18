import numpy as np
import funcsandder


def calcRHS(phi, pi, delx, xsteps):
    result = np.zeros((2, xsteps + 1), dtype=np.double)
    result[2, 1:-1] = (phi[2:] - 2 * phi[1:-1] + phi[:-2]) / (delx * delx)
    result[1, 1:-1] = pi[1:-1]

    # Boundary Conditions by copying:
    result[:, 0] = result[:, 1]
    result[:, -1] = result[:, -2]

    return result


def solving(xsteps=100, fileName="calculateddata"):
    t = 0
    tsteps = 100
    tmax = 1
    delt = tmax / tsteps
    x0 = 0
    xmax = 1
    delx = (xmax - x0) / xsteps

    f = open(fileName, "a")
    f.write("\n\nNew Run")

    # Set initial values to one of the function s,g. 0 so far.
    pi = funcsandder.s1(np.array(range(0, 1, delx)))
    phi = funcsandder.s1(np.array(range(0, 1, delx)))

    while t < tmax:
        rhs = np.zeros((2, xsteps + 2), dtype=np.double)
        rhs = calcRHS(phi, pi, delx, xsteps)

        xstep = 0
        while xstep <= xsteps:
            # Berechnung der Ki:
            kphi1 = rhs[1, xstep]
            kpi1 = rhs[2, xstep]
            kphi2 = rhs[1, xstep] + delt * kphi1 / 2
            kpi2 = rhs[2, xstep] + delt / 2 * 2  # TODO
            kphi3 = rhs[1, xstep] + delt * kphi2 / 2
            kpi3 = rhs[2,]
            kphi4 = rhs[1, xstep] + delt * kphi3
            kpi4 = rhs[2,]

            pi[xstep + 1] = pi[xstep + 1] + delt * (kpi1 / 6 + kpi2 / 3 + kpi3 / 3 + kpi4 / 6)
            phi[xstep + 1] = phi[xstep + 1] + delt * (kphi1 / 6 + kphi2 / 3 + kphi3 / 3 + kphi4 / 6)
            x = x + delx

        t = t + delt
