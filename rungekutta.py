from typing import Callable


def RK4(f: Callable, u0, delt, t0, tmax, fileName: str):
    f = open(fileName, "a")
    f.write("New Run")
    iter = (tmax - t0) / delt

    u = u0
    t = t0

    for it in range(iter):
        k1 = f(t, u)
        k2 = f(t + delt / 2, u0 + delt * k1 / 2)
        k3 = f(t + delt / 2, u0 + delt * k2 / 2)
        k4 = f(t + delt, u0 + k3)
        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)
        t += delt

        f.write("Iteration: " + it + "   Time: " + t + "  u: " + u)

    f.close()
