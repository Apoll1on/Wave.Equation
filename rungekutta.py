from typing import Callable, String


def RK4(f: Callable, u0, delt, iter, fileName: String):
    f = open(filename, "a")
    f.write("New Run")
    t = 0
    u = u0

    for _ in range(iter):
        k1 = f(t, u)
        k2 = f(t + delt / 2, u0 + delt * k1 / 2)
        k3 = f(t + delt / 2, u0 + delt * k2 / 2)
        k4 = f(t + delt, u0 + k3)
        u = u + delt * (k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6)

    function.close()
    return u
