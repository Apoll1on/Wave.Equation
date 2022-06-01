import numpy as np
import funcsandder as fad
from matplotlib import pyplot as plt

n = 100
xaxis = np.linspace(0, 1, n + 1)
xaxis2 = np.linspace(0, 1, 2 * n + 1)
xaxis4 = np.linspace(0, 1, 4 * n + 1)

def d1_s1(x, x2):
    fh = fad.DUDX(fad.s1(x))
    fh2 = fad.DUDX(fad.s1(x2))[::2]
    return np.abs(fad.Ds1DX(x) - fh) / np.abs(fad.Ds1DX(x) - fh2)

def nr2(x, x2):
    fh = fad.DUDX(fad.s1(x))
    fh2 = fad.DUDX(fad.s1(x2))[::2]
    return fad.Ds1DX(x) - fh

def plot():
    fh = fad.DUDX(fad.s1(xaxis))
    fh2 = fad.DUDX(fad.s1(xaxis2))[::2]
    fh4 = fad.DUDX(fad.s1(xaxis4))[::4]
    Q = np.abs(fad.Ds1DX(xaxis) - fh) / np.abs(fad.Ds1DX(xaxis) - fh2)
    # plt.plot(xaxis[1:-1], Q[1:-1], 'ro', xaxis[1:-1], fad.s1(xaxis)[1:-1], markersize = 3)
    # plt.show()
    # plt.plot(xaxis[1:-1], fad.Ds1DX(xaxis)[1:-1] - fh[1:-1], xaxis[1:-1], 4 * (fad.Ds1DX(xaxis)[1:-1] - fh2[1:-1]), 'ro', markersize = 1, linewidth = 7)
    # plt.show()
    plt.plot(xaxis[1:-1], fh[1:-1] - fh2[1:-1], xaxis[1:-1], 4 * (fh2[1:-1] - fh4[1:-1]), 'ro', markersize=1, linewidth=7)
    plt.show()
