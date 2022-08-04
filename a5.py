import numpy as np
import funcsandder as fad
from matplotlib import pyplot as plt

n = 22
xaxis = np.linspace(0, 1, n + 1)
xaxis2 = np.linspace(0, 1, 2 * n + 1)
xaxis4 = np.linspace(0, 1, 4 * n + 1)


def nr2(x, x2):
    fh = fad.DUDX(fad.s1(x))
    fh2 = fad.DUDX(fad.s1(x2))[::2]
    return fad.Ds1DX(x) - fh

def plot():
    # fh = fad.DUDX(fad.s1(xaxis))
    # fh2 = fad.DUDX(fad.s1(xaxis2))[::2]
    # fh4 = fad.DUDX(fad.s1(xaxis4))[::4]
    # Q = np.abs(fad.Ds1DX(xaxis) - fh) / np.abs(fad.Ds1DX(xaxis) - fh2)
    # # plt.plot(xaxis[1:-1], Q[1:-1], 'ro', xaxis[1:-1], fad.s1(xaxis)[1:-1], markersize = 3)
    # # plt.show()
    # # plt.plot(xaxis[1:-1], fad.Ds1DX(xaxis)[1:-1] - fh[1:-1], xaxis[1:-1], 4 * (fad.Ds1DX(xaxis)[1:-1] - fh2[1:-1]), 'ro', markersize = 1, linewidth = 7)
    # # plt.show()
    # plt.plot(xaxis[1:-1], fh[1:-1] - fh2[1:-1], xaxis[1:-1], 4 * (fh2[1:-1] - fh4[1:-1]), 'ro', markersize=1, linewidth=7)
    # plt.show()

    cos = fad.Ds1DX(xaxis)
    der = fad.DUDX(fad.s1(xaxis))

    sin2 = fad.D2s1DX2(xaxis)
    der2 = fad.D2UDX2(fad.s1(xaxis))

    fig1, axs1 = plt.subplots()
    axs1.plot(xaxis[1:-1], cos[1:-1], label='analytic derivative')
    axs1.plot(xaxis[1:-1], der[1:-1], label='numeric derivative')
    axs1.legend()
    axs1.set_title('First derivative')
    axs1.set_ylabel('Analytic and numeric derivative')
    plt.show()

    fig2, axs2 = plt.subplots()
    axs2.set_ylabel('Analytic - numeric derivative')
    axs2.legend()
    plt.show()

    fig3,axs3 = plt.subplots()
    axs3.plot(xaxis[1:-1], sin2[1:-1], label='analytic derivative')
    axs3.plot(xaxis[1:-1], der2[1:-1], label='numeric derivative')
    axs3.legend()
    plt.show()

    fig4, axs4 = plt.subplots()
    axs4.legend()
    plt.show()

def convergence():
    fig, axs = plt.subplots()
    axs.plot(xaxis[1:-1], (fad.Ds1DX(xaxis) - fad.DUDX(fad.s1(xaxis)))[1:-1], label="$f' - f'^{(h)}$")
    axs.plot(xaxis[1:-1], 4 * (fad.Ds1DX(xaxis) - fad.DUDX(fad.s1(xaxis2))[::2])[1:-1], label='$4 \cdot (f\' - f\'^{ (h/2)})$')

    axs.legend()
    plt.show()

def selfconvergence():
    fig, axs = plt.subplots()
    axs.plot(xaxis[1:-1], (fad.DUDX(fad.s1(xaxis)) - fad.DUDX(fad.s1(xaxis2))[::2])[1:-1], label="$f\' - f\'^{ (h/2)}$")
    axs.plot(xaxis[1:-1], 4 * (fad.DUDX(fad.s1(xaxis2))[::2] - fad.DUDX(fad.s1(xaxis4))[::4])[1:-1], label="$4 \cdot (f\'^{ (h/2)} - f\'^{ (h/4)})$")

    axs.legend()
    plt.show()