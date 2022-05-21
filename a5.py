import numpy as np
import funcsandder as fad
from matplotlib import pyplot as plt

n = 101
xaxis = np.linspace(0, 1, n)

def d1_s1():
    return (fad.Ds1DX(xaxis) - fad.DUDX(fad.s1(xaxis), h)) / (fad.Ds1DX(xaxis) - fad.DUDX(fad.s1(xaxis), h/2))

def plot():
    # print(fad.Ds1DX(xaxis))
    # print(fad.DUDX(fad.s1(xaxis)))
    # print(fad.Ds1DX(xaxis))
    # print(fad.DUDX(fad.s1(xaxis), h/2))
    #plt.plot(xaxis[1:-1], d1_s1()[1:-1])
    #plt.show()
    plt.plot(xaxis, fad.DUDX(fad.s1(xaxis)))
    plt.show()