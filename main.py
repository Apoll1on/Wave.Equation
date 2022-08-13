import matplotlib.pyplot as plt

import funcsandder
import numpy as np
import a5

import a8

# xsteps
x0 = 0
xmax = 1
xpoints = 4001

# timesteps
timesteps = 500  # number of timesteps, factor has to be the same as in " #read data"
tcount = 10 # Number of lines to read out/plot
t0 = 0  # starting time

# alpha
alpha = 1

# Initial conditions
xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)
phiinit = funcsandder.gausswave(xarray, 0.5, 0.05)
piinit = funcsandder.dergaus(xarray, 0.5, 0.05)

# phiinit = funcsandder.s1(xarray)
# piinit = funcsandder.Ds1DX(xarray)
# piinit = np.zeros(xpoints)
# phiinit[int(-x0 * (xpoints - 1) / (xmax - x0)): int((1-x0) * (xpoints - 1) / (xmax - x0) + 1)] = funcsandder.s1(np.linspace(0, 1, int((xpoints - 1) / (xmax - x0) + 1)))

# Boundary condition
boundaryCondition = "periodic"

# read data
fileName = "calculateddata.txt"
linestoread = [0]
for i in range(1, tcount + 1):
    linestoread.append(i * timesteps / tcount)
print(linestoread)

a8. conv_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)

# testa10new.calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)
# testa10new.plotten(xpoints,xarray,linestoread,fileName)

# a5.plot()