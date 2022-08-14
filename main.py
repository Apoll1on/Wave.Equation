import matplotlib.pyplot as plt

import funcsandder
import numpy as np
import a5

import a8

# xsteps
x0 = 0
xmax = 15
xpoints = 4001

# timesteps
timesteps = 20000  # number of timesteps, factor has to be the same as in " #read data"
tcount = 2000  # Number of lines to read out/plot
t0 = 0  # starting time

# alpha
alpha = 1

# Initial conditions
xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)
phiinit = funcsandder.gausswave(xarray, 0.5, 0.05)
piinit = -funcsandder.dergaus(xarray, 0.5, 0.05)

# phiinit = funcsandder.s1(xarray)
# piinit = funcsandder.Ds1DX(xarray)
# piinit = np.zeros(xpoints)
# phiinit[int(-x0 * (xpoints - 1) / (xmax - x0)): int((1-x0) * (xpoints - 1) / (xmax - x0) + 1)] = funcsandder.s1(np.linspace(0, 1, int((xpoints - 1) / (xmax - x0) + 1)))

# Boundary condition
boundaryCondition = "FDstencil"

# read data
fileName = "calculateddata.txt"
linestoread = [0]
for i in range(1, tcount + 1):
    linestoread.append(int(i * timesteps / tcount))
print(linestoread)

# nur conv fdstencil
# a8.conv_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)
# a8.plotconvfd(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)


# für den plot mit allen 3
# a8.compareconv(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)
# a8.plotcompconv(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)

# a8.stabtest(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread)

# a5.plot()

# for selfconv
a8.selfconv_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, "FDstencil", fileName, linestoread)
# a8.selfconv_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, "advection", fileName, linestoread)
# a8.selfconv_FD(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, "extrapolation", fileName, linestoread)
