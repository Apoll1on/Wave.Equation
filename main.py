import funcsandder
import numpy as np

import a8

# xsteps
x0 = 0
xmax = 1
xpoints = 501

# timesteps
timesteps = 1000  # number of timesteps, factor has to be the same as in " #read data"
tcount = 10  # Number of lines to read out/plot
t0 = 0  # starting time

# alpha
alpha = 1

# Initial conditions
xarray = np.linspace(x0, xmax, xpoints)
# phiinit=funcsandder.gausswave(xarray,-20.,0.05)#phi
# piinit=funcsandder.dergaus(xarray,-20,0.05)#pi
phiinit = funcsandder.s1(xarray)
piinit = np.zeros(xpoints, dtype=np.double)

# Boundary condition
boundaryCondition = "FDstencil"

# read data
fileName = "calculateddata.txt"
linestoread = [0]
for i in range(1, tcount + 1):
    linestoread.append(i * timesteps / tcount)

a8.convergence(x0, xmax, xpoints, t0, timesteps, alpha, phiinit, piinit, boundaryCondition, fileName, linestoread, 15)

# testa10new.calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)
# testa10new.plotten(xpoints,xarray,linestoread,fileName)
# a8.stabtest()
