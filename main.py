import funcsandder
import numpy as np
import a8

#xsteps
x0=0
xmax=1
xpoints=1001


#timesteps
tcount=10 # Number of lines to read out/plot
timesteps= tcount * 100 # number of timesteps, factor has to be the same as in " #read data"
t0=0 #starting time


#alpha
alpha=1

#Initial conditions
xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)
# phiinit=funcsandder.gausswave(xarray,-20.,0.05)#phi
# piinit=funcsandder.dergaus(xarray,-20,0.05)#pi
phiinit = funcsandder.s1(xarray)
piinit=np.zeros(xpoints,dtype=np.double)

#Boundary condition
boundaryCondition="advection"

#read data
fileName="calculateddata.txt"
linestoread = [0]
for i in range(1, tcount + 1):
    linestoread.append(i * 100)


a8.stabtest(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)

#testa10new.calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)
#testa10new.plotten(xpoints,xarray,linestoread,fileName)
#a8.stabtest()