import a8
import a10
import funcsandder
import numpy as np
import testa10ram
import a8



#xsteps
x0=-200
xmax=200
xpoints=2000


#timesteps
periods=200
timesteps=periods*1000#number of timesteps
t0=0 #starting time


#alpha
alpha=1

#Initial conditions
xarray = np.array(np.linspace(x0, xmax, xpoints), dtype=np.double)
phiinit=funcsandder.gausswave(xarray,150,3)#phi
piinit=-funcsandder.dergaus(xarray,150,3)#pi
# phiinit = funcsandder.s1(xarray)
# piinit=np.zeros(xpoints,dtype=np.double)

#Boundary condition
boundaryCondition="FDstencil"

#read data
fileName="C:\\Users\\Programmieren\\PycharmProjects\\Wave.Equation\\calculateddata.txt"
linestoread = [0]
for i in range(1, periods):
    linestoread.append(i * 1000)




# a8.stabtest(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)

# testa10ram.calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)
# testa10ram.plotten(xpoints,xarray,linestoread,fileName)

a8.convtestadvfd(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)