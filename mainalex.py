import a8
import a10
import funcsandder
import numpy as np
import testa10



#xsteps
x0=-30
xmax=10
xpoints=4001


#timesteps
periods=20
timesteps=periods*1000#number of timesteps
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
boundaryCondition="extrapolation"

#read data
fileName="calculateddata.txt"
linestoread = [0]
for i in range(1, periods):
    linestoread.append(i * 1000)




# a8.stabtest(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)

testa10.calcplot(x0,xmax,xpoints,t0,timesteps,alpha,phiinit,piinit,boundaryCondition,fileName,linestoread)
#testa10.plotten(xpoints,xarray,linestoread,fileName)
#a8.stabtest()
