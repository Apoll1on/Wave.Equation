import numpy as np



def calcpiasarea(sides):
    pass




def calcintegral(delx):
    xarray=np.linspace(0,1-delx,int(1/delx))
    fvalue=delx/(1-xarray*xarray)
    sum=np.sum(fvalue)
    print(sum)



calcintegral(0.000001)