import numpy as np

a = 0
b = 1


# n = 101
# h = 1 / (n - 1)


def DUDX(u):
    diff = np.zeros_like(u)
    h = 1 / (diff.size - 1)
    diff[1:-1] = (u[2:] - u[:-2]) / (2 * h)
    diff[0] = diff[1]
    diff[-1] = diff[-2]
    return diff


def D2UDX2(u):
    diff = np.zeros_like(u)
    h = 1 / (diff.size - 1)
    diff[1:-1] = (u[2:] - 2 * u[1:-1] + u[:-2]) / (h * h)
    return diff


def f1(x):
    return (x - 1 / 2) ** 2 + x


def f2(x):
    return (x - 1 / 2) ** 3 + (x - 1 / 2) ** 2 + x


def f3(x):
    return np.sqrt(x)


def s1(x):
    return np.sin(12 * np.pi * x)


def s2(x):
    return np.sin(12 * np.pi * x) ** 4


def g(x, a):
    return np.exp(-a * x ** 2)


def Df1DX(x):
    return 2 * x


def Df2DX(x):
    return 3 * (x) ** 2 - x + 3 / 4


def Df3DX(x):
    return 0.5 / np.sqrt(x)


def Ds1DX(x):
    return 12 * np.pi * np.cos(12 * np.pi * x)

def D2s1DX2(x):
    return -144 * np.pi * np.pi * np.sin(12 * np.pi * x)


def D2s1D2X(x):
    return -12 * 12 * np.pi * np.pi * np.sin(12 * np.pi * x)


def Ds2DX(x):
    return 48 * np.pi * np.sin(12 * np.pi * x) ** 3 * np.cos(12 * np.pi * x)


def DgDX(x, a):
    return -2 * a * x * np.exp(-a * (x) ** 2)


def gausswave(xarray, mu, sigma):
    return np.exp(-(xarray - mu) * (xarray - mu) / (2 * sigma * sigma)) / (np.sqrt(2 * np.pi) * sigma)


def dergaus(xarray, mu, sigma):
    return - (xarray - mu) * np.exp(-(xarray - mu) * (xarray - mu) / (2 * sigma * sigma)) / (np.sqrt(2 * np.pi) * sigma * sigma * sigma)

def gausswave_verschoben(xarray, mu, sigma, t):
    return np.exp(-(xarray + t - mu) * (xarray + t - mu) / (2 * sigma * sigma)) / (np.sqrt(2 * np.pi) * sigma)
