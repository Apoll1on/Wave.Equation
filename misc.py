import numpy as np


def savedata(fileObject, data):
    """Put data in following order: time,phi,pi"""
    f = fileObject
    f.write(str(data[0]) + "#")
    for parts in data[1:]:
        for pieces in parts:
            f.write(str(pieces) + "#")
    f.write("\n")


def readdata(fileName, xpoints, lines):
    readLines = []
    f = open(fileName, "r")
    for pos, line in enumerate(f):
        if pos in lines:
            readLines.append(line)
    #print(readLines)
    a = len(lines)
    phiarray = np.zeros((len(lines), xpoints), dtype=np.double)
    piarray = np.zeros((len(lines), xpoints), dtype=np.double)
    times = []

    timeindex = 0
    phiindex = range(2, 2 + xpoints)
    piindex = range(2 + xpoints + 1, 2 + xpoints + 1 + xpoints)

    index = 0
    for line in readLines:

        spl = line.split("#")
        times.append(spl[timeindex])
        index2 = 0
        index3 = 0
        for i in phiindex:
            phiarray[index, index2] = spl[i]
            index2 = index2 + 1

        for i in piindex:
            piarray[index, index3] = spl[i]
            index3 = index3 + 1

        index = index + 1
    # print(times, phiarray, piarray)
    return (times, phiarray, piarray)
