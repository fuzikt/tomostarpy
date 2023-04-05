#!/usr/bin/env python3
import numpy as np
import struct
import argparse
import sys
import os
from lib.newstlist import *

def readMrcSizeApix(mrcFileName):
    with open(mrcFileName, "rb") as mrcFile:
        imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
        mrcFile.seek(28)
        mx = int(struct.unpack('i', mrcFile.read(4))[0])
        mrcFile.seek(40)
        xlen = float(struct.unpack('f', mrcFile.read(4))[0])
        apix = xlen / mx
        mrcFile.seek(196)
        originX = float(struct.unpack('f', mrcFile.read(4))[0])
        originY = float(struct.unpack('f', mrcFile.read(4))[0])
        originZ = float(struct.unpack('f', mrcFile.read(4))[0])
    return imageSizeX, imageSizeY, imageSizeZ, apix, originX, originY, originZ


def readMrcData(mrcFileName):
    with open(mrcFileName, "rb") as mrcFile:
        imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
        mrcData = np.fromfile(mrcFile, dtype=np.dtype(np.float32), count=(imageSizeX * imageSizeY * imageSizeZ),
                              offset=1024 - 12)
    return mrcData

def main(tiltSeriesMrcFile, outputTxtFile, thresholdValue):
    tiltSeriesFile = tiltSeriesMrcFile
    avgThreshold = thresholdValue
    outputFile = outputTxtFile

    print("Processing:  %s...." % tiltSeriesFile)
    imageSizeX, imageSizeY, imageSizeZ, apix, originX, originY, originZ = readMrcSizeApix(tiltSeriesFile)
    tiltSeriesData = readMrcData(tiltSeriesFile)

    imageSizeXY = imageSizeX * imageSizeY

    sliceAverages = np.mean(tiltSeriesData.reshape(-1, imageSizeXY), axis=1)

    tiltCounter = 1

    if os.path.exists(outputFile):
        print("Previous %s file found. The excluded tilts present in it will be extended by the new ones...." % outputFile)
        with open(outputFile, 'r') as tiltListFile:
            rangeList = tiltListFile.read()
        excludeTilts =  extractMembersOfRangeList(rangeList)
    else:
        excludeTilts = []

    for sliceAvg in sliceAverages:
        if sliceAvg <= avgThreshold:
            excludeTilts.append(tiltCounter)
        tiltCounter += 1

    with open(outputFile, 'w') as outTiltListFile:
        outTiltListFile.write("%s\n" % rangeListCreate(excludeTilts))

    print("Tilts with average under the threshold: %s" % rangeListCreate(excludeTilts))
    print("%s written..." % outputFile)
    print("All done! Have fun!")


if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Analyze the input tilt-series mrc stack and write out text file with tilt-numbers that have average signal under the threshold. Tilts are numbered from 1. Useful to remove dark tilts.")
    add = parser.add_argument
    add('--i', help="Input tilt-series mrc-stack.")
    add('--o', help="Output text file with tilt numbers under the threshold.")
    add('--threshold', type=float, default=0.5,
        help="Threshold value for the average value of the tilt. (Default: 0.5)")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()
    main(args.i, args.o, args.threshold)
