#!/usr/bin/env python3
import argparse
import sys
import os
from lib.newstlist import *

def getCtfFind4Results(diagnosticTxtFile):
    ctffind_results = []
    with open(diagnosticTxtFile, 'r') as diagnosticFile:
        for line in diagnosticFile.readlines():
            if line.split()[0] != "#":
                ctffind_results.append([float(i) for i in line.split()])
    return  ctffind_results

def writeCtfFind4ResultsSumary(ctfFind4Results, outputFile):
    with open(outputFile, 'w') as ctfFind4outFile:
        for ctfFind4Result in ctfFind4Results:
            ctfFind4outFile.write("%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n" % tuple(ctfFind4Result))

def main(diagnosticTxtFile, outputTxtFile, thresholdValue):

    print("Processing:  %s...." % diagnosticTxtFile)

    ctfFind4Data = getCtfFind4Results(diagnosticTxtFile)

    tiltCounter = 1


    if os.path.exists(outputTxtFile):
        print("Previous %s file found. The excluded tilts present in it will be extended by the new ones...." % outputTxtFile)
        with open(outputTxtFile, 'r') as tiltListFile:
            rangeList = tiltListFile.read()
        excludeTilts =  extractMembersOfRangeList(rangeList)
    else:
        excludeTilts = []

    for ctfFindTilt in ctfFind4Data:
        if ctfFindTilt[6] >= thresholdValue:
            excludeTilts.append(tiltCounter)
        tiltCounter += 1

    ctfFind4DataOut = [ ctfFindTilt for ctfFindTilt in ctfFind4Data if ctfFindTilt[6] < thresholdValue]

    with open(outputTxtFile, 'w') as outTiltListFile:
        outTiltListFile.write("%s\n" % rangeListCreate(excludeTilts))

    writeCtfFind4ResultsSumary(ctfFind4DataOut, outputTxtFile+".ctf")

    print("Tilts with CTF ring success fit over the threshold: %s" % rangeListCreate(excludeTilts))
    print("%s written " % outputTxtFile)
    print("All done! Have fun!")


if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Analyze the Ctffind4 output txt file and write out text file with tilt-numbers that have CTF rings fitted over the threshold. Tilts are numbered from 1.")
    add = parser.add_argument
    add('--i', help="Input Ctffind4 diagnostic txt file.")
    add('--o', help="Output text file with tilt numbers under the threshold.")
    add('--threshold', type=float, default=100,
        help="Threshold (in Angstroms) value for the resolution up to which CTF rings were fit successfully. (Default: 100)")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()
    main(args.i, args.o, args.threshold)
