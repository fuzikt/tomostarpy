#!/usr/bin/env python3
import argparse
import sys


def main(inputTomoStar, inputTltDose, outputTomoStar):
    tomoStarLines = []
    with open(inputTomoStar, 'r') as tomoStarFile:
        for line in tomoStarFile.readlines():
            if len(line.split()) > 2:
                tomoStarLines.append(line.split())

    tltDoseLines = []
    with open(inputTltDose, 'r') as tltDoseFile:
        for line in tltDoseFile.readlines():
            tltDoseLines.append(line.split())

    for tomoStarLine in tomoStarLines:
        for doseLine in tltDoseLines:
            if round(float(tomoStarLine[1]) * 10) == round(-float(doseLine[0]) * 10):
                tomoStarLine[5] = doseLine[1]
                break

    with open(outputTomoStar, 'w') as outTomoStarFile:
        outTomoStarFile.write("\n")
        outTomoStarFile.write("data_\n")
        outTomoStarFile.write("\n")
        outTomoStarFile.write("loop_\n")
        outTomoStarFile.write("_wrpMovieName #1\n")
        outTomoStarFile.write("_wrpAngleTilt #2\n")
        outTomoStarFile.write("_wrpAxisAngle #3\n")
        outTomoStarFile.write("_wrpAxisOffsetX #4\n")
        outTomoStarFile.write("_wrpAxisOffsetY #5\n")
        outTomoStarFile.write("_wrpDose #6\n")
        for tomoStarLine in tomoStarLines:
            outTomoStarFile.write("%s %s %s %s %s %s\n" % (
            tomoStarLine[0], tomoStarLine[1], tomoStarLine[2], tomoStarLine[3], tomoStarLine[4], tomoStarLine[5]))


if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Replace dose in warp .tomostar by value from .tltdose file.")
    add = parser.add_argument
    add('--i', help="Input tomostar file.")
    add('--itlt', help="Input dosetlt file.")
    add('--o', help="Output tomostar file.")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()

    main(args.i, args.itlt, args.o)
