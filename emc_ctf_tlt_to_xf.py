#!/usr/bin/env python3
import os
import sys
import argparse
from math import *
from argparse import RawTextHelpFormatter


class ctftltToXf:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Converts emClarity *_aliX_ctf.tlt from tomoCPR to IMOD format *.xf + *.tlt file and Ctffind4 format diag. output *.txt file.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input emClarity *_aliX_ctf.tlt file.")
        add('--o', help="Output prefix for IMOD style xf/tlt file and Ctffind4 style diag. output file.")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):
        if len(sys.argv) == 1:
            self.error("No input file given.")

        if not os.path.exists(args.i):
            self.error("Input file '%s' not found."
                       % args.i)

    class TiltClass:
        index = 0
        tx = 0.0
        ty = 0.0
        r11 = 0.0
        r12 = 0.0
        r21 = 0.0
        r22 = 0.0
        angle = 0.0
        defoc = 0.0
        deltaDefoc = 0.0
        astigAngle = 0.0

    def readCtfTlt(self, filename):
        tilts = []
        # tilt = self.TiltClass()
        with open(filename) as file:
            for line in file:
                splitLine = line.split()
                tilt = self.TiltClass()
                tilt.index = int(splitLine[0])
                tilt.tx = float(splitLine[1])
                tilt.ty = float(splitLine[2])
                tilt.angle = float(splitLine[3])
                tilt.r11 = float(splitLine[6])
                tilt.r21 = float(splitLine[7])
                tilt.r12 = float(splitLine[8])
                tilt.r22 = float(splitLine[9])
                tilt.deltaDefoc = float(splitLine[11])
                tilt.astigAngle = float(splitLine[12])
                tilt.defoc = float(splitLine[14])
                tilts.append(tilt)
        return tilts

    def sortByTilt(self, tilts):
        sortedTilts = []
        for tiltIndex in range(len(tilts)):
            for tilt in tilts:
                if tiltIndex + 1 == tilt.index:
                    sortedTilts.append(tilt)
                    continue
        return sortedTilts

    def writeXf(self, filename, tilts):
        with open(filename, 'w') as file:
            for tilt in tilts:
                file.write("  %0.7f  %0.7f  %0.7f  %0.7f  %0.3f  %0.3f" % (
                    tilt.r11, tilt.r21, tilt.r12, tilt.r22, tilt.tx, tilt.ty,))
                file.write('\n')
        print("%s written succesfully." % filename)

    def writeCtffind4Diag(self, filename, tilts):
        with open(filename, 'w') as file:
            tiltCount = 1
            for tilt in tilts:
                file.write(" %0.7f  %0.7f  %0.7f  %0.7f  0.0000  0.0000  0.0000" % (tiltCount,
                                                                                    -(
                                                                                            tilt.defoc - tilt.deltaDefoc) * 10e9,
                                                                                    -(
                                                                                            tilt.defoc + tilt.deltaDefoc) * 10e9,
                                                                                    degrees(tilt.astigAngle)))
                file.write('\n')
                tiltCount += 1
        print("%s written succesfully." % filename)

    def writeTlt(self, filename, tilts):
        with open(filename, 'w') as file:
            for tilt in tilts:
                file.write("%0.2f" % tilt.angle)
                file.write('\n')
        print("%s written succesfully." % filename)

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        print("Reading in %s ....." % args.i)

        tilts = self.readCtfTlt(args.i)
        print("%s tilts present in file." % len(tilts))

        sortedTilts = self.sortByTilt(tilts)
        self.writeXf(args.o + ".xf", sortedTilts)
        self.writeCtffind4Diag(args.o + ".txt", sortedTilts)
        self.writeTlt(args.o + ".tlt", sortedTilts)

        print("All done. Have fun!")


if __name__ == "__main__":
    ctftltToXf().main()
