#!/usr/bin/env python3

import os
import sys
import argparse
from argparse import RawTextHelpFormatter


class doseSymOrderCSV:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Generates order_list.csv for RELION from *.tlt file according to dose-symmetric tilting scheme. (Default 2 pos, 2 neg generates: 0,1,-1,-2,2,3,-3,4...)",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input *.tlt file.")
        add('--o', help="Output csv file.")
        add('--nr_pos', type=str, default="2",
            help="Number of consecutive positive tilts in output. Default: 2")
        add('--nr_neg', type=str, default="2",
            help="Number of consecutive negative tilts in output. Default: 2")
        add('--pacetomo', dest='pacetomo', action='store_true', default=False,
            help="Pace-tomo style dose symmetric (0, ++, --, ++, --.....)")

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

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        inFileName = args.i
        outFileName = args.o
        nrPosDirection = int(args.nr_pos)
        nrNegDirection = int(args.nr_neg)

        with open(inFileName) as file:
            lines = file.readlines()
            all_tilts = []
            negative_tilt = []
            positive_tilt = []

            for line in lines:
                all_tilts.append(float(line))

            nr_of_tilts = len(all_tilts)
            for i in range(round(nr_of_tilts / 2)):
                negative_tilt.append(all_tilts.pop(0))
            if args.pacetomo:
                zerotilt = all_tilts.pop(0)
            for i in range(len(all_tilts)):
                positive_tilt.append(all_tilts.pop(0))

        dosesym_ordered = []

        if args.pacetomo:
            dosesym_ordered.append(zerotilt)

        while len(negative_tilt) > 0 or len(positive_tilt) > 0:
            for i in range(nrPosDirection):
                if len(positive_tilt) > 0:
                    dosesym_ordered.append(positive_tilt.pop(0))
            for i in range(nrNegDirection):
                if len(negative_tilt) > 0:
                    dosesym_ordered.append(negative_tilt.pop(-1))

        with open(outFileName, 'w') as outfile:
            i = 1
            for dosesym_item in dosesym_ordered:
                outfile.write("%d, %.2f\n" % (i, dosesym_item))
                i += 1

        print("File %s created." % outFileName)
        print("Have fun!")


if __name__ == "__main__":
    doseSymOrderCSV().main()
