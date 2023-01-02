#!/usr/bin/env python3

import os
import sys
import math
import argparse
import lib.matrix2 as matrix2


class rotXFfile:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Set the Euler angles (tilt, psi) of the particles along a filament according to the angles of the between the neighboring particle centers.")
        add = self.parser.add_argument
        add('--i', help="Input IMOD xf file.")
        add('--o', help="Output IMOD xf file.")
        add('--ang', type=str, default="0.00",
            help="Rotation angle in degrees. (Default 0.0)")

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

    def rotateElements(self, xfLine, angle):
        xfLine = [float(lineElement) for lineElement in xfLine.split()]
        m = matrix2.Matrix2([xfLine[0], xfLine[1], xfLine[2], xfLine[3]])
        x = xfLine[4]
        y = xfLine[5]
        rot_m = matrix2.matrix_from_angle(math.radians(angle))
        rot_m_t = matrix2.matrix_transpose(rot_m)
        m2 = matrix2.matrix_multiply(m, rot_m_t)
        x = x * rot_m.m[0][0] + y * rot_m.m[0][1]
        y = x * rot_m.m[1][0] + y * rot_m.m[1][1]

        return "%0.7f %0.7f %0.7f %0.7f %0.3f %0.3f" % (m2.m[0][0], m2.m[0][1], m2.m[1][0], m2.m[1][1], x, y)

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        inputXF = args.i
        outputXF = args.o
        angle = float(args.ang)

        with open(inputXF) as f:
            xfLines = [line.rstrip() for line in f]

        outXFlines = []

        for line in xfLines:
            outXFlines.append(self.rotateElements(line, angle))

        with open(outputXF, 'w') as f:
            f.write('\n'.join(outXFlines))

        print("New xf file %s created. Have fun!" % args.o)

if __name__ == "__main__":
    rotXFfile().main()
