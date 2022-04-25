#!/usr/bin/env python3
import argparse
import math
from argparse import RawTextHelpFormatter
from lib.matrix3 import *
from lib.vector3 import *
from lib.euler import *


class getZYZEulersFromMatrix():
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Calculates ZYZ convention Euler angles from a 3x3 rotation matrix.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Rotation matrix members in format: m[1,1],m[1,2],m[1,3],m[2,1],m[2,2]")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):
        if args.i == None:
            self.error("No input matrix.")
        if len((str(args.i).split(","))) != 9:
            self.error("Input must have 9 matrix members separated by comma.")

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        rotationMatrix = Matrix3([float(i) for i in str(args.i).split(",")])
        rot, tilt, psi = euler_from_matrix(rotationMatrix)
        print("Euler angles in degrees:\nRot: %0.2f \nTilt: %0.2f \nPsi: %0.2f" % (
            math.degrees(rot), math.degrees(tilt), math.degrees(psi)))


if __name__ == "__main__":
    getZYZEulersFromMatrix().main()
