#!/usr/bin/env python3
import argparse
import math
from lib.matrix3 import *
from lib.metadata import MetaData
from argparse import RawTextHelpFormatter
from lib.chimcmm import cmmData
from lib.chimbild import bildArrowData


class placebackArrow:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Performs placeback of oriented arrows (*.bild format) according to the coordinates and euler angles defined in a star file.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input star file file.")
        add('--o', help="Output prefix")
        add('--length', default="200.0", type=float, help="Length of the arrow in Angstroms. (default: 200)")
        add('--thickness', default="0", type=float,
            help="Radius of the arrow base in Angstroms. (default: 1/20 of the arrow length)")
        add('--cmm', dest='cmm', action='store_true', default=False,
            help="Generate cmm file as well.")
        add('--color_lb', type=str, default="",
            help="Label from the star file that will be used for rainbow coloring of the arrows and cmm markers. (default: empty)")
        add('--invert', dest='invert', action='store_true', default=False,
            help="Invert the pointing direction of the arrow.")
        add('--tomo_name', type=str, default="",
            help="Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back")

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

        if args.o == "":
            self.error("No output prefix was specified.")

    def rotateArrow(self, X, Y, Z, rot, tilt, psi):

        rotation_matrix = matrix_from_euler(rot, tilt, psi)
        Xrot = (rotation_matrix.m[0][0] * X + rotation_matrix.m[0][1] * Y + rotation_matrix.m[0][2] * Z)
        Yrot = (rotation_matrix.m[1][0] * X + rotation_matrix.m[1][1] * Y + rotation_matrix.m[1][2] * Z)
        Zrot = (rotation_matrix.m[2][0] * X + rotation_matrix.m[2][1] * Y + rotation_matrix.m[2][2] * Z)

        return Xrot, Yrot, Zrot

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        md = MetaData(args.i)
        optic_groups = []
        for optic_group in md.data_optics:
            optic_groups.append(optic_group)

        # get unbinned apix from star file
        apix = float(optic_groups[0].rlnTomoTiltSeriesPixelSize)

        if args.thickness == 0:
            thickness = args.length / 20
        else:
            thickness = args.thickness

        if args.invert:
            length = - args.length
        else:
            length = args.length

        tomoName = args.tomo_name

        outCmmItems = cmmData(tomoName)
        outBildItems = bildArrowData()

        for particle in md:
            if tomoName == "" or tomoName == particle.rlnTomoName:
                rot = math.radians(particle.rlnAngleRot)
                tilt = math.radians(particle.rlnAngleTilt)
                psi = math.radians(particle.rlnAnglePsi)

                rotX, rotY, rotZ = self.rotateArrow(0, 0, length, rot, tilt, psi)

                shiftX = particle.rlnCoordinateX * apix - particle.rlnOriginXAngst
                shiftY = particle.rlnCoordinateY * apix - particle.rlnOriginYAngst
                shiftZ = particle.rlnCoordinateZ * apix - particle.rlnOriginZAngst

                if args.color_lb == "":
                    value = 0
                else:
                    value = getattr(particle, args.color_lb)

                if args.cmm:
                    outCmmItems.addItem(shiftX, shiftY, shiftZ, value, particle.rlnCoordinateX, particle.rlnCoordinateY,
                                        particle.rlnCoordinateZ, args.length / 5)
                outBildItems.addItem(rotX, rotY, rotZ, value, shiftX, shiftY, shiftZ, thickness)

        outBildItems.writeBuilsArrowFile(args.o + ".bild")

        if args.cmm:
            outCmmItems.writeCmmFile(args.o + ".cmm")


if __name__ == "__main__":
    placebackArrow().main()
