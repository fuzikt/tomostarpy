#!/usr/bin/env python3
import argparse
import math
from lib.matrix3 import *
from lib.metadata import MetaData
from argparse import RawTextHelpFormatter
from lib.chimcmm import cmmData
from lib.chimbild import bildPolygonData


class placebackPolygon:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Performs placeback of oriented arrows according to the coordinates and euler angles defined in a star file.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input star file file.")
        add('--o', help="Output prefix")
        add('--n', default="3", type=int,
            help="Order of the polygon (3 - triangle; 4 - square; 5 - pentagon etc. (default: 3)")
        add('--size', default="200", type=float,
            help="Distance between the center and the vertex in Angstroms. (default: 200)")
        add('--pre_rot', default="0.0", type=float,
            help="Pre rotation (in degrees) applied to the polygon before placing back. (default: 0)")
        add('--cmm', dest='cmm', action='store_true', default=False,
            help="Generate cmm file as well.")
        add('--color_lb', type=str, default="",
            help="Label from the star file that will be used for rainbow coloring of the arrows and cmm markers. (default: empty")
        add('--tomo_name', type=str, default="",
            help="Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back")
        add('--single', dest='single', action='store_true', default=False,
            help="Creates a single polygon with applied ---pre_rot and of size --size. Useful to find the --pre_rot parameter. No input star file needed.")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):
        if not args.single and len(sys.argv) == 1:
            self.error("No input file given.")

        if not args.single and not os.path.exists(args.i):
            self.error("Input file '%s' not found."
                       % args.i)

        if args.o == "":
            self.error("No output prefix was specified.")

    def rotatePolygon(self, X, Y, Z, rot, tilt, psi):

        rotation_matrix = matrix_from_euler(rot, tilt, psi)
        Xrot = (rotation_matrix.m[0][0] * X + rotation_matrix.m[0][1] * Y + rotation_matrix.m[0][2] * Z)
        Yrot = (rotation_matrix.m[1][0] * X + rotation_matrix.m[1][1] * Y + rotation_matrix.m[1][2] * Z)
        Zrot = (rotation_matrix.m[2][0] * X + rotation_matrix.m[2][1] * Y + rotation_matrix.m[2][2] * Z)

        return Xrot, Yrot, Zrot

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        # create polygon vertex coordinates by rotating a vector od "size" length
        polygonVertices = []
        rotStep = radians(360 / args.n)

        for i in range(args.n):
            rotAngle = (i * rotStep) + radians(args.pre_rot)
            Xrot, Yrot, Zrot = self.rotatePolygon(args.size, 0, 0, rotAngle, 0, 0)
            polygonVertices.append([Xrot, Yrot, Zrot])

        if args.single:
            outBildPolygonItems = bildPolygonData()
            outBildPolygonItems.addItem(polygonVertices, 0, 0, 0, 0)
            outBildPolygonItems.writeBuildPolygonFile(args.o + ".bild")
            exit()

        md = MetaData(args.i)
        optic_groups = []
        for optic_group in md.data_optics:
            optic_groups.append(optic_group)

        # get unbinned apix from star file
        apix = float(optic_groups[0].rlnTomoTiltSeriesPixelSize)

        tomoName = args.tomo_name

        outCmmItems = cmmData(tomoName)
        outBildPolygonItems = bildPolygonData()

        for particle in md:
            if tomoName == "" or tomoName == particle.rlnTomoName:
                shiftX = particle.rlnCoordinateX * apix - particle.rlnOriginXAngst
                shiftY = particle.rlnCoordinateY * apix - particle.rlnOriginYAngst
                shiftZ = particle.rlnCoordinateZ * apix - particle.rlnOriginZAngst

                rotatedPolygonVertices = []

                for vertex in polygonVertices:
                    rot = math.radians(particle.rlnAngleRot)
                    tilt = math.radians(particle.rlnAngleTilt)
                    psi = math.radians(particle.rlnAnglePsi)

                    rotX, rotY, rotZ = self.rotatePolygon(vertex[0], vertex[1], vertex[2], rot, tilt, psi)
                    rotatedPolygonVertices.append([rotX, rotY, rotZ])


                if args.color_lb == "":
                    value = 0
                else:
                    value = getattr(particle, args.color_lb)

                if args.cmm:
                    outCmmItems.addItem(shiftX, shiftY, shiftZ, value, particle.rlnCoordinateX, particle.rlnCoordinateY,
                                        particle.rlnCoordinateZ, args.size / 5)
                outBildPolygonItems.addItem(rotatedPolygonVertices, value, shiftX, shiftY, shiftZ)

        outBildPolygonItems.writeBuildPolygonFile(args.o + ".bild")

        if args.cmm:
            outCmmItems.writeCmmFile(args.o + ".cmm")


if __name__ == "__main__":
    placebackPolygon().main()
