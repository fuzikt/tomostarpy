#!/usr/bin/env python3
import argparse
import math
from lib.matrix3 import *
from lib.metadata import MetaData
from argparse import RawTextHelpFormatter


class placebackPDB:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Performs placeback of a PDB strcture into full tomogram volume according to the coordinates and euler angles defined in a star file.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input star file file.")
        add('--i_pdb', help="Input star file file.")
        add('--o', help="Output prefix")
        add('--center_offset', default="0,0,0", help="PDB center offset. (default: 0,0,0)")
        add('--backbone', dest='backbone', action='store_true', default=False,
            help="Store only backbone atoms in the output.")
        add('--tomo_name', type=str, default="",
            help="Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back")
        add('--separate_pdb', dest='separate_pdb', action='store_true', default=False,
            help="Write separate PDB file per particle")

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

        if not os.path.exists(args.i_pdb):
            self.error("Input PDB file '%s' not found."
                       % args.i_pdb)

        if args.o == "":
            self.error("No output prefix was specified.")

    def rotateShiftPDBcoordiantes(self, inputPDB, rot, tilt, psi, transX, transY, transZ, offset_x, offset_y, offset_z, backboneOnly):
        rotatedPDB = []

        rotation_matrix = matrix_from_euler(rot, tilt, psi)

        for inputPDBline in inputPDB:
            if inputPDBline.split()[0] == "ATOM":
                if (inputPDBline[13:15] == "N " or inputPDBline[13:15] == "H " or   inputPDBline[13:15] == "CA" or  inputPDBline[13:15] == "HA" or inputPDBline[13:15] == "C " or inputPDBline[13:15] == "O ") or not backboneOnly:
                    X = float(inputPDBline[30:38]) - offset_x
                    Y = float(inputPDBline[38:46]) - offset_y
                    Z = float(inputPDBline[46:54]) - offset_z

                    Xrot = (rotation_matrix.m[0][0] * X + rotation_matrix.m[0][1] * Y + rotation_matrix.m[0][2] * Z) + transX
                    Yrot = (rotation_matrix.m[1][0] * X + rotation_matrix.m[1][1] * Y + rotation_matrix.m[1][2] * Z) + transY
                    Zrot = (rotation_matrix.m[2][0] * X + rotation_matrix.m[2][1] * Y + rotation_matrix.m[2][2] * Z) + transZ

                    newPDBline = inputPDBline[:30] + "{:8.3f}".format(Xrot) + "{:8.3f}".format(Yrot) + "{:8.3f}".format(
                        Zrot) + inputPDBline[55:]
                    rotatedPDB.append(newPDBline)
            else:
                rotatedPDB.append(inputPDBline)

        return rotatedPDB

    def readPDB(self, pdbFileName):
        with open(pdbFileName) as pdbFile:
            pdbLines = [line for line in pdbFile]
        return pdbLines

    def writePDB(self, pdbFileName, pdbContent):
        with open(pdbFileName, 'w') as pdbFile:
            for line in pdbContent:
                pdbFile.write(line)

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

        offset_x, offset_y, offset_z = [float(i) for i in args.center_offset.split(",")]

        backboneOnly = args.backbone
        tomoName = args.tomo_name

        # read in pdb file
        inputPDB = self.readPDB(args.i_pdb)

        newPDB = []

        particle_count = 1

        for particle in md:
            if tomoName == "" or tomoName == particle.rlnTomoName:
                rot = math.radians(particle.rlnAngleRot)
                tilt = math.radians(particle.rlnAngleTilt)
                psi = math.radians(particle.rlnAnglePsi)

                transX = particle.rlnCoordinateX * apix - particle.rlnOriginXAngst
                transY = particle.rlnCoordinateY * apix - particle.rlnOriginYAngst
                transZ = particle.rlnCoordinateZ * apix - particle.rlnOriginZAngst

                newPDB += self.rotateShiftPDBcoordiantes(inputPDB, rot, tilt, psi, transX, transY, transZ, offset_x,
                                                         offset_y, offset_z, backboneOnly)
                if args.separate_pdb:
                    self.writePDB(args.o+"_"+str(particle_count)+".pdb", newPDB)
                    particle_count += 1
                    newPDB = []

        if not args.separate_pdb:
            self.writePDB(args.o+".pdb", newPDB)


if __name__ == "__main__":
    placebackPDB().main()