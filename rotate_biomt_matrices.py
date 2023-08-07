#!/usr/bin/env python3

import os
import sys
import argparse
import lib.matrix3 as matrix3


class rotBiomtMat:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Rotate and shift REMARK 350 BIOMT matrices from PDB file, using a given matrix or ZYZ convention euler angles.")
        add = self.parser.add_argument
        add('--i', help="Input PDB file.")
        add('--o', help="Output file.")
        add('--i_mat', help="Input file with rotation and translation matrix in BIOMT format used for the rotation. If defined --rot, --tilt and --psi are ignored.")
        add('--rot', type=float, default="0.00",
            help="Euler angle ROT in degrees.  (Default 0.0)")
        add('--tilt', type=float, default="0.00",
            help="Euler angle TILT in degrees.  (Default 0.0)")
        add('--psi', type=float, default="0.00",
            help="Euler angle PSI in degrees.  (Default 0.0)")
        add('--x', type=float, default="0.00",
            help="X shift in Angstroms.  (Default 0.0)")
        add('--y', type=float, default="0.00",
            help="Y shift in Angstroms.  (Default 0.0)")
        add('--z', type=float, default="0.00",
            help="Z shift in Angstroms.  (Default 0.0)")

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

    def readBIOMTmatrices(self, inputPDB):
        rotMatrices = []
        transMatrices = []
        with open(inputPDB, 'r') as pdbFile:
            for pdbLine in pdbFile.readlines():
                if "REMARK 350" in pdbLine:
                    remarkLine = pdbLine.split()
                    if remarkLine[2] == "BIOMT1":
                        m00,m01,m02 = (float(remarkLine[4]),float(remarkLine[5]),float(remarkLine[6]))
                        x = float(remarkLine[7])
                    elif remarkLine[2] == "BIOMT2":
                        m10, m11, m12 = (float(remarkLine[4]), float(remarkLine[5]), float(remarkLine[6]))
                        y = float(remarkLine[7])
                    elif remarkLine[2] == "BIOMT3":
                        m20, m21, m22 = (float(remarkLine[4]), float(remarkLine[5]), float(remarkLine[6]))
                        z = float(remarkLine[7])
                        rotMatrices.append(matrix3.Matrix3([m00,m01,m02,m10, m11, m12,m20, m21, m22]))
                        transMatrices.append([x, y, z])
        return rotMatrices, transMatrices

    def writeBIOMTmatrices(self, outputPDB, rotMatrices, transMatrices):
        with open(outputPDB, 'w') as outputPDBFile:
            for i in range(len(rotMatrices)):
                outputPDBFile.write("REMARK 350   BIOMT1 %3d%10.6f%10.6f%10.6f%15.5f\n" % (i+1, rotMatrices[i].m[0][0], rotMatrices[i].m[0][1], rotMatrices[i].m[0][2], transMatrices[i][0]))
                outputPDBFile.write("REMARK 350   BIOMT2 %3d%10.6f%10.6f%10.6f%15.5f\n" % (
                i+1, rotMatrices[i].m[1][0], rotMatrices[i].m[1][1], rotMatrices[i].m[1][2], transMatrices[i][1]))
                outputPDBFile.write("REMARK 350   BIOMT3 %3d%10.6f%10.6f%10.6f%15.5f\n" % (
                i+1, rotMatrices[i].m[2][0], rotMatrices[i].m[2][1], rotMatrices[i].m[2][2], transMatrices[i][2]))

    def rotateSymMatrices(self, rotSymMatrices, transSymMatrices, rotMatrix, transMatrix):
        rotatedMatrices = []
        translatedMatrices = []
        for i in range(len(rotSymMatrices)):
            #for rotSymMatrix in rotSymMatrices:
            rotatedMatrices.append(matrix3.matrix_multiply(rotSymMatrices[i], rotMatrix))

     #   for transSymMatrix in transSymMatrices:
            xnew = (rotSymMatrices[i].m[0][0] * transMatrix[0]) + (rotSymMatrices[i].m[0][1] * transMatrix[1]) + (rotSymMatrices[i].m[0][2] * transMatrix[2])
            ynew = (rotSymMatrices[i].m[1][0] * transMatrix[0]) + (rotSymMatrices[i].m[1][1] * transMatrix[1]) + (rotSymMatrices[i].m[1][2] * transMatrix[2])
            znew = (rotSymMatrices[i].m[2][0] * transMatrix[0]) + (rotSymMatrices[i].m[2][1] * transMatrix[1]) + (rotSymMatrices[i].m[2][2] * transMatrix[2])
            if i == 0:
                print(xnew, ynew, znew)
            translatedMatrices.append([transSymMatrices[i][0]+xnew,transSymMatrices[i][1]+ynew,transSymMatrices[i][2]+znew])
        return rotatedMatrices, translatedMatrices

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        biomtSymMatrices, biomtSymTrans = self.readBIOMTmatrices(args.i)

        if args.i_mat == "":
            rotMatrix = matrix3.matrix_from_euler(args.rot, args.tilt, args.psi)
            rotTrans = [args.x, args.y, args.z]
        else:
            rotMatrix, rotTrans = self.readBIOMTmatrices(args.i_mat)

        rotatedBiomtMatrices, rotatedBiomtTrans = self.rotateSymMatrices(biomtSymMatrices,biomtSymTrans,rotMatrix[0], rotTrans[0])

        self.writeBIOMTmatrices(args.o, rotatedBiomtMatrices, rotatedBiomtTrans)

        print("File %s created. Have fun!" % args.o)

if __name__ == "__main__":
    rotBiomtMat().main()