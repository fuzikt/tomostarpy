#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import struct
import numpy as np
import subprocess
from lib.matrix3 import *
import shutil


class convmapModGetParticles:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Get particle parameters (cc, coordinates, Euler angles) from a convmap file at positions defined by a mod file. Output is writen in emClarity csv format.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--convmap', help="Input _convmap.mrc file.")
        add('--mod', type=str, default="", help="Input mod file.")
        add('--o', type=str, default="", help="Output emClarity csv file")

    def usage(self):
        self.parser.print_help()

    def error(self, *msgs):
        self.usage()
        print("Error: " + '\n'.join(msgs))
        print(" ")
        sys.exit(2)

    def validate(self, args):
        if shutil.which("model2point") == None:
            self.error("model2point program not found. Please source IMOD to make model2point available.")

        if len(sys.argv) == 1:
            self.error("No input file given.")

        if not os.path.exists(args.convmap):
            self.error("Input file '%s' not found."
                       % args.convmap)

        if not os.path.exists(args.mod):
            self.error("Input file '%s' not found."
                       % args.mod)

        if args.convmap == "":
            self.error("No input _convmap.mrc was given.")

        if args.mod == "":
            self.error("No input mod file was given.")

        if args.o == "":
            self.error("No output csv file was specified.")

    def readMrcSizeApix(self, mrcFileName):
        with open(mrcFileName, "rb") as mrcFile:
            imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
            mrcFile.seek(28)
            mx = int(struct.unpack('i', mrcFile.read(4))[0])
            mrcFile.seek(40)
            xlen = float(struct.unpack('f', mrcFile.read(4))[0])
            apix = xlen / mx
            mrcFile.seek(196)
            originX = float(struct.unpack('f', mrcFile.read(4))[0])
            originY = float(struct.unpack('f', mrcFile.read(4))[0])
            originZ = float(struct.unpack('f', mrcFile.read(4))[0])
        return imageSizeX, imageSizeY, imageSizeZ, apix, originX, originY, originZ

    def readMrcData(self, mrcFileName):
        with open(mrcFileName, "rb") as mrcFile:
            imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
            mrcData = np.fromfile(mrcFile, dtype=np.dtype(np.float32), count=(imageSizeX * imageSizeY * imageSizeZ),
                                  offset=1024 - 12)
        return mrcData

    def readAngleListFile(self, angleListFile):
        anglesList = []
        with open(angleListFile) as file:
            for line in file:
                anglesList.append([float(line.split()[0]), float(line.split()[1]),
                                   float(line.split()[2])])
        return anglesList

    def readSelectedParticlesMODFile(self, selectedParticleMODFile):
        # convert mod file to coordinates txt
        sys.stdout.flush()
        command = "model2point -inp %s -ou tmp_model_coords.txt > /dev/null " % selectedParticleMODFile
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

        selectedParticles = []
        with open("tmp_model_coords.txt") as file:
            for line in file:
                selectedParticles.append([float(line.split()[0]), float(line.split()[1]),
                                          float(line.split()[2])])

        os.remove("tmp_model_coords.txt")

        return selectedParticles

    def writeEmClarity(self, outputEmClarityCSV, emCparticles):
        with open(outputEmClarityCSV, 'w') as file:
            for particle in emCparticles:
                file.write(particle)
                file.write('\n')

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        binning = int(args.convmap.split("_")[-2].replace("bin", ""))

        imageSizeX, imageSizeY, imageSizeZ, apix, originX, originY, originZ = self.readMrcSizeApix(args.convmap)

        convmapPath = "/".join(args.convmap.split("/")[:-1])
        convmapName = args.convmap.split("/")[-1]

        convmapMrc = self.readMrcData(convmapPath + "/" + convmapName)
        anglesMrc = self.readMrcData(convmapPath + "/" + convmapName.replace("convmap", "angles"))
        anglesList = self.readAngleListFile(convmapPath + "/" + convmapName.replace("convmap.mrc", "angles.list"))

        particles = self.readSelectedParticlesMODFile(args.mod)

        outputParticleList = []
        particleCounter = 1

        for particle in particles:
            # data[x][y][z] = data[x + ymaxX + zmaxX*maxY]
            matrixIndex = int(round(particle[0])) - 1 + (int(round(particle[1])) - 1) * imageSizeX + (
                    int(round(particle[2])) - 1) * imageSizeX * imageSizeY
            CC = convmapMrc[matrixIndex]
            angleListIndex = int(anglesMrc[matrixIndex]) - 1
            eulers = anglesList[angleListIndex]
            rotMatrix = matrix_from_euler_zxz(radians(eulers[0]), radians(eulers[1]), radians(eulers[2]))
            x = particle[0] * binning
            y = particle[1] * binning
            z = particle[2] * binning

            # matrix elements negated because of different emClarity ZXZ convention
            outputParticleList.append(
                "%.2f %.0d 0 %0.d 1 1 1 1 1 0 %.6f %.6f %.6f %.2f %.2f %.2f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f 1" % (
                    CC, binning, particleCounter, x, y, z, eulers[0], eulers[1], eulers[2], rotMatrix.m[0][0],
                    -rotMatrix.m[0][1], rotMatrix.m[0][2], -rotMatrix.m[1][0], rotMatrix.m[1][1], -rotMatrix.m[1][2],
                    rotMatrix.m[2][0], -rotMatrix.m[2][1], rotMatrix.m[2][2]))

            particleCounter += 1

        self.writeEmClarity(args.o, outputParticleList)

        print("%s particles processed" % (particleCounter-1))

        print("%s file written" % args.o)

        print("All done. Have fun!")


if __name__ == "__main__":
    convmapModGetParticles().main()
