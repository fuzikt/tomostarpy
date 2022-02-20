#!/usr/bin/env python3
import os
import sys
import re
import subprocess
from math import *
from lib.metadata import MetaData
from lib.metadata import Item as starItem
import argparse
from argparse import RawTextHelpFormatter


class emcCSVtoStar:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Perform conversion from emClarity template matching csv file to Relion star file. It also recalculates the particle coordinates from emClarity partial tomogram to full tomogram.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input emClarity csv file.")
        add('--o', help="Output prefix. Prefix of the files generated by the script.")
        add('--reconsh', type=str, default="",
            help="emClarity *_recon.sh file of the input csv file. (located in emc_prject/recon/)")
        add('--tiltcom', type=str, default="",
            help="tilt.com file (generated by IMOD) of the original non-splitted tomogram used by Relion")
        add('--mod', type=str, default="",
            help="Modfile used for filtering particles form the input csv. Note: Set the \"--modbin\" binnig factor used for template matching.")
        add('--modbin', type=str, default="1.00",
            help="Binnig factor used for template matching. Needed only if modfile filtering is enabled. (Default 1)")
        add('--outbin', type=str, default="1.00",
            help="Binnig factor used for output coordinates and mod files. NOT used for star file (always unbinned)! (Default 1)")
        add('--cs', type=str, default="2.7",
            help="Cs value of the microscope. Used in opticsgroup in the output star file. (Default 2.7)")
        add('--kv', type=str, default="300.0",
            help="Accerelation voltage of the microscope. Used in opticsgroup in the output star file. (Default 300.0)")
        add('--apix', type=str, default="1.00",
            help="Apix of the unbinned tomogram. Used in opticsgroup in the output star file.")

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

        if args.mod != "" and args.modbin == "1.00":
            print("!!!WARNING: MOD file used for csv filtering while --modbin not set. Using default value 1.0.")

        if args.apix == "1.00":
            print("!!!WARNING: --apix not not set. Using default value 1.0.")

        if args.outbin == "1.00":
            print("!!!INFO: --outbin not not set. The output coordinate file and mod file will be unbinned.")

    def readEmClarityCSV(self, filename):
        particles = []
        with open(filename) as file:
            for line in file:
                particles.append(line.split())
        return particles

    def readSelectedParticlesMODFile(self, selectedParticleMODFile):
        # convert mod file to coordinates txt
        sys.stdout.flush()
        command = "model2point -inp %s -ou tmp_model_coords.txt > /dev/null " % selectedParticleMODFile
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

        selectedParticles = []
        with open("tmp_model_coords.txt") as file:
            for line in file:
                selectedParticles.append([int(round(float(line.split()[0]))), int(round(float(line.split()[1]))),
                                          int(round(float(line.split()[2])))])

        os.remove("tmp_model_coords.txt")

        return selectedParticles

    def readTiltCom(self, filename):
        xaxistilt = 0
        width = 0
        slice1 = 0
        slice2 = 0
        with open(filename) as file:
            for line in file:
                if "WIDTH" in line:
                    width = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[0])
                if "SLICE" in line:
                    slice1 = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[0])
                    slice2 = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[1])
                if "THICKNESS" in line:
                    thickness = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[0])
                if "SHIFT" in line:
                    shift1 = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[0])
                    shift2 = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[1])
                if "XAXISTILT" in line:
                    xaxistilt = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[0])
                if "FULLIMAGE" in line:
                    fullimage1 = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[0])
                    fullimage2 = float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line)[1])
        if width == 0:
            width = fullimage1

        return width, slice1, slice2, thickness, shift1, shift2, xaxistilt

    def euler_from_matrix(self, matrix):
        """converts a matrix to Eulers - as in Relion euler.cpp"""
        FLT_EPSILON = 1.19209e-07
        sign = lambda x: x and (1, -1)[x < 0]
        abs_sb = sqrt(matrix[0][2] * matrix[0][2] + matrix[1][2] * matrix[1][2])
        if (abs_sb > 16 * FLT_EPSILON):
            psi = atan2(matrix[1][2], -matrix[0][2])
            rot = atan2(matrix[2][1], matrix[2][0])
            if (abs(sin(psi)) < FLT_EPSILON):
                sign_sb = sign(-matrix[0][2] / cos(psi))
            else:
                sign_sb = sign(matrix[1][2]) if (sin(psi) > 0) else -sign(matrix[1][2])
            tilt = atan2(sign_sb * abs_sb, matrix[2][2])
        else:
            if (sign(matrix[2][2]) > 0):
                rot = 0
                tilt = 0
                psi = atan2(-matrix[1][0], matrix[0][0])
            else:
                rot = 0
                tilt = pi
                psi = atan2(matrix[1][0], -matrix[0][0])
        return rot, tilt, psi

    def writeEmClarity(self, outputEmClarityCSV, emCparticles):
        with open(outputEmClarityCSV, 'w') as file:
            for particle in emCparticles:
                file.write(" ".join(str(element) for element in particle))
                file.write('\n')

    def writeCoordinateFile(self, filename, particles, binning, x_ind, y_ind, z_ind):
        with open(filename, 'w') as file:
            for particle in particles:
                file.write(
                    "%s %s %s" % (particle[x_ind] / binning, particle[y_ind] / binning, particle[z_ind] / binning))
                file.write('\n')

    def writeModFile(self, outputCoordinateModFile, outputCoordinateFile):
        # convert mod file to coordinates txt
        sys.stdout.flush()
        command = "point2model -sphere 3 -scat -inp %s -ou %s > /dev/null" % (
            outputCoordinateFile, outputCoordinateModFile)
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

    def writeStarFile(self, outputStarFile, csValue, acVolatage, apix, tomoName, emCparticles, x_ind, y_ind, z_ind,
                      r11):
        # create output star file
        mdOut = MetaData()
        mdOut.version = "3.1"
        mdOut.addDataTable("data_optics")
        mdOut.addLabels("data_optics", "rlnOpticsGroup", "rlnOpticsGroupName", "rlnSphericalAberration", "rlnVoltage",
                        "rlnTomoTiltSeriesPixelSize")
        mdOut.addDataTable("data_particles")
        mdOut.addLabels("data_particles", ["rlnTomoName", "rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ",
                                           "rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi", "rlnOriginXAngst",
                                           "rlnOriginYAngst", "rlnOriginZAngst", "rlnCtfFigureOfMerit"])

        # create optics groups
        opticGroup = starItem()
        opticsGroups = []

        opticGroup.rlnOpticsGroup = 1
        opticGroup.rlnOpticsGroupName = "opticsGroup1"
        opticGroup.rlnSphericalAberration = csValue
        opticGroup.rlnVoltage = acVolatage
        opticGroup.rlnTomoTiltSeriesPixelSize = apix
        opticsGroups.append(opticGroup)
        mdOut.addData("data_optics", opticsGroups)

        # create particles table
        new_particles = []
        for emCparticle in emCparticles:
            particle = starItem()
            particle.rlnTomoName = tomoName
            particle.rlnCoordinateX = emCparticle[x_ind]
            particle.rlnCoordinateY = emCparticle[y_ind]
            particle.rlnCoordinateZ = emCparticle[z_ind]
            rot_matrix = [[float(emCparticle[r11]), float(emCparticle[r11 + 3]), float(emCparticle[r11 + 6])],
                          [float(emCparticle[r11 + 1]), float(emCparticle[r11 + 4]), float(emCparticle[r11 + 7])],
                          [float(emCparticle[r11 + 2]), float(emCparticle[r11 + 5]), float(emCparticle[r11 + 8])]]
            particle.rlnCtfFigureOfMerit = float(emCparticle[0])
            rot, tilt, psi = self.euler_from_matrix(rot_matrix)
            particle.rlnAngleRot = degrees(rot)
            particle.rlnAngleTilt = degrees(tilt)
            particle.rlnAnglePsi = degrees(psi)
            particle.rlnOriginXAngst = 0.0
            particle.rlnOriginYAngst = 0.0
            particle.rlnOriginZAngst = 0.0
            new_particles.append(particle)
        mdOut.addData("data_particles", new_particles)
        mdOut.write(outputStarFile)

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        inputEmClarityCSV = args.i
        selectedParticleMODFile = args.mod
        modFileBinning = float(args.modbin)
        outputEmClarityCSV = args.o + ".csv"
        outputCoordinateFile = args.o + "_coords.txt"
        outputCoordinateModFile = args.o + "_coords.mod"
        outputStarFile = args.o + ".star"
        originalTiltCom = args.reconsh
        newTiltCom = args.tiltcom
        csValue = float(args.cs)
        acVolatage = float(args.kv)
        apix = float(args.apix)

        outputCoordinatesBinning = float(args.outbin)

        # indices of x,y,z in emClarity csv file
        x_ind = 10
        y_ind = 11
        z_ind = 12
        # indices for rotation matrix element R11 in emClarity csv file
        r11 = 16

        emCparticles = self.readEmClarityCSV(inputEmClarityCSV)

        (origWIDTH, origSLICE1, origSLICE2, origTHICKNESS, origSHIFT1, origSHIFT2, origXAXISTILT) = self.readTiltCom(
            originalTiltCom)
        (newWIDTH, newSLICE1, newSLICE2, newTHICKNESS, newSHIFT1, newSHIFT2, newXAXISTILT) = self.readTiltCom(
            newTiltCom)

        offset_x = ((newWIDTH - origWIDTH) / 2 - origSHIFT1 + newSHIFT1)
        offset_y = origSLICE1 - newSLICE1
        offset_z = ((newTHICKNESS - origTHICKNESS) / 2 + origSHIFT2 - newSHIFT2)

        if selectedParticleMODFile != "":
            print("Selecting particles according to input mod file....")
            selectedEmCparticles = []
            selectedParticles = self.readSelectedParticlesMODFile(selectedParticleMODFile)

            for emCparticle in emCparticles:
                if [int(round(float(emCparticle[x_ind]) / modFileBinning)),
                    int(round(float(emCparticle[y_ind]) / modFileBinning)),
                    int(round(float(emCparticle[z_ind]) / modFileBinning))] in selectedParticles:
                    selectedEmCparticles.append(emCparticle)
            emCparticles = selectedEmCparticles

        print("Transforming particle coordinates...")
        for emCparticle in emCparticles:
            emCparticle[x_ind] = float(emCparticle[x_ind]) + offset_x
            emCparticle[y_ind] = float(emCparticle[y_ind]) + offset_y
            emCparticle[z_ind] = float(emCparticle[z_ind]) + offset_z

        self.writeCoordinateFile(outputCoordinateFile, emCparticles, outputCoordinatesBinning, x_ind, y_ind, z_ind)
        print("%s file written" % outputCoordinateFile)

        self.writeModFile(outputCoordinateModFile, outputCoordinateFile)
        print("%s file written" % outputCoordinateModFile)

        self.writeStarFile(outputStarFile, csValue, acVolatage, apix, args.o, emCparticles, x_ind, y_ind, z_ind, r11)
        print("%s file written" % outputStarFile)

        self.writeEmClarity(outputEmClarityCSV, emCparticles)
        print("%s file written" % outputEmClarityCSV)

        print("All done. Have fun!")


if __name__ == "__main__":
    emcCSVtoStar().main()
