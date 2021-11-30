#!/usr/bin/env python

import os
import sys
import re
import subprocess
from lib.metadata import MetaData
import argparse


class SelCmmModStar:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Select particles from star file according to matching particle coordinates listed in Chimera cmm or IMOD mod file.")
        add = self.parser.add_argument
        add('--i', help="Input STAR filename with particles.")
        add('--o', help="Output STAR filename.")
        add('--cmm', type=str, default='',
            help="Chimera CMM file with desires coordinates.")
        add('--mod', type=str, default='',
            help="IMOD mod file with desired coordinates.")
        add('--bin', type=float, default=0,
            help="Binnig factor of the IMOD mod file.")

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

        if args.cmm == "" and args.mod == "":
            self.error("At least one of the --cmm or --mod has to be specified.")

        if args.cmm != "" and args.mod != "":
            self.error("You cannot specify both --cmm and --mod at the same time.")

        if args.mod != '' and args.bin == 0:
            args.bin = 1
            print("!!!Warning: Using mod file and binning factor was not specified. Using default value: 1")

    def readCmmFile(self, cmmFileName):
        selectedParticlesCoords = []
        with open(cmmFileName) as file:
            for line in file:
                if "id=" in line:
                    for line_element in line.split():
                        if "coordX" in line_element:
                            coordX = int(round(float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line_element)[0])))
                        elif "coordY" in line_element:
                            coordY = int(round(float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line_element)[0])))
                        elif "coordZ" in line_element:
                            coordZ = int(round(float(re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", line_element)[0])))
                    selectedParticlesCoords.append([coordX, coordY, coordZ])
        return selectedParticlesCoords

    def readModFile(self, modFileName, binning):
        sys.stdout.flush()
        command = "model2point -inp %s -ou tmp_model_coords.txt > /dev/null " % modFileName
        proc = subprocess.Popen(command, shell=True)
        proc.wait()

        selectedParticlesCoords = []

        with open("tmp_model_coords.txt") as file:
            for line in file:
                selectedParticlesCoords.append(
                    [int(round(float(line.split()[0]) * binning)), int(round(float(line.split()[1]) * binning)),
                     int(round(float(line.split()[2]) * binning))])

        os.remove("tmp_model_coords.txt")

        return selectedParticlesCoords

    def get_particles(self, md):
        particles = []
        for particle in md:
            particles.append(particle)
        return particles

    def selParticles(self, particles, selectedParticlesCoords):
        newParticles = []
        for particle in particles:
            if [int(round(particle.rlnCoordinateX)), int(round(particle.rlnCoordinateY)),
                int(round(particle.rlnCoordinateZ))] in selectedParticlesCoords:
                newParticles.append(particle)
        print(str(len(newParticles)) + " particles included in selection.")
        return newParticles

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        md = MetaData(args.i)

        mdOut = MetaData()

        new_particles = []

        particles = self.get_particles(md)

        optic_groups = []
        for optic_group in md.data_optics:
            optic_groups.append(optic_group)

        if args.cmm != '':
            print("Selecting particles from star file according to matching particle coordinates listed in cmm file...")
            selectedParticlesCoords = self.readCmmFile(args.cmm)
        else:
            print("Selecting particles from star file according to matching particle coordinates listed in mod file...")
            selectedParticlesCoords = self.readModFile(args.mod, args.bin)

        new_particles.extend(self.selParticles(particles, selectedParticlesCoords))

        if md.version == "3.1":
            mdOut.version = "3.1"
            mdOut.addDataTable("data_optics")
            mdOut.addLabels("data_optics", md.getLabels("data_optics"))
            mdOut.addData("data_optics", getattr(md, "data_optics"))
            particleTableName = "data_particles"
        else:
            particleTableName = "data_"

        mdOut.addDataTable(particleTableName)
        mdOut.addLabels(particleTableName, md.getLabels(particleTableName))
        mdOut.addData(particleTableName, new_particles)
        mdOut.write(args.o)

        print("New star file %s created. Have fun!" % args.o)


if __name__ == "__main__":
    SelCmmModStar().main()
