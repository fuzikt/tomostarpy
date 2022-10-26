#!/usr/bin/env python3

import os
import sys
from copy import deepcopy
from lib.metadata import MetaData
import argparse
from argparse import RawTextHelpFormatter


class FindInterstarDuplicates:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Remove duplicate particles appearing both in --i1 --i2 form --i2 and write result into --o.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i1', help="Input1 STAR filename.")
        add('--i2', help="Input2 STAR filename.")
        add('--o', help="Output STAR filename.")
        add('--tol', type=float, default="3.0",
            help="Tolerance in pixels of considering the coordinates of particles to be the same (Default: 3).")

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

        if not os.path.exists(args.i1):
            self.error("Input1 file '%s' not found." % args.i)

        if not os.path.exists(args.i2):
            self.error("Input2 file '%s' not found." % args.i)

    def get_particles(self, md, dataTableName):
        particles = []
        for particle in getattr(md, dataTableName):
            particles.append(particle)
        return particles

    def remove_interstar_duplicates(self, particles1, particles2, tolerance):
        outParticles = []
        foundDuplicate = False

        for particle2 in particles2:
            for particle1 in particles1:
                if (abs(particle2.rlnCoordinateX - particle1.rlnCoordinateX) <= tolerance) and (
                        abs(particle2.rlnCoordinateY - particle1.rlnCoordinateY) <= tolerance) and (
                        abs(particle2.rlnCoordinateZ - particle1.rlnCoordinateZ) <= tolerance):
                    foundDuplicate = True
                    break
            if not foundDuplicate:
                outParticles.append(particle2)
            foundDuplicate = False
        return outParticles

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        print("Selecting particles from star file...")

        tolerance = float(args.tol)

        md1 = MetaData(args.i1)
        md2 = MetaData(args.i2)

        dataTableName = "data_particles"

        particles1 = self.get_particles(md1, dataTableName)
        particles2 = self.get_particles(md2, dataTableName)

        print("%s particles in --i1" % str((len(particles1))))
        print("%s particles in --i2" % str((len(particles2))))

        print(
            "Finding duplicate particles between Input1 and Input2....")

        if md2.version == "3.1":
            mdOut = deepcopy(md2)
            mdOut.removeDataTable(dataTableName)
        else:
            mdOut = MetaData()
            dataTableName = "data_"

        mdOut.addDataTable(dataTableName)

        outParticles = self.remove_interstar_duplicates(particles1, particles2, tolerance)

        mdOut.addLabels(dataTableName, md2.getLabels(dataTableName))
        mdOut.addData(dataTableName, outParticles)

        print("%s particles were non-duplicates..." % str((len(outParticles))))
        mdOut.write(args.o)

        print("New star file %s created. Have fun!" % args.o)

if __name__ == "__main__":
    FindInterstarDuplicates().main()
