#!/usr/bin/env python3

import os
import sys
from lib.metadata import MetaData
import copy
import argparse


class avgSmoothCoordinates:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Smooth particle positions by linear interpolation of the coordinates between the neighboring particles. Correction is stored in rlnOriginAngst.")
        add = self.parser.add_argument
        add('--i', help="Input STAR filename with particles.")
        add('--o', help="Output STAR filename.")
        add('--it', type=str, default="1",
            help="Number of smoothing iterations to perform (Default 1).")

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

    def get_particles(self, md):
        particles = []
        for particle in md:
            particles.append(particle)
        return particles

    def smoothParticleCoords(self, particles, apix, iterations):

        for iteration in range(iterations):
            if iteration > 0:
                particles = copy.deepcopy(newParticles)

            newParticles = []

            for i in range(len(particles)):
                newParticles.append(particles[i])
                if i > 0 and i < len(particles) - 1:
                    newParticles[i].rlnOriginXAngst = particles[i].rlnCoordinateX * apix - ((particles[
                                                                                                 i - 1].rlnCoordinateX * apix -
                                                                                             particles[
                                                                                                 i - 1].rlnOriginXAngst +
                                                                                             particles[
                                                                                                 i + 1].rlnCoordinateX * apix -
                                                                                             particles[
                                                                                                 i + 1].rlnOriginXAngst) / 2)
                    newParticles[i].rlnOriginYAngst = particles[i].rlnCoordinateY * apix - ((particles[
                                                                                                 i - 1].rlnCoordinateY * apix -
                                                                                             particles[
                                                                                                 i - 1].rlnOriginYAngst +
                                                                                             particles[
                                                                                                 i + 1].rlnCoordinateY * apix -
                                                                                             particles[
                                                                                                 i + 1].rlnOriginYAngst) / 2)
                    newParticles[i].rlnOriginZAngst = particles[i].rlnCoordinateZ * apix - ((particles[
                                                                                                 i - 1].rlnCoordinateZ * apix -
                                                                                             particles[
                                                                                                 i - 1].rlnOriginZAngst +
                                                                                             particles[
                                                                                                 i + 1].rlnCoordinateZ * apix -
                                                                                             particles[
                                                                                                 i + 1].rlnOriginZAngst) / 2)

        return newParticles

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        iterations = int(args.it)

        md = MetaData(args.i)

        new_particles = []

        print("Reading in input star file.....")

        optic_groups = []
        for optic_group in md.data_optics:
            optic_groups.append(optic_group)

        # get unbinned apix from star file
        apix = float(optic_groups[0].rlnTomoTiltSeriesPixelSize)

        particles = self.get_particles(md)

        print(
            "Total %s particles in input star file. \nSelecting one orientation per particle according to the greatest value of rlnMaxValueProbDistribution." % str(
                len(particles)))

        new_particles.extend(self.smoothParticleCoords(particles, apix, iterations))

        mdOut = MetaData()
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
    avgSmoothCoordinates().main()
