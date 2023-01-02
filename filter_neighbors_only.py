#!/usr/bin/env python3

import os
import sys
import math
from copy import deepcopy
from lib.metadata import MetaData
import argparse


class sortByParticleDistance:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Select only particles that have certain amount of neighbors in a particular distance. Useful for filtering out surface features.")
        add = self.parser.add_argument
        add('--i', help="Input STAR file name with particles.")
        add('--o', help="Output STAR file name.")
        add('--dist', type=str, default="1.00",
            help="Distance in agstroms that consider particles as neighbors. (Default 1.0)")
        add('--min_neigh', type=str, default="1",
            help="Minimum number of neighbors at --dist particle has to be kept in seection! (Default 1)")
        add('--min_corr', type=str, default="0",
            help="Minimum cross-correlation value of the neighbour to be considered as a true neighbour! (Default 0)")

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

    def sortParticleByDistance(self, particles, apix, minNeighNr, maxDist, minXcorr):
        newParticles = []
        counter = 0
        repeat = True

        while repeat:
            for particle in particles:
                for compParticle in particles:
                    if (math.sqrt((particle.rlnCoordinateX * apix - particle.rlnOriginXAngst - compParticle.rlnCoordinateX * apix + compParticle.rlnOriginXAngst) ** 2 + (
                            particle.rlnCoordinateY * apix - particle.rlnOriginYAngst - compParticle.rlnCoordinateY * apix + compParticle.rlnOriginYAngst) ** 2 + (
                            particle.rlnCoordinateZ * apix - particle.rlnOriginZAngst - compParticle.rlnCoordinateZ * apix + compParticle.rlnOriginZAngst) ** 2) <= maxDist) and compParticle.rlnCtfFigureOfMerit >= minXcorr:
                        counter += 1
                        if counter == minNeighNr:
                            counter = 0
                            newParticles.append(particle)
                            break
            if len(newParticles) < len(particles):
                repeat = True
                particles = deepcopy(newParticles)
                newParticles = []
            else:
                repeat = False

        print("Total %s particles were selected." % str(len(newParticles)))

        return newParticles

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        md = MetaData(args.i)

        new_particles = []

        print("Reading in input star file.....")

        particles = self.get_particles(md)

        optic_groups = []
        for optic_group in md.data_optics:
            optic_groups.append(optic_group)

        # get unbinned apix from star file
        apix = float(optic_groups[0].rlnTomoTiltSeriesPixelSize)

        new_particles.extend(self.sortParticleByDistance(particles, apix, int(args.min_neigh), float(args.dist), float(args.min_corr)))

        if md.version == "3.1":
            mdOut = md.clone()
            dataTableName = "data_particles"
            mdOut.removeDataTable(dataTableName)
        else:
            mdOut = MetaData()
            dataTableName = "data_"

        mdOut.addDataTable(dataTableName, md.isLoop(dataTableName))
        mdOut.addLabels(dataTableName, md.getLabels(dataTableName))
        mdOut.addData(dataTableName, new_particles)

        mdOut.write(args.o)

        print("New star file %s created. Have fun!" % args.o)


if __name__ == "__main__":
    sortByParticleDistance().main()
