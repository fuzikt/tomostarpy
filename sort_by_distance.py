#!/usr/bin/env python3

import os
import sys
from lib.metadata import MetaData
import argparse


class sortByParticleDistance:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Sort particles in star file to order based on the minimal distance between consecutive particles. In general, this should form a filament path.")
        add = self.parser.add_argument
        add('--i', help="Input STAR file name with particles.")
        add('--o', help="Output STAR file name.")

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

    def sortParticleByDistance(self, particles, apix):
        newParticles = []
        distxyz = 0
        counter = 0

        # find center of all points
        sumX = 0
        sumY = 0
        sumZ = 0

        for particle in particles:
            sumX += particle.rlnCoordinateX * apix - particle.rlnOriginXAngst
            sumY += particle.rlnCoordinateX * apix - particle.rlnOriginXAngst
            sumZ += particle.rlnCoordinateX * apix - particle.rlnOriginXAngst

        avgX = sumX / len(particles)
        avgY = sumY / len(particles)
        avgZ = sumZ / len(particles)

        # find the particle most distant form the center of all points and assume that this is the begining of the filament
        for particle in particles:
            if distxyz < (particle.rlnCoordinateX * apix - particle.rlnOriginXAngst - avgX) ** 2 + (
                    particle.rlnCoordinateY * apix - particle.rlnOriginYAngst - avgY) ** 2 + (
                    particle.rlnCoordinateZ * apix - particle.rlnOriginZAngst - avgZ) ** 2:
                distxyz = (particle.rlnCoordinateX * apix - particle.rlnOriginXAngst - avgX) ** 2 + (
                            particle.rlnCoordinateY * apix - particle.rlnOriginYAngst - avgY) ** 2 + (
                                      particle.rlnCoordinateZ * apix - particle.rlnOriginZAngst - avgZ) ** 2
                maxDistParticleIndex = counter
            counter += 1
        newParticles.append(particles.pop(maxDistParticleIndex))

        while len(particles) > 0:
            lastParticleInRow = newParticles[-1]
            lastMinDist = ((lastParticleInRow.rlnCoordinateX * apix - lastParticleInRow.rlnOriginXAngst - (
                        particles[0].rlnCoordinateX * apix - particles[0].rlnOriginXAngst)) ** 2 + (
                                       lastParticleInRow.rlnCoordinateY * apix - lastParticleInRow.rlnOriginYAngst - (
                                           particles[0].rlnCoordinateY * apix - particles[0].rlnOriginYAngst)) ** 2 + (
                                       lastParticleInRow.rlnCoordinateZ * apix - lastParticleInRow.rlnOriginZAngst - (
                                           particles[0].rlnCoordinateZ * apix - particles[0].rlnOriginZAngst)) ** 2)

            minDistIndex = 0

            if len(particles) != 0:
                particlePos = 0
                while particlePos < len(particles):
                    if ((lastParticleInRow.rlnCoordinateX * apix - lastParticleInRow.rlnOriginXAngst - (
                            particles[particlePos].rlnCoordinateX * apix - particles[
                        particlePos].rlnOriginXAngst)) ** 2 + (
                                lastParticleInRow.rlnCoordinateY * apix - lastParticleInRow.rlnOriginYAngst - (
                                particles[particlePos].rlnCoordinateY * apix - particles[
                            particlePos].rlnOriginYAngst)) ** 2 + (
                                lastParticleInRow.rlnCoordinateZ * apix - lastParticleInRow.rlnOriginZAngst - (
                                particles[particlePos].rlnCoordinateZ * apix - particles[
                            particlePos].rlnOriginZAngst)) ** 2) < lastMinDist:
                        lastMinDist = ((lastParticleInRow.rlnCoordinateX * apix - lastParticleInRow.rlnOriginXAngst - (
                                    particles[particlePos].rlnCoordinateX * apix - particles[
                                particlePos].rlnOriginXAngst)) ** 2 + (
                                                   lastParticleInRow.rlnCoordinateY * apix - lastParticleInRow.rlnOriginYAngst - (
                                                       particles[particlePos].rlnCoordinateY * apix - particles[
                                                   particlePos].rlnOriginYAngst)) ** 2 + (
                                                   lastParticleInRow.rlnCoordinateZ * apix - lastParticleInRow.rlnOriginZAngst - (
                                                       particles[particlePos].rlnCoordinateZ * apix - particles[
                                                   particlePos].rlnOriginZAngst)) ** 2)
                        minDistIndex = particlePos
                    particlePos += 1
            newParticles.append(particles.pop(minDistIndex))

        print("Total %s particles were sorted." % str(len(newParticles)))

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

        new_particles.extend(self.sortParticleByDistance(particles, apix))

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
