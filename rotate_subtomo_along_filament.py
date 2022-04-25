#!/usr/bin/env python3

import os
import sys
from lib.metadata import MetaData
import math
import argparse


class rotParticlesAlongFilament:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Set the Euler angles (tilt, psi) of the particles along a filament according to the angles of the between the neighboring particle centers.")
        add = self.parser.add_argument
        add('--i', help="Input STAR filename with particles.")
        add('--o', help="Output STAR filename.")

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

    def setEulerAngles(self, particles, apix):
        newParticles = []

        for i in range(len(particles)):
            newParticles.append(particles[i])
            if i > 0 and (i < len(particles) - 1):
                dx = (particles[i + 1].rlnCoordinateX * apix - particles[i + 1].rlnOriginXAngst) - (
                            particles[i].rlnCoordinateX * apix - particles[i -1].rlnOriginXAngst)
                dy = (particles[i + 1].rlnCoordinateY * apix - particles[i + 1].rlnOriginYAngst) - (
                            particles[i].rlnCoordinateY * apix - particles[i - 1].rlnOriginYAngst)
                dz = (particles[i + 1].rlnCoordinateZ * apix - particles[i + 1].rlnOriginZAngst) - (
                            particles[i].rlnCoordinateZ * apix - particles[i - 1].rlnOriginZAngst)
                vlength = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                dx /= vlength
                dy /= vlength
                dz /= vlength

                if vlength < 1000:
                    newParticles[i].rlnAngleRot = 0.0
                    newParticles[i].rlnAngleTilt = math.degrees(math.acos(dz))
                    newParticles[i].rlnAnglePsi = math.degrees(math.atan2(dx, dy)) + 90
                else:
                    newParticles[i].rlnAngleRot = 0.0
                    newParticles[i].rlnAngleTilt = newParticles[i-1].rlnAngleTilt
                    newParticles[i].rlnAnglePsi = newParticles[i-1].rlnAnglePsi

        # set the last the same orientation as the last befor last (as it does not have a following neighbor
        newParticles[-1].rlnAngleRot = 0.0
        newParticles[-1].rlnAngleTilt = newParticles[-2].rlnAngleTilt
        newParticles[-1].rlnAnglePsi = newParticles[-2].rlnAnglePsi

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

        print("Total %s particles in input star file." % str(len(particles)))

        new_particles.extend(self.setEulerAngles(particles, apix))

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
    rotParticlesAlongFilament().main()
