#!/usr/bin/env python3

import os
import sys
import math
import argparse
from lib.metadata import MetaData
import lib.matrix2 as matrix2


class rotateCoordinates:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Rotate the inplane coordinates X, Y of the subvolume.")
        add = self.parser.add_argument
        add('--i', help="Input star file.")
        add('--o', help="Output star file.")
        add('--ang', type=str, default="0.00",
            help="Rotation angle in degrees. (Default 0.0)")
        add('--xmax', type=float, default="4092",
            help="X-size of the unbinned tomogram. (Default 4092)")
        add('--ymax', type=float, default="5760",
            help="X-size of the unbinned tomogram. (Default  5760)")

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

    def rotatePrticleCoordXY(self, particles, angle, xmax, ymax):
        newParticles = []
        for particle in particles:
            x = particle.rlnCoordinateX - xmax / 2
            y = particle.rlnCoordinateY - ymax / 2
            rot_m = matrix2.matrix_from_angle(math.radians(angle))
            x = x * rot_m.m[0][0] + y * rot_m.m[0][1]
            y = x * rot_m.m[1][0] + y * rot_m.m[1][1]
            particle.rlnCoordinateX = x + xmax / 2
            particle.rlnCoordinateY = y + ymax / 2
            newParticles.append(particle)
        return newParticles

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        md = MetaData(args.i)

        new_particles = []

        print("Reading in input star file.....")

        particles = self.get_particles(md)

        print("Total %s particles in input star file." % str(len(particles)))

        new_particles.extend(self.rotatePrticleCoordXY(particles, args.ang, args.xmax, args.ymax))

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
    rotateCoordinates().main()
