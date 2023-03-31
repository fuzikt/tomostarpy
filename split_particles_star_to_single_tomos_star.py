#!/usr/bin/env python3
import argparse
import sys
import os
from lib.metadata import MetaData

class SplitParticleStarPerTomos:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
        description="Split a particle star file into separate files per tomoName containing the corresponding particles.")
        add = self.parser.add_argument
        add('--i', help="Input star file.")
        add('--o', default="", help="Output directory name.")

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

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        if not os.path.exists(args.o):
            os.makedirs(args.o)

        mdparticles = MetaData(args.i)

        particles = self.get_particles(mdparticles)

        tomoNames = []
        # get all cluster identifiers
        for particle in particles:
            if particle.rlnTomoName not in tomoNames:
                tomoNames.append(particle.rlnTomoName)

        # sort and unique
        tomoNames = sorted(set(tomoNames))

        for tomoName in tomoNames:
            perTomoParticles = []
            for particle in particles:
                if tomoName == particle.rlnTomoName:
                    perTomoParticles.append(particle)
            #write star
            mdOut = mdparticles.clone()
            mdOut.removeDataTable("data_particles")

            mdOut.addDataTable("data_particles", mdparticles.isLoop("data_particles"))
            mdOut.addLabels("data_particles", mdparticles.getLabels("data_particles"))
            mdOut.addData("data_particles", perTomoParticles)
            mdOut.write(args.o+"/"+tomoName+".star")



if __name__ == "__main__":
    SplitParticleStarPerTomos().main()