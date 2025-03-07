#!/usr/bin/env python3

import os
import sys
from lib.metadata import MetaData
import argparse
from copy import deepcopy


class Relion5ToCryodrgnStar:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Set the Euler angles (tilt, psi) of the particles along a filament according to the angles of the between the neighboring particle centers.")
        add = self.parser.add_argument
        add('--i', help="Input particles STAR filename.")
        add('--t', help="Input tomograms STAR filename.")
        add('--project_dir', help="Relion project directory path.")
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

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        print("Reading in input star file.....")

        md_particles = MetaData(args.i)
        md_tomograms = MetaData(args.t)

        if md_particles.data_general[0].rlnTomoSubTomosAre2DStacks != True:
            print("ERROR: The input particles star file does not contain 2D stacks.")
            sys.exit(2)

        voltage = float(md_particles.data_optics[0].rlnVoltage)
        cs = float(md_particles.data_optics[0].rlnSphericalAberration)
        ac = float(md_particles.data_optics[0].rlnAmplitudeContrast)
        apix = float(md_particles.data_optics[0].rlnImagePixelSize)

        new_particles = []

        particles = self.get_particles(md_particles)

        tomoTiltSeriesStarFileDic = {}
        for tomoNameRecord in md_tomograms.data_global:
            tomoTiltSeriesStarFileDic[tomoNameRecord.rlnTomoName] = tomoNameRecord.rlnTomoTiltSeriesStarFile

        print("Converting particles.....")
        # progress bar initialization
        progress_step = max(int(len(particles) / 20), 1)
        i = 0

        for particle in particles:
            md_tiltSeries = MetaData(f"{args.project_dir}/{tomoTiltSeriesStarFileDic[particle.rlnTomoName]}")
            visibleFrames = particle.rlnTomoVisibleFrames.strip("[]").split(",")
            counter = 1
            for tiltNr, tilt in enumerate(getattr(md_tiltSeries, "data_" + particle.rlnTomoName), start=1):
                if visibleFrames[tiltNr - 1] != "0":
                    newParticle = deepcopy(particle)
                    newParticle.rlnImageName = f"{counter:2d}@{particle.rlnImageName}"
                    newParticle.rlnDefocusU = tilt.rlnDefocusU
                    newParticle.rlnDefocusV = tilt.rlnDefocusV
                    newParticle.rlnDefocusAngle = tilt.rlnDefocusAngle
                    if hasattr(tilt, "rlnCtfScalefactor"):
                        newParticle.rlnCtfScalefactor = tilt.rlnCtfScalefactor
                    else:
                        newParticle.rlnCtfScalefactor = 1.0
                    newParticle.rlnMicrographName = particle.rlnTomoName
                    newParticle.rlnGroupName = particle.rlnTomoName
                    newParticle.rlnVoltage = voltage
                    newParticle.rlnSphericalAberration = cs
                    newParticle.rlnAmplitudeContrast = ac
                    newParticle.rlnPhaseShift = 0.0
                    newParticle.rlnMagnification = 10000
                    newParticle.rlnDetectorPixelSize = apix
                    new_particles.append(newParticle)
                    counter += 1
            i += 1
             # a simple progress bar
            sys.stdout.write('\r')
            progress = int(i / progress_step)
            sys.stdout.write("[%-20s] %d%%" % ('=' * progress, 5 * progress))
            sys.stdout.flush()
        sys.stdout.write('\n')

        print("Total %s particles in input star file." % str(len(particles)))

        mdOut = MetaData()
        dataTableName = "data_"
        mdOut.addDataTable(dataTableName, True)
        mdOut.addLabels(dataTableName,
                        ["rlnMagnification", "rlnDetectorPixelSize", "rlnVoltage", "rlnSphericalAberration",
                         "rlnAmplitudeContrast", "rlnPhaseShift", "rlnDefocusU", "rlnDefocusV", "rlnDefocusAngle","rlnCtfScalefactor",
                         "rlnImageName", "rlnMicrographName", "rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi", "rlnGroupName"])
        mdOut.addData(dataTableName, new_particles)
        mdOut.write(args.o)

        print("New star file %s created. Have fun!" % args.o)


if __name__ == "__main__":
    Relion5ToCryodrgnStar().main()
