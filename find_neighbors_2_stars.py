#!/usr/bin/env python3

import os
import sys
import math
from lib.metadata import MetaData
from lib.euler import euler_from_vector, euler_from_matrix
from lib.vector3 import Vector3, vector_from_two_eulers, dot_product
import argparse


class sortByParticleDistance:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Select particles from --i2 that are in max distance --dist from particles in --i1. Useful for filtering out subparticles around particle centers.")
        add = self.parser.add_argument
        add('--i1', help="Input STAR file name with particles. Defines the particles around which the neighbors from --i2 are matched.")
        add('--i2', help="Input STAR file name with particles.")
        add('--o', help="Output STAR file name.")
        add('--star_apix1', type=float, default="0.0",
            help="Apix of the coordinates in the --i1 star file. Autodetected form star file if set to 0. (Default: 0 - autodetect)")
        add('--star_apix2', type=float, default="0.0",
            help="Apix of the coordinates in the --i2 star file. Autodetected form star file if set to 0. (Default: 0 - autodetect)")
        add('--dist', type=str, default="1.00",
            help="Distance in Angstroms that consider particles as neighbors. (Default 1.0)")
        add('--tilt_diff', type=str, default="360.00",
            help="Max +-titlAngle difference in degrees between the vector from the particle center in --i1 and particle rlnAngleTilt from --i2  (Default 360)")


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
            self.error("Input file '%s' not found."
                       % args.i1)

        if not os.path.exists(args.i2):
            self.error("Input file '%s' not found."
                       % args.i2)

    def get_particles(self, md):
        particles = []
        for particle in md:
            particles.append(particle)
        return particles

    def findNeigborsByDistance(self, particles1, particles2, apix1, apix2, maxDist, tilt_diff):
        newParticles = []

        counter = 1
        for particle in particles1:
            if hasattr(particle, "rlnCenteredCoordinateXAngst"):
                # switch to Relion5 coordinate definition for all coordinates
                particleCoordXAngst = particle.rlnCenteredCoordinateXAngst
                particleCoordYAngst = particle.rlnCenteredCoordinateYAngst
                particleCoordZAngst = particle.rlnCenteredCoordinateZAngst
            else:
                particleCoordXAngst = particle.rlnCoordinateX * apix1
                particleCoordYAngst = particle.rlnCoordinateY * apix1
                particleCoordZAngst = particle.rlnCoordinateZ * apix1

            try:
                particleCenterX = particleCoordXAngst + particle.rlnOriginXAngst
            except:
                particleCenterX = particleCoordXAngst
            try:
                particleCenterY = particleCoordYAngst + particle.rlnOriginYAngst
            except:
                particleCenterY = particleCoordYAngst
            try:
                particleCenterZ = particleCoordZAngst + particle.rlnOriginZAngst
            except:
                particleCenterZ = particleCoordZAngst

            indicesToRemove = []

            for compParticle in particles2:
                if hasattr(compParticle, "rlnCenteredCoordinateXAngst"):
                    #switch to Relion5 coordinate definition for all coordinates
                    compParticleCoordXAngst = compParticle.rlnCenteredCoordinateXAngst
                    compParticleCoordYAngst = compParticle.rlnCenteredCoordinateYAngst
                    compParticleCoordZAngst = compParticle.rlnCenteredCoordinateZAngst
                else:
                    compParticleCoordXAngst = compParticle.rlnCoordinateX * apix2
                    compParticleCoordYAngst = compParticle.rlnCoordinateY * apix2
                    compParticleCoordZAngst = compParticle.rlnCoordinateZ * apix2
                try:
                    compParticleX = compParticleCoordXAngst + compParticle.rlnOriginXAngst
                except:
                    compParticleX = compParticleCoordXAngst
                try:
                    compParticleY = compParticleCoordYAngst + compParticle.rlnOriginYAngst
                except:
                    compParticleY = compParticleCoordYAngst
                try:
                    compParticleZ = compParticleCoordZAngst + compParticle.rlnOriginZAngst
                except:
                    compParticleZ = compParticleCoordZAngst

                try:
                    compParticleTiltAngle = compParticle.rlnAngleTilt
                except:
                    compParticleTiltAngle = 0.0

                try:
                    compParticlePsiAngle = compParticle.rlnAnglePsi
                except:
                    compParticlePsiAngle = 0.0

                try:
                    compParticleRotAngle = compParticle.rlnAngleRot
                except:
                    compParticleRotAngle = 0.0

                vector = Vector3(
                    [compParticleX - particleCenterX, compParticleY - particleCenterY, compParticleZ - particleCenterZ])

                if (vector.length() <= maxDist):
                            rot, tilt, psi = euler_from_vector(vector)
                            v2 = vector_from_two_eulers(math.radians(compParticleRotAngle), math.radians(compParticlePsiAngle))
                            v2.scale(-1)
                            compParticleRotAngle, compParticleTiltAngle, compParticlePsiAngle = euler_from_vector(v2)

                            dp = math.sin(tilt)*math.sin(compParticleTiltAngle)+math.cos(tilt)*math.cos(compParticleTiltAngle)*math.cos(rot-compParticleRotAngle)

                            angle = math.acos(dp)

                            if angle <= math.radians(tilt_diff):
                                newParticles.append(compParticle)
                                indicesToRemove.append(particles2.index(compParticle))
                                #print(compParticlePsiAngle)
                                counter += 1

            for i in sorted(indicesToRemove, reverse=True):
                del particles2[i]

        print("Total %s particles were selected." % str(len(newParticles)))

        return newParticles

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()

        self.validate(args)

        md1 = MetaData(args.i1)
        md2 = MetaData(args.i2)

        new_particles = []

        print("Reading in input star file.....")

        particles1 = self.get_particles(md1)
        particles2 = self.get_particles(md2)

        if args.star_apix1 == 0:
            if hasattr(md1, "data_optics"):
                apix1 = float(md1.data_optics[0].rlnTomoTiltSeriesPixelSize)
                print("Apix of star got form the first optic group. Apix = %f0.2" % apix1)
            elif hasattr(md1, "data_"):
                apix1 = float(md1.data_[0].rlnDetectorPixelSize)
                print(
                    "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % apix1)
            elif hasattr(md1, "data_particles"):
                apix1 = float(md1.data_[0].rlnDetectorPixelSize)
                print(
                    "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % apix1)
            else:
                print(
                    "Could not get the apix of the particles from the star file. Define it using --star_apix parameter.")
                exit()
        else:
            apix1 = args.star_apix1

        if args.star_apix2 == 0:
            if hasattr(md2, "data_optics"):
                apix2 = float(md2.data_optics[0].rlnTomoTiltSeriesPixelSize)
                print("Apix of star got form the first optic group. Apix = %f0.2" % apix2)
            elif hasattr(md2, "data_"):
                apix2 = float(md2.data_[0].rlnDetectorPixelSize)
                print(
                    "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % apix2)
            elif hasattr(md2, "data_particles"):
                apix2 = float(md2.data_[0].rlnDetectorPixelSize)
                print(
                    "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % apix2)
            else:
                print(
                    "Could not get the apix of the particles from the star file. Define it using --star_apix parameter.")
                exit()
        else:
            apix2 = args.star_apix2

        new_particles.extend(self.findNeigborsByDistance(particles1, particles2, apix1, apix2, float(args.dist), float(args.tilt_diff)))

        if md2.version == "3.1":
            mdOut = md2.clone()
            dataTableName = "data_particles"
            mdOut.removeDataTable(dataTableName)
        else:
            mdOut = MetaData()
            dataTableName = "data_"

        mdOut.addDataTable(dataTableName, md2.isLoop(dataTableName))
        mdOut.addLabels(dataTableName, md2.getLabels(dataTableName))
        mdOut.addData(dataTableName, new_particles)

        mdOut.write(args.o)

        print("New star file %s created. Have fun!" % args.o)


if __name__ == "__main__":
    sortByParticleDistance().main()
