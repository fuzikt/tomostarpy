#!/usr/bin/env python3

import os
import sys
import math
from copy import deepcopy
from lib.metadata import MetaData
import numpy as np
import argparse


class sortByParticleDistance:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Select only particles that have certain amount of neighbors in a particular distance. Useful for filtering out template matched particles.")
        add = self.parser.add_argument
        add('--i', help="Input STAR file name with particles.")
        add('--o', help="Output STAR file name.")
        add('--dist', type=str, default="1.00",
            help="Distance in Angstroms that consider particles as neighbors. (Default 1.0)")
        add('--min_neigh', type=str, default="1",
            help="Minimum number of neighbors at --dist particle has to be kept in selection! (Default 1)")
        add('--min_corr', type=str, default="0",
            help="Minimum cross-correlation value of the neighbour to be considered as a true neighbour! (Default 0)")
        add('--lb_corr', type=str, default="rlnLCCmax",
            help="Label of the cross-correlation value in the star file! (Default rlnLCCmax)")
        add('--max_ang_dist', type=float, default="360",
            help="Maximum angular distance in degrees between neighbors to be considered as a true neighbour! (Default 360)")

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

    def angDistance(self, p1, p2):
        e1 = (p1.rlnAngleRot, p1.rlnAngleTilt, p1.rlnAnglePsi)
        e2 = (p2.rlnAngleRot, p2.rlnAngleTilt, p2.rlnAnglePsi)

        #        return self.angular_distance_matrix(e1, e2, degrees=True)
        return self.angular_distance_quat(e1, e2, degrees=True)

    def rot_z(self, angle):
        c, s = np.cos(angle), np.sin(angle)
        return np.array([
            [c, -s, 0],
            [s, c, 0],
            [0, 0, 1]
        ])

    def rot_y(self, angle):
        c, s = np.cos(angle), np.sin(angle)
        return np.array([
            [c, 0, s],
            [0, 1, 0],
            [-s, 0, c]
        ])

    def euler_zyz_to_matrix(self, rot, tilt, psi, degrees=True):
        if degrees:
            rot, tilt, psi = np.radians([rot, tilt, psi])
        return self.rot_z(rot) @ self.rot_y(tilt) @ self.rot_z(psi)

    def angular_distance_matrix(self, e1, e2, degrees=True):
        R1 = self.euler_zyz_to_matrix(*e1, degrees=degrees)
        R2 = self.euler_zyz_to_matrix(*e2, degrees=degrees)

        R_rel = R1.T @ R2
        trace = np.trace(R_rel)

        # Numerical safety
        cos_theta = (trace - 1) / 2
        cos_theta = np.clip(cos_theta, -1.0, 1.0)

        theta = np.arccos(cos_theta)
        return np.degrees(theta) if degrees else theta

    def angDistance(self,p1, p2):
        p1ang = (p1.rlnAngleRot, p1.rlnAngleTilt, p1.rlnAnglePsi)
        p2ang = (p2.rlnAngleRot, p2.rlnAngleTilt, p2.rlnAnglePsi)

        return self.angular_distance_matrix(p1ang, p2ang, degrees=True)

    def sortParticleByDistance(self, particles, apix, minNeighNr, maxDist, minXcorr, lbCorr, maxAngDist):
        newParticles = []
        counter = 0
        repeat = True
        while repeat:
            for particle in particles:
                for compParticle in particles:
                    if compParticle != particle:
                        if hasattr(particle, 'rlnCenteredCoordinateXAngst'):
                            partX = particle.rlnCenteredCoordinateXAngst
                            partY = particle.rlnCenteredCoordinateYAngst
                            partZ = particle.rlnCenteredCoordinateZAngst
                            compPartX = compParticle.rlnCenteredCoordinateXAngst
                            compPartY = compParticle.rlnCenteredCoordinateYAngst
                            compPartZ = compParticle.rlnCenteredCoordinateZAngst
                        else:
                            partX = particle.rlnCoordinateX * apix - particle.rlnOriginXAngst
                            partY = particle.rlnCoordinateY * apix - particle.rlnOriginYAngst
                            partZ = particle.rlnCoordinateZ * apix - particle.rlnOriginZAngst
                            compPartX = compParticle.rlnCoordinateX * apix - compParticle.rlnOriginXAngst
                            compPartY = compParticle.rlnCoordinateY * apix - compParticle.rlnOriginYAngst
                            compPartZ = compParticle.rlnCoordinateZ * apix - compParticle.rlnOriginZAngst
                        if (math.sqrt((partX-compPartX) ** 2 + (partY-compPartY) ** 2 + (partZ-compPartZ) ** 2) <= maxDist)\
                            and float(getattr(compParticle, lbCorr)) >= minXcorr and (abs(self.angDistance(particle,compParticle)) <= maxAngDist):
                            counter += 1
                            if counter >= minNeighNr:
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

        # get unbinned apix from star file
        if hasattr(md, "data_optics"):
            apix = float(md.data_optics[0].rlnTomoTiltSeriesPixelSize)
            print("Apix of star got form the first optic group. Apix = %f0.2" % apix)
        elif hasattr(md, "data_"):
            apix = float(md.data_[0].rlnDetectorPixelSize)
            print(
                "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % apix)
        elif hasattr(md, "data_particles"):
            apix = float(md.data_[0].rlnDetectorPixelSize)
            print(
                "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % apix)
        else:
            print("Could not get the apix of the particles from the star file. Define it using --star_apix parameter.")
            exit()

        print(f"Sorting {len(particles)} particles by distance.....")
        new_particles.extend(self.sortParticleByDistance(particles, apix, int(args.min_neigh), float(args.dist), float(args.min_corr), args.lb_corr, args.max_ang_dist))

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
