#!/usr/bin/env python3
import argparse
import math
from lib.matrix3 import *
from lib.vector3 import *
import random
from lib.metadata import MetaData
from lib.euler import *
from argparse import RawTextHelpFormatter


class subSubvolume:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Creates a star file with sub-subtomo coordinates similar to sub-particle approach in SPA.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input star file.")
        add('--o', help="Output star file.")
        add('--cmm',
              help="A CMM file defining the location(s) of the subparticle(s) "
                   "(use instead of --vector). Coordinates should be in Angstrom.")
        add('--vector', default="0,0,1", help="Vector defining the location of the subparticle. (default: 0,0,1)")
        add('--length',
              help="Alternative length of the vector. Use to adjust the "
                   "subparticle center (A). (default: length of the given vector)")
        add('--sym', default="C1", help="Symmetry of the particle. (default: C1)")
        add('--align_subparticles', action='store_true',
              help="Align subparticles to the standard orientation.")
        add('--randomize', action='store_true',
              help="Randomize the order of the symmetry matrices. "
                   "Useful for preventing preferred orientations.")
        add('--unique', type=float, default=-1,
            help="Keep only unique subparticles within angular distance "
                 "(useful to remove overlapping subparticles on symmetry axis).")
        add('--mindist', type=float, default=-1,
            help="Minimum distance between the subparticles in the image "
                 "(all overlapping ones will be discarded; pixels).")
        add('--side', type=float, default=-1,
            help="Keep only particles within specified angular distance from "
                 "side views (all others will be discarded; degrees).")
        add('--top', type=float, default=-1,
            help="Keep only particles within specified angular distance from "
                 "top views (all others will be discarded; degrees).")
        add('--library_path', type=str, default='',
              help="define LD_LIBRARY_PATH used. Default: empty")

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
        if args.o == "":
            self.error("No output prefix was specified.")

    def within_mindist(self, p1, p2, mindist):
        """ Returns True if two particles are closer to each other
        than the given distance in the projection. """

        x1 = p1.rlnCoordinateX
        y1 = p1.rlnCoordinateY
        z1 = p1.rlnCoordinateZ
        x2 = p2.rlnCoordinateX
        y2 = p2.rlnCoordinateY
        z2 = p2.rlnCoordinateZ
        distance_sqr = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
        mindist_sqr = mindist ** 2

        return distance_sqr < mindist_sqr

    def within_unique(self, p1, p2, unique):
        """ Returns True if two particles are closer to each other
        than the given angular distance. """

        v1 = vector_from_two_eulers(radians(p1.rlnAnglePsi), radians(p1.rlnAngleTilt))
        v2 = vector_from_two_eulers(radians(p2.rlnAnglePsi), radians(p2.rlnAngleTilt))

        dp = dot_product(v1, v2) / (v1.length() * v2.length())

        if dp < -1:
            dp = -1.000

        if dp > 1:
            dp = 1.000

        angle = math.acos(dp)

        return angle <= math.radians(unique)

    def filter_unique(self, subparticles, subpart, unique):
        """ Return True if subpart is not close to any other subparticle
            by unique (angular distance).
            For this function we assume that subpart is not contained
            inside."""
        for sp in subparticles:
            if self.within_unique(sp, subpart, unique):
                return False
        return True

    def filter_mindist(self, subparticles, subpart, mindist):
        """ Return True if subpart is not close to any other subparticle
        by mindist. """
        for sp in subparticles:
            # TODO redefien for unique subtomo
            if (sp.rlnImageName[:6] != subpart.rlnImageName[:6] and
                    self.within_mindist(sp, subpart, mindist)):
                return False
        return True

    def filter_side(self, subpart, side):
        return (abs(abs(radians(subpart.rlnAngleTilt)) - radians(90)) < side)

    def filter_top(self, subpart, top):
        return (abs(abs(radians(subpart.rlnAngleTilt)) - radians(180)) < top)

    def filter_subparticles(self, subparticles, filters):
        return [sp for sp in subparticles
                if all(f(subparticles, sp) for f in filters)]

    def create_subparticles(self, particle, symmetry_matrices, subparticle_vector_list,
                            randomize, unique, align_subparticles, filters, ang_pix):
        """ Obtain all subparticles from a given particle and set
        the properties of each such subparticle. """

        # Euler angles that take particle to the orientation of the model

        rot = math.radians(particle.rlnAngleRot)
        tilt = math.radians(particle.rlnAngleTilt)
        psi = math.radians(particle.rlnAnglePsi)

        matrix_particle = matrix_from_euler(rot, tilt, psi)

        subparticles = []
        # subpart_id = 1
        # subparticles_total += 1

        symmetry_matrix_ids = range(1, len(symmetry_matrices) + 1)

        if randomize:
            # randomize the order of symmetry matrices, prevents preferred views
            random.shuffle(symmetry_matrix_ids)

        for subparticle_vector in subparticle_vector_list:
            matrix_from_subparticle_vector = subparticle_vector.matrix()

            for symmetry_matrix_id in symmetry_matrix_ids:
                # symmetry_matrix_id can be later written out to find out
                # which symmetry matrix created this subparticle
                symmetry_matrix = symmetry_matrices[symmetry_matrix_id - 1]

                subpart = particle.clone()

                m = matrix_multiply(matrix_particle, (matrix_multiply(matrix_transpose(symmetry_matrix),
                                                                      matrix_transpose(
                                                                          matrix_from_subparticle_vector))))

                if align_subparticles:
                    rotNew, tiltNew, psiNew = euler_from_matrix(m)
                else:
                    m2 = matrix_multiply(matrix_particle, matrix_transpose(symmetry_matrix))
                    rotNew, tiltNew, psiNew = euler_from_matrix(m2)

                # save Euler angles that take the model to the orientation of the subparticle

                subpart.rlnAngleRot = math.degrees(rotNew)
                subpart.rlnAngleTilt = math.degrees(tiltNew)
                subpart.rlnAnglePsi = math.degrees(psiNew)

                # subparticle origin
                d = subparticle_vector.distance()

                x = -m.m[0][2] * d + particle.rlnOriginXAngst / ang_pix
                y = -m.m[1][2] * d + particle.rlnOriginYAngst / ang_pix
                z = -m.m[2][2] * d + particle.rlnOriginZAngst / ang_pix

                # save the subparticle coordinates (integer part) relative to the
                # user given image size and as a small shift in the origin (decimal part)
                x_d, x_i = math.modf(x)
                y_d, y_i = math.modf(y)
                z_d, z_i = math.modf(z)

                subpart.rlnCoordinateX = particle.rlnCoordinateX - x_i
                subpart.rlnCoordinateY = particle.rlnCoordinateY - y_i
                subpart.rlnCoordinateZ = particle.rlnCoordinateZ - z_i

                subpart.rlnOriginXAngst = x_d * ang_pix
                subpart.rlnOriginYAngst = y_d * ang_pix
                subpart.rlnOriginZAngst = z_d * ang_pix

                overlaps = (unique >= 0 and not self.filter_unique(subparticles, subpart, unique))

                if not overlaps:
                # TODO possibility to filter particles outside of tomograms (then indent the next 3 lines behind the if statement
                # if not ((subpart.rlnCoordinateX<0) or (subpart.rlnCoordinateY<0) or (subpart.rlnCoordinateX>part_image_sizeX) or (subpart.rlnCoordinateY>part_image_sizeY)):
                    subparticles.append(subpart)
                # subpart_id += 1
                # subparticles_total += 1

        if filters:
            subparticles = self.filter_subparticles(subparticles, filters)

        return subparticles

    def load_vectors(self, cmm_file, vectors_str, distances_str, angpix):
        """ Load subparticle vectors either from Chimera CMM file or from
        a vectors string. Distances can also be specified for each vector
        in the distances_str. """

        if cmm_file:
            subparticle_vector_list = vectors_from_cmm(cmm_file, angpix)
        else:
            subparticle_vector_list = vectors_from_string(vectors_str)

        if distances_str:
            # Change distances from A to pixel units
            subparticle_distances = [float(x) / angpix for x in
                                     distances_str.split(',')]

            if len(subparticle_distances) != len(subparticle_vector_list):
                raise Exception("Error: The number of distances does not match "
                                "the number of vectors!")

            for vector, distance in zip(subparticle_vector_list,
                                         subparticle_distances):
                if distance > 0:
                    vector.set_distance(distance)
                else:
                    vector.compute_distance()
        else:
            for vector in subparticle_vector_list:
                vector.compute_distance()

        print("Using vectors:")

        for subparticle_vector in subparticle_vector_list:
            print("Vector: ")
            subparticle_vector.normalize()
            subparticle_vector.compute_matrix()
            subparticle_vector.print_vector()
            print("")
            print("Length: %.2f pixels" % subparticle_vector.distance())
        print("")

        return subparticle_vector_list

    def load_filters(self, side, top, mindist):
        """ Create some filters depending on the conditions imposed by the user.
        Each filter will return True if the subparticle will be kept in the
        subparticles list.
        """
        filters = []

        if side > 0:
            filters.append(lambda x, y: self.filter_side(y, side))

        if top > 0:
            filters.append(lambda x, y: self.filter_top(y, top))

        if mindist > 0:
            filters.append(lambda x, y: self.filter_mindist(x, y, mindist))

        return filters

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        md = MetaData(args.i)
        optic_groups = []
        for optic_group in md.data_optics:
            optic_groups.append(optic_group)

        # get unbinned apix from star file
        apix = float(optic_groups[0].rlnTomoTiltSeriesPixelSize)

        subparticle_vector_list = self.load_vectors(args.cmm, args.vector,
                                               args.length, apix)

        print("Creating subparticles...")

        # Generate symmetry matrices with Relion convention
        symmetry_matrices = matrix_from_symmetry(args.sym, args.library_path)

        # Define some conditions to filter subparticles
        filters = self.load_filters(radians(args.side), radians(args.top), args.mindist)

        mdOut = MetaData()
        mdOut.version = "3.1"
        mdOut.addDataTable("data_optics", True)
        mdOut.addLabels("data_optics", md.getLabels("data_optics"))
        mdOut.addData("data_optics", getattr(md, "data_optics"))
        particleTableName = "data_particles"
        mdOut.addDataTable(particleTableName, True)
        md.removeLabels(particleTableName, 'rlnTomoParticleName')
        mdOut.addLabels(particleTableName, md.getLabels(particleTableName))

        for particle in md:
            subparticles = self.create_subparticles(particle,
                                               symmetry_matrices,
                                               subparticle_vector_list,
                                               args.randomize,
                                               args.unique,
                                               args.align_subparticles,
                                               filters,
                                               apix)

            mdOut.addData(particleTableName, subparticles)

        mdOut.write(args.o)

if __name__ == "__main__":
    subSubvolume().main()
