#!/usr/bin/env python3
import os
import sys
import argparse
from lib.placeback_subvolume_cyt import *
from lib.placeback_subvolume_cuda import *
from argparse import RawTextHelpFormatter


class placebackSubvolume:
    def define_parser(self):
        self.parser = argparse.ArgumentParser(
            description="Performs placeback of a subvolume into full tomogram volume according to the coordinates and euler angles defined in a star file.",
            formatter_class=RawTextHelpFormatter)
        add = self.parser.add_argument
        add('--i', help="Input star file file.")
        add('--isub', type=str, default="", help="Input file of sub-volume to be placed.")
        add('--itomo', type=str, default="",
            help="Input file of the stencil tomogram. Size and origin is taken from this file and applied on the output.")
        add('--o', help="Output prefix")
        add('--cmm', dest='cmm', action='store_true', default=False,
            help="Create Chimera cmm file with the coordinates of the placed sub-volumes.")
        add('--tomo_name', type=str, default="",
            help="Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles used in place-back")
        add('--star_apix', type=float, default="0.0",
            help="Apix of the coordinates in the star file. Autodetected form star file if set to 0. (Default: 0 - autodetect)")
        add('--bin', type=str, default="0",
            help="Binning factor of the stencil tomogram. If not provided calculated from MRC header and star file apix")
        add('--no_partial', dest='no_partial', action='store_true', default=False,
            help="Do not place partial volumes. If the subvolume is partially out of the output volume, it is not placed at all. This avoids \"half-cut\" sub-volumes in the output.")
        add('--recenter', dest='recenter', action='store_true', default=False,
            help="Recenter the particles by subtracting rlnOriginX/Y/ZAngst from X/Y/Z coordinates.")
        add('--color_lb', type=str, default="",
            help="Label from the star file that will be used for rainbow coloring of the cmm markers.")
        add('--color_map', dest='color_map', action='store_true', default=False,
            help="Create map with coloring values stored as pixel value.")
        add('--radial_color', dest='radial_color', action='store_true', default=False,
            help="Create map with coloring values storing the value of radial distance form the center of the box. Useful for radial coloring.")
        add('--color_map_threshold', type=str, default="0.01",
            help="Threshold value at which the contour of the --isub in the output color map should contain values")
        add('--color_map_extend', type=str, default="2",
            help="Extend the border around the threshold meeting value by this amount of pixels.")
        add('--xtilt', type=str, default="0.00",
            help="Tomo X axis tilt in degrees. (Default: 0)")
        add('--ytilt', type=str, default="0.00",
            help="Tomo Y axis tilt in degrees. (Default: 0)")
        add('--gpu', dest='gpu', action='store_true', default=False,
            help="Use GPU acceleration.")

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

        if args.isub == "":
            self.error("No input subvolume to place was given. Please, use the --isub to specify one")

        if args.itomo == "":
            self.error("No input stencil tomogram was given. Please, use the --itomo to specify one.")

        if args.o == "":
            self.error("No output prefix was specified.")

        if float(args.color_map_extend) < 1:
            self.error("--color_map_extend has to be a positive number.")

    def main(self):
        self.define_parser()
        args = self.parser.parse_args()
        self.validate(args)

        inputStarFile = args.i
        inputVolumeToPlace = args.isub
        outputMapStencil = args.itomo
        outputPrefix = args.o
        outputCmm = args.cmm
        tomoName = args.tomo_name
        starApix = float(args.star_apix)
        binning = float(args.bin)
        if args.no_partial == False:
            placePartialVolumes = True
        else:
            placePartialVolumes = False

        recenter = args.recenter

        coloringLabel = args.color_lb

        outputColorMap = args.color_map

        colorMapTreshold = float(args.color_map_threshold)

        colorMapExtend = int(float(args.color_map_extend))

        radial_color = args.radial_color

        Xtilt = float(args.xtilt)
        Ytilt = float(args.ytilt)

        if args.gpu:
            placeSubvolumes_gpu(inputStarFile, inputVolumeToPlace, outputMapStencil, outputPrefix, outputCmm, tomoName,
                            binning,
                            placePartialVolumes, recenter, coloringLabel, outputColorMap, colorMapTreshold, colorMapExtend,
                            radial_color, Xtilt, Ytilt)
        else:
            placeSubvolumes(inputStarFile, inputVolumeToPlace, outputMapStencil, outputPrefix, outputCmm, tomoName, starApix,
                            binning, placePartialVolumes, recenter, coloringLabel, outputColorMap, colorMapTreshold,
                            colorMapExtend, radial_color, Xtilt, Ytilt)


if __name__ == "__main__":
    placebackSubvolume().main()
