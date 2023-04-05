#!/usr/bin/env python3
import argparse
import sys
from lib.newstlist import *

def excludeFromTiltListFile(inputFile, outputFile, excludeList):

    with open(inputFile, 'r') as tiltListFile:
        tiltList =  tiltListFile.readlines()

    outputTiltList = []

    for pos, tiltItem in enumerate(tiltList, 1):
        if pos not in excludeList:
            outputTiltList.append(tiltItem)

    with open(outputFile, 'w') as outTiltListFile:
        for tiltItem in outputTiltList:
            outTiltListFile.write("%s" % tiltItem)

def main(inputFile, outputFile, excludeFile):
    with open(excludeFile, 'r') as tiltListFile:
        rangeList =  tiltListFile.read()

    if len(rangeList.strip()) == 0:
        print("No tilts listed to be excluded...")
        exit()

    tiltList = extractMembersOfRangeList(rangeList)

    excludeFromTiltListFile(inputFile, outputFile, tiltList)

    print("%s written..." % outputFile)
    print("All done. Have fun!")

if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Removes lines defined in --exclude_file from input text file. ")
    add = parser.add_argument
    add('--i', help="Input file to exclude lines from.")
    add('--o', help="Output file.")
    add('--exclude_file', type=str, default="",
        help="Textfile with comma separated list of excluded views (or range of views as in newstack).")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()

    main(args.i, args.o, args.exclude_file)