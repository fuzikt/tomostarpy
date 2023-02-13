#!/usr/bin/env python3
import os
import subprocess
import sys
import argparse

def runNewstackToExclude(inputFile, outputFile, excludeList):
    newstackcmd = "newstack"
    newStackCommand = [newstackcmd, "-in", inputFile, "-ou", outputFile, "-exc", ",".join([str(i) for i in excludeList]), "-fr"]
    #if self.versbosity > 1: print("Running: %s" % " ".join(newStackCommand))
    newstackProcess = subprocess.run(newStackCommand, stdout=subprocess.DEVNULL)
    if newstackProcess.returncode != 0:
        print("Something bad happened during running the newstack command.")
        print(newstackProcess.stderr)
        return False

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

def excludeFromMdocFile(inputFile, outputFile, excludeList):
    header = True
    headerLines = []
    zValueLines = []
    zValueBlocks = []
    outputZvalueBlocks = []

    with open(inputFile, 'r') as mdocFile:
        for mdocLine in mdocFile.readlines():
            if "ZValue" in mdocLine:
                header = False
                if len(zValueLines) > 0:
                    zValueBlocks.append(zValueLines)
                    zValueLines = []
            elif header:
                headerLines.append(mdocLine)
            else:
                zValueLines.append(mdocLine)

    for pos, zValueBlockItem in enumerate(zValueBlocks, 1):
        if pos not in excludeList:
            outputZvalueBlocks.append(zValueBlockItem)

    with open(outputFile, 'w') as outMdocFile:
        for headerLine in headerLines:
            outMdocFile.write("%s" % headerLine)
        for counter, zValueBlockItem in enumerate(outputZvalueBlocks):
            outMdocFile.write("[ZValue = %s ]\n" % counter)
            for zValueBlockLine in zValueBlockItem:
                outMdocFile.write("%s" % zValueBlockLine)


def extractMembersOfRangeList(rangeList):
    listItems = []
    for item in rangeList.split(","):
        if "-" in item:
            start = int(item.split("-")[0])
            stop = int(item.split("-")[1])+1
            listItems.extend(range(start,stop))
        else:
            listItems.append(int(item))
    return listItems

def main(inputPrefix, outputPrefix, excludeFile):
    with open(excludeFile, 'r') as tiltListFile:
        rangeList =  tiltListFile.read()

    if len(rangeList) == 0:
        print("No tilts listed to be excluded...")
        exit()

    tiltList = extractMembersOfRangeList(rangeList)

    if os.path.exists(inputPrefix + ".mrc"):
        print("Excluding tilts %s from %s..." % (rangeList, inputPrefix+".mrc"))
        runNewstackToExclude(inputPrefix+".mrc", outputPrefix+".mrc", tiltList)
    else:
        print("%s not found!!" % (inputPrefix+".rawtlt"))
        exit()

    if os.path.exists(inputPrefix+".rawtlt"):
        print("Excluding tilts from %s..." % (inputPrefix + ".rawtlt"))
        excludeFromTiltListFile(inputPrefix+".rawtlt", outputPrefix+".rawtlt", tiltList)
    else:
        print("%s not found skipping processing of it..." % (inputPrefix+".rawtlt"))

    if os.path.exists(inputPrefix + ".mdoc"):
        print("Excluding tilts from %s..." % (inputPrefix + ".mdoc"))
        excludeFromMdocFile(inputPrefix + ".mdoc", outputPrefix + ".mdoc", tiltList)
    else:
        print("%s not found skipping processing of it..." % (inputPrefix+".mdoc"))

    print("All done! Have fun!")

if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Removes tilts defined in --exclude_file from tilt series mrc-stack and corresponding *.rawtlt, *.tltdose, *.tltorder, *.mdoc. ")
    add = parser.add_argument
    add('--i', help="Input prefix mrc-stack.")
    add('--o', help="Output prefix.")
    add('--exclude_file', type=str, default="",
        help="Textfile with comma separated list of excluded views (or range of views as in newstack).")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()

    main(args.i, args.o, args.exclude_file)
