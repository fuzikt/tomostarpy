#!/usr/bin/env python3

import os
import sys
import argparse

class MergeMdoc:

    def __init__(self, tilts_dir, out_dir, startTilt=0, startTiltSeriesNr=1, tomo_prefix='tomo', mdoc_suffix=".tif.mdoc", pre_dose=0.0, dose=0.0, verbosity = 1):
        self.tilts_dir = tilts_dir
        self.out_dir = out_dir
        self.startTilt = startTilt
        self.lastTiltSeriesNr = startTiltSeriesNr
        self.tomo_prefix = tomo_prefix
        self.mdoc_suffix = mdoc_suffix
        self.pre_dose = pre_dose
        self.dose = dose
        self.versbosity = verbosity

    def processMdocFiles(self, mdocTiltList, tomoNr):
        #generate rawtlt file
        rawTltName = self.out_dir + "/" + self.tomo_prefix + str(tomoNr) +"/"+ self.tomo_prefix + str(tomoNr)+".rawtlt"
        with open(rawTltName, 'w') as rawTltFile:
            for mdocName in mdocTiltList:
                if not os.path.isfile(mdocName):
                    print("ERROR: File %s does not exist. Skipping mdoc processing." % mdocName)
                    return
                with open(mdocName, 'r') as mdocFile:
                    for mdocLine in mdocFile.readlines():
                        if "TiltAngle" in mdocLine:
                            tiltAngle = float(mdocLine.split()[2])
                rawTltFile.write("%0.2f\n" % tiltAngle)
        if self.versbosity > 0: print(" -> %s created." % rawTltName)

        #generate joined mdoc file
        tomoMdocName = self.out_dir + "/" + self.tomo_prefix + str(tomoNr) +"/"+ self.tomo_prefix + str(tomoNr) + ".mdoc"
        with open(tomoMdocName, 'w') as tomoMdocFile:
            firstMdoc = True
            frameSetFound = False
            tomoCounter = 0
            for mdocName in mdocTiltList:
                with open(mdocName, 'r') as mdocFile:
                    for mdocLine in mdocFile.readlines():
                        if "FrameSet" in mdocLine:
                            frameSetFound = True
                            firstMdoc = False
                            tomoMdocFile.write("\n[ZValue = %s ]\n" % tomoCounter)
                            continue
                        if not frameSetFound and firstMdoc:
                            tomoMdocFile.write(mdocLine)
                            continue
                        if frameSetFound and mdocLine.strip():
                            tomoMdocFile.write(mdocLine)
                frameSetFound = False
                tomoCounter += 1
        if self.versbosity > 0: print(" -> %s created." % tomoMdocName)

        #generate tiltOrderList + dosePerTiltList
        mdocTiltList.sort(key=lambda x: x.split("_")[1])
        tiltCounter = 0
        tiltDose = self.pre_dose
        tiltDoseList = []
        tltOrderName = self.out_dir + "/" + self.tomo_prefix + str(tomoNr) + "/" + self.tomo_prefix + str(
            tomoNr) + ".tltorder"
        with open(tltOrderName, 'w') as tltOrderFile:
            for mdocName in mdocTiltList:
                tiltCounter += 1
                tiltDose += self.dose
                with open(mdocName, 'r') as mdocFile:
                    for mdocLine in mdocFile.readlines():
                        if "TiltAngle" in mdocLine:
                            tiltAngle = float(mdocLine.split()[2])
                            break
                    tltOrderFile.write("%s, %0.2f\n" % (tiltCounter, tiltAngle))
                    tiltDoseList.append("%0.2f %0.2f\n" % (tiltAngle, tiltDose))
        if self.versbosity > 0: print(" -> %s created." % tltOrderName)

        tltDoseName = self.out_dir + "/" + self.tomo_prefix + str(tomoNr) + "/" + self.tomo_prefix + str(
            tomoNr) + ".tltdose"
        tiltDoseList.sort(key=lambda x: float(x.split(" ")[0]))
        with open(tltDoseName, 'w') as tltDoseFile:
            tltDoseFile.write("".join(tiltDoseList))
        if self.versbosity > 0: print(" -> %s created." % tltDoseName)


    def writePerPointTitlSeriesListFile(self, perPointTiltSeries):
        os.chdir(self.tilts_dir)

        for i in range(len(perPointTiltSeries)):
            if self.versbosity > 0 : print("Processing %s%s..." % (self.tomo_prefix, self.lastTiltSeriesNr + i))

            perPointTiltSeries[i].sort(key=lambda x: float(x.split("_")[2].replace(self.mdoc_suffix, "")))

            # create tomo dir
            if not os.path.exists("%s/%s%s" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i)):
                if self.versbosity > 1 : print("Creating %s directory in output path..." % self.tomo_prefix)
                try:
                    os.makedirs("%s/%s%s" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i))
                except:
                    print("ERROR: Could not create output directory: %s/%s%s" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i))
                    return False

            # Process mdoc files
            self.processMdocFiles(perPointTiltSeries[i], self.lastTiltSeriesNr + i)

        self.lastTiltSeriesNr += (i + 1)
        return True

    def processFilesInList(self, tilt_files):
        perPointTiltSeries = []
        lookForNewStartTilt = True

        for fileName in tilt_files:
            if float(fileName.split("_")[2].replace(self.mdoc_suffix, "")) == self.startTilt:
                if lookForNewStartTilt:
                    if len(perPointTiltSeries) > 0:
                        if self.writePerPointTitlSeriesListFile(perPointTiltSeries):
                            # reset to 0, to let the "watchdog" know that files were removed
                            lastAmountOfFiles = 0
                        else:
                            # something wrong happened during the tilt-series file creation process
                            break
                    perPointTiltSeries = []
                    perPointTiltSeries.append([fileName])
                    lookForNewStartTilt = False
                else:
                    perPointTiltSeries.append([fileName])
                lastTilt = float(fileName.split("_")[2].replace(self.mdoc_suffix, ""))
            else:
                lookForNewStartTilt = True
                if lastTilt != float(fileName.split("_")[2].replace(self.mdoc_suffix, "")):
                    perPointTiltCounter = 0
                perPointTiltSeries[perPointTiltCounter].append(fileName)
                perPointTiltCounter += 1
                lastTilt = float(fileName.split("_")[2].replace(self.mdoc_suffix, ""))

        self.writePerPointTitlSeriesListFile(perPointTiltSeries)

    def validateInputs(self):

        if not os.path.isabs(self.tilts_dir):
            self.tilts_dir = os.getcwd() +"/"+ self.tilts_dir

        if not os.path.isabs(self.out_dir):
            self.out_dir = os.getcwd() +"/"+ self.out_dir

        if not os.path.exists(self.out_dir):
            if self.versbosity > 1 : print("Creating output directory: %s" % self.out_dir)
            os.makedirs(self.out_dir)

    def main(self):
        self.validateInputs()

        tilt_files = os.listdir(self.tilts_dir)

        # remove all non-mdoc
        tilt_files = [x for x in tilt_files if "mdoc" in x]
        # sort by acquisition serial no
        tilt_files.sort(key=lambda x: x.split("_")[1])

        self.processFilesInList(tilt_files)

        if self.versbosity > 0 : print("All done! Have fun!")


if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Performs merging of separate movie MDOC files from PACE-tomo or dose symmetric SerialEM script into tilt-series MDOC files.")
    add = parser.add_argument
    add('--i', help="Input mdoc directory.")
    add('--o', help="Output directory of merged mdocs.")
    add('--start_tilt', type=int, default=0,
        help="Angle of the \"zero-tilt\" - at which the tilt-series begin. Useful for pretilted tilt-series. (Default: 0)" )
    add('--start_series_nr', type=int, default=1,
        help="The sequential number of the start nr. of tilt series (used for continuing the process). (Default: 1)")
    add('--prefix', type=str, default="tomo",
        help="Prefix of the output tilt series file. (Default: tomo)")
    add('--mdoc_suffix', type=str, default=".tif.mdoc",
        help="Suffix of the *.mdoc files after the mic name. (Default: .tif.mdoc)")
    add('--pre_dose', type=float, default=0.0,
        help="Dose applied before acquiring the first tilt. [e-/A^2]. (Default: 0.0)")
    add('--dose', type=float, default=0.0,
        help="Dose applied per tilt. [e-/A^2]. (Default: 0.0)")
    add('--verb', type=int, default="1",
        help="Verbosity level (0,1,2). (Default: 1)")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()

    mergeMdocs = MergeMdoc(args.i, args.o, args.start_tilt, args.start_series_nr, args.prefix, args.mdoc_suffix, args.pre_dose, args.dose, args.verb)

    mergeMdocs.main()