#!/usr/bin/env python3

import signal
import time
import os
import sys
import subprocess
import argparse
from shutil import which

class MergeTiltsFlow:

    def __init__(self, tilts_dir, out_dir, startTilt=0, startTiltSeriesNr=1, tomo_prefix='tomo', mdoc_dir="", mdoc_suffix=".tif.mdoc", pre_dose=0.0, dose=0.0, apix=1.0, newstackcmd = 'newstack', alterheadercmd = 'alterheader', path='', watch_interval=1, verbosity = 1):
        self.tilts_dir = tilts_dir
        self.out_dir = out_dir
        self.startTilt = startTilt
        self.lastTiltSeriesNr = startTiltSeriesNr
        self.tomo_prefix = tomo_prefix
        self.mdoc_dir = mdoc_dir
        self.mdoc_suffix = mdoc_suffix
        self.pre_dose = pre_dose
        self.dose = dose
        self.apix = apix
        self.newstackcmd = newstackcmd
        self.alterheadercmd = alterheadercmd
        self.watch_interval = watch_interval
        self.versbosity = verbosity
        self.my_env = os.environ.copy()
        self.my_env["PATH"] = path + self.my_env["PATH"]

    class GracefulKiller:
        # from https://stackoverflow.com/questions/18499497
        kill_now = False

        def __init__(self):
            signal.signal(signal.SIGINT, self.exit_gracefully)
            signal.signal(signal.SIGTERM, self.exit_gracefully)

        def exit_gracefully(self, *args):
            self.kill_now = True

    def processMdocFiles(self, mrcTiltList, tomoNr):
        #generate rawtlt file
        rawTltName = self.out_dir + "/" + self.tomo_prefix + str(tomoNr) +"/"+ self.tomo_prefix + str(tomoNr)+".rawtlt"
        with open(rawTltName, 'w') as rawTltFile:
            for mrcTiltName in mrcTiltList:
                mdocName = self.mdoc_dir+"/"+mrcTiltName[:-4]+self.mdoc_suffix
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
            for mrcTiltName in mrcTiltList:
                mdocName = self.mdoc_dir+"/"+mrcTiltName[:-4]+self.mdoc_suffix
                with open(mdocName, 'r') as mdocFile:
                    for mdocLine in mdocFile.readlines():
                        if "FrameSet" in mdocLine:
                            frameSetFound = True
                            firstMdoc = False
                            tomoMdocFile.write("[ZValue = %s ]\n" % tomoCounter)
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
        mrcTiltList.sort(key=lambda x: x.split("_")[1])
        tiltCounter = 0
        tiltDose = self.pre_dose
        tiltDoseList = []
        tltOrderName = self.out_dir + "/" + self.tomo_prefix + str(tomoNr) + "/" + self.tomo_prefix + str(
            tomoNr) + ".tltorder"
        with open(tltOrderName, 'w') as tltOrderFile:
            for mrcTiltName in mrcTiltList:
                tiltCounter += 1
                tiltDose += self.dose
                mdocName = self.mdoc_dir + "/" + mrcTiltName[:-4] + self.mdoc_suffix
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

            perPointTiltSeries[i].sort(key=lambda x: float(x.split("_")[2].replace(".tif", "")))
            # write txt file with sorted tilt files
            with open('%s%s.txt' % (self.tomo_prefix, self.lastTiltSeriesNr + i), 'w') as f:
                f.write("%s\n" % len(perPointTiltSeries[i]))
                for tiltMrcFileName in perPointTiltSeries[i]:
                    f.write(tiltMrcFileName + "\n" + "0\n")
                if self.versbosity > 1 : print('%s%s.txt written' % (self.tomo_prefix, self.lastTiltSeriesNr + i))

            # create merged tilt-series mrc by IMOD newstack
            # create tomo dir
            if not os.path.exists("%s/%s%s" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i)):
                if self.versbosity > 1 : print("Creating %s directory in output path..." % self.tomo_prefix)
                try:
                    os.makedirs("%s/%s%s" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i))
                except:
                    print("ERROR: Could not create output directory: %s/%s%s" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i))
                    return False

            newStackCommand =[self.newstackcmd, "-mode", "2", "-filei", "%s%s.txt" % (self.tomo_prefix, self.lastTiltSeriesNr + i), "-ou", "%s/%s%s/%s%s.mrc" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i, self.tomo_prefix, self.lastTiltSeriesNr + i)  ]
            if self.versbosity > 1 : print("Running: %s" % " ".join(newStackCommand))
            newstackProcess = subprocess.run(newStackCommand, env=self.my_env, stdout=subprocess.DEVNULL)
            if newstackProcess.returncode != 0:
                print("Something bad happened during running the newstack command.")
                print(newstackProcess.stderr)
                return False

            alterHeaderCommand =[self.alterheadercmd, "-d", "%0.3f,%0.3f,%0.3f" % (self.apix, self.apix, self.apix), "%s/%s%s/%s%s.mrc" % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i, self.tomo_prefix, self.lastTiltSeriesNr + i)  ]
            if self.versbosity > 1 : print("Running: %s" % " ".join(alterHeaderCommand))
            alterHeaderProcess = subprocess.run(alterHeaderCommand, env=self.my_env, stdout=subprocess.DEVNULL)
            if alterHeaderProcess.returncode != 0:
                print("Something bad happened during running the alterheader command.")
                print(alterHeaderProcess.stderr)
                return False


            if self.versbosity > 0 : print(" -> %s/%s%s/%s%s.mrc created." % (self.out_dir, self.tomo_prefix, self.lastTiltSeriesNr + i, self.tomo_prefix, self.lastTiltSeriesNr + i))

            # remove original tilts and the text file
            for tiltMrcFileName in perPointTiltSeries[i]:
                try:
                    os.remove(tiltMrcFileName)
                except:
                    print("ERROR: Could not remove tilt file: %s" % tiltMrcFileName)
                    return False

            try:
                os.remove('%s%s.txt' % (self.tomo_prefix, self.lastTiltSeriesNr + i))
            except:
                print("ERROR: Could not remove tilt-series description file: %s%s.txt" % (self.tomo_prefix, self.lastTiltSeriesNr + i))
                return False

            # Process mdoc files
            if self.mdoc_dir != "":
                self.processMdocFiles(perPointTiltSeries[i], self.lastTiltSeriesNr + i)

        self.lastTiltSeriesNr += (i + 1)
        return True

    def processFilesInList(self, tilt_files):
        lastAmountOfFiles = len(tilt_files)
        perPointTiltSeries = []
        lookForNewStartTilt = True

        for fileName in tilt_files:
            if float(fileName.split("_")[2].replace(".tif", "")) == self.startTilt:
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
                lastTilt = float(fileName.split("_")[2].replace(".tif", ""))
            else:
                lookForNewStartTilt = True
                if lastTilt != float(fileName.split("_")[2].replace(".tif", "")):
                    perPointTiltCounter = 0
                perPointTiltSeries[perPointTiltCounter].append(fileName)
                perPointTiltCounter += 1
                lastTilt = float(fileName.split("_")[2].replace(".tif", ""))

        return lastAmountOfFiles, perPointTiltSeries

    def validateInputs(self):
        if which(self.newstackcmd) is None:
            print("IMOD newstack command not found! Make sure it is in path.")
            exit()

        if not self.mdoc_dir == "" and not os.path.isabs(self.mdoc_dir):
            self.mdoc_dir = os.getcwd() +"/"+ self.mdoc_dir

        if not self.mdoc_dir == "" and not os.path.exists(self.mdoc_dir):
            print("Mdoc dir not found! Check if %s exists!" % self.mdoc_dir)
            exit()

        if not os.path.isabs(self.tilts_dir):
            self.tilts_dir = os.getcwd() +"/"+ self.tilts_dir

        if not os.path.isabs(self.out_dir):
            self.out_dir = os.getcwd() +"/"+ self.out_dir

        if not os.path.exists(self.out_dir):
            if self.versbosity > 1 : print("Creating output directory: %s" % self.out_dir)
            os.makedirs(self.out_dir)

    def main(self):
        self.validateInputs()

        killer = self.GracefulKiller()
        lastAmountOfFiles = 0

        # watchdog for new files
        if self.versbosity > 0 : print("Starting the processing watchdog processing... Stop by sending SIGTERM (Ctrl-C).")

        while not killer.kill_now:
            tilt_files = os.listdir(self.tilts_dir)

            if lastAmountOfFiles < len(tilt_files):
                # remove all non-mrc
                tilt_files = [x for x in tilt_files if "mrc" in x]
                # sort by acquisition serial no
                tilt_files.sort(key=lambda x: x.split("_")[1])

                lastAmountOfFiles, perPointTiltSeries = self.processFilesInList(tilt_files)

            if self.watch_interval > 0:
                time.sleep(self.watch_interval)
            else:
                killer.kill_now = True

        # write the last set of tilt-series whe SIGTERM received
        if self.versbosity > 0 : print("\nSIGTERM received. Writing the last tilt-series...")
        self.writePerPointTitlSeriesListFile(perPointTiltSeries)
        if self.versbosity > 0 : print("All done! Have fun!")


if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Performs merging of separate tilt files from PACE-tomo or dose symmetric SerialEM script into tilt-series files.")
    add = parser.add_argument
    add('--i', help="Input tilts directory.")
    add('--o', help="Output directory of merged tilt series.")
    add('--start_tilt', type=int, default=0,
        help="Angle of the \"zero-tilt\" - at which the tilt-series begin. Useful for pretilted tilt-series. (Default: 0)" )
    add('--start_series_nr', type=int, default=1,
        help="The sequential number of the start nr. of tilt series (used for continuing the process). (Default: 1)")
    add('--prefix', type=str, default="tomo",
        help="Prefix of the output tilt series file. (Default: tomo)")
    add('--mdoc_dir', type=str, default="tomo",
        help="Directory where *.mdoc files are located. If empty, no *.rawtlt and merged *.mdoc are generated (Default: empty)")
    add('--mdoc_suffix', type=str, default=".tif.mdoc",
        help="Suffix of the *.mdoc files after the mic name. (Default: .tif.mdoc)")
    add('--pre_dose', type=float, default=0.0,
        help="Dose applied before acquiring the first tilt. [e-/A^2]. (Default: 0.0)")
    add('--dose', type=float, default=0.0,
        help="Dose applied per tilt. [e-/A^2]. (Default: 0.0)")
    add('--apix', type=float, default=1.0,
        help="Pixel size of the merged tilt-series. [A/px]. (Default: 1.0)")
    add('--newstack_cmd', type=str, default="newstack",
        help="Name of the IMOD newstack command. (Default newstack)")
    add('--newstack_cmd', type=str, default="alterheader",
        help="Name of the IMOD alterheader command. (Default alterheader)")
    add('--path', type=str, default="",
        help="Path to be included in env PATH for IMOD! (Default: empty)")
    add('--watch_int', type=int, default="1",
        help="Time interval in seconds for the watchdog. If set to -1 the watchdog is deactivated and runs a single-batch. (Default: 1)")
    add('--verb', type=int, default="1",
        help="Verbosity level (0,1,2). (Default: 1)")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()

    mergeTiltsWatchdog = MergeTiltsFlow(args.i, args.o, args.start_tilt, args.start_series_nr, args.prefix, args.mdoc_dir, args.mdoc_suffix, args.pre_dose, args.dose, args.apix, args.newstack_cmd, args.alterheader_cmd, args.path, args.watch_int, args.verb)

    mergeTiltsWatchdog.main()