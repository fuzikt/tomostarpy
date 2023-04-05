#!/usr/bin/env python3
import struct
import os
import numpy as np
import subprocess
import time
import shutil
import math
import argparse
import sys
from lib.metadata import *

class CtfFitTiltSeriesCtFind4:
    def __init__(self, tiltsFile, tltFile, outDir, defocus_window = 10000, apix = 1, voltage = 300, cs_val = 2.7, amp_cont = 0.07, spectr_size = 512, min_res = 30, max_res = 8, min_defoc = 5000, max_defoc = 75000, step_defoc = 100, threads = 10, tmpDir = "ctffind_tmp", ctffind_cmd = "ctffind", gctf = False, gpu = 0, gctf_cmd = "gctf", path='', verbosity = 0):
        self.tiltsFile = tiltsFile
        self.outDir = outDir
        self.tltFile = tltFile
        self.defocus_window = defocus_window
        self.apix = apix
        self.voltage = voltage
        self.cs_val = cs_val
        self.amp_cont = amp_cont
        self.spectr_size = spectr_size
        self.min_res = min_res
        self.max_res = max_res
        self.min_defoc = min_defoc
        self.max_defoc = max_defoc
        self.step_defoc = step_defoc
        self.threads = threads
        self.tmpDir = tmpDir
        self.ctffindcmd = ctffind_cmd
        self.gctfcmd = gctf_cmd
        self.doGctf = gctf
        self.gpuID = gpu
        self.versbosity = verbosity

        self.tiltcounter = 1
        self.ctfFind4Results = []
        self.outFilePreffix = os.path.basename(self.tiltsFile)[:-4]
        self.my_env = os.environ.copy()
        self.my_env["PATH"] = path + self.my_env["PATH"]

    def readMrcSizeApix(self, mrcFileName):
        with open(mrcFileName, "rb") as mrcFile:
            imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
            mrcFile.seek(28)
            mx = int(struct.unpack('i', mrcFile.read(4))[0])
            mrcFile.seek(40)
            xlen = float(struct.unpack('f', mrcFile.read(4))[0])
            apix = xlen / mx
            mrcFile.seek(196)
            originX = float(struct.unpack('f', mrcFile.read(4))[0])
            originY = float(struct.unpack('f', mrcFile.read(4))[0])
            originZ = float(struct.unpack('f', mrcFile.read(4))[0])
        return imageSizeX, imageSizeY, imageSizeZ, apix, originX, originY, originZ


    def readMrcData(self, mrcFileName):
        with open(mrcFileName, "rb") as mrcFile:
            imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
            imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
            mrcData = np.fromfile(mrcFile, dtype=np.dtype(np.float32), count=(imageSizeX * imageSizeY * imageSizeZ),
                                  offset=1024 - 12)
        return mrcData


    def writeMrcFile(self, mrcData, stencilFile, zDim, outFile):
        with open(stencilFile, "rb") as mrcStencilFile:
            mrcHeader = mrcStencilFile.read(1024)
        with open(outFile, 'wb+') as mrcFile:
            mrcFile.write(mrcHeader)
            # change Z dimension to 1
            mrcFile.seek(8, 0)
            mrcFile.write(np.array([zDim], dtype=np.int32))
            mrcFile.seek(1024, 0)
            mrcData.astype('float32').tofile(mrcFile)


    def runCtfFind4(self, inputFile, outputPreffix, apix, voltage, cs_val, amp_cont, spectr_size, min_res, max_res, min_defoc,
                    max_defoc, step_defoc, threads):
        if self.versbosity > 0:
            STDOUT = None
        else:
            STDOUT = subprocess.DEVNULL

        p = subprocess.Popen(self.ctffindcmd+" <<EOF\n" + inputFile + "\n"
                             + outputPreffix + "\n"
                             + str(apix) + "\n"
                             + str(voltage) + "\n"
                             + str(cs_val) + "\n"
                             + str(amp_cont) + "\n"
                             + str(spectr_size) + "\n"
                             + str(min_res) + "\n"
                             + str(max_res) + "\n"
                             + str(min_defoc) + "\n"
                             + str(max_defoc) + "\n"
                             + str(step_defoc) + "\n"
                             + "no\n"  # Do you know what astigmatism is present?
                             + "no\n"  # Slower, more exhaustive search?
                             + "no\n"  # Use a restraint on astigmatism?
                             + "no\n"  # Find additional phase shift?
                             + "yes\n"  # Do you want to set expert options?Do you want to set expert options?
                             + "yes\n"  # Resample micrograph if pixel size too small?
                             + "no\n"  # Do you already know the defocus?
                             + str(threads) + "\n"
                             + "\nEOF", shell=True, stdout=STDOUT, stderr=subprocess.STDOUT)
        p.wait()

    def runGctf(self, inputFile, outputStar, apix, voltage, cs_val, amp_cont, spectr_size, min_res, max_res, min_defoc,
                max_defoc, step_defoc):
        if self.versbosity > 0:
            STDOUT = None
        else:
            STDOUT = subprocess.DEVNULL

        p = subprocess.Popen(self.gctfcmd
                             + " --apix " + str(apix)
                             + " --kv " + str(voltage)
                             + " --cs " + str(cs_val)
                             + " --ac " + str(amp_cont)
                             + " --do_EPA 1"
                             + " --boxsize " + str(spectr_size)
                             + " --resL " + str(min_res)
                             + " --resH " + str(max_res)
                             + " --defL " + str(min_defoc)
                             + " --defH " + str(max_defoc)
                             + " --defS " + str(step_defoc)
                             + " --gid " + str(self.gpuID)
                             + " --ctfstar " + outputStar
                             + " " + inputFile, shell=True, stdout=STDOUT, stderr=subprocess.STDOUT)
        p.wait()


    def runCtfFindOnSeriesOfTilts(self, previousTiltAngleDefocus, tiltValues, rangeOfTilts):
        for tiltIndex in rangeOfTilts:
            min_defoc_restr = max(0, previousTiltAngleDefocus - self.defocus_window)
            max_defoc_restr = previousTiltAngleDefocus + self.defocus_window
            max_res_weighted = self.max_res * 1 / math.cos(math.radians(tiltValues[tiltIndex]))
            print("processing tilt %s/%s" % (self.tiltcounter, len(tiltValues)))
            self.tiltcounter += 1
            if not self.doGctf:
                self.runCtfFind4("%s/tilt_%s.mrc" % (self.tmpDir, tiltIndex),
                                 "%s/tilt_%s_ctf.mrc" % (self.tmpDir, tiltIndex), self.apix, self.voltage, self.cs_val,
                                 self.amp_cont, self.spectr_size,
                                 self.min_res, max_res_weighted, min_defoc_restr, max_defoc_restr, self.step_defoc,
                                 self.threads)
                ctfFind4Result = self.getCtfFind4Results("%s/tilt_%s_ctf.txt" % (self.tmpDir, tiltIndex))
            else:
                self.runGctf("%s/tilt_%s.mrc" % (self.tmpDir, tiltIndex),
                                "%s/tilt_%s_ctf.star" % (self.tmpDir, tiltIndex), self.apix, self.voltage, self.cs_val,
                                 self.amp_cont, self.spectr_size,
                                 self.min_res, max_res_weighted, min_defoc_restr, max_defoc_restr, self.step_defoc)
                ctfFind4Result = self.getGctfResults(
                    "%s/tilt_%s_ctf.star" % (self.tmpDir, tiltIndex))
                shutil.move("%s/tilt_%s.ctf" % (self.tmpDir, tiltIndex),
                            "%s/tilt_%s_ctf.mrc" % (self.tmpDir, tiltIndex))

            # set correct tilt nr
            ctfFind4Result[0] = tiltIndex + 1
            self.ctfFind4Results.append(ctfFind4Result)

            previousTiltAngleDefocus = (ctfFind4Result[1] + ctfFind4Result[2]) / 2

            if previousTiltAngleDefocus < 0:
                previousTiltAngleDefocus = 0.0

    def findSmallestTiltFromTlt(self, tiltAngles):
        minTiltAngle = 1000
        for tiltAngle in tiltAngles:
            if abs(tiltAngle) < abs(minTiltAngle):
                minTiltAngle = tiltAngle
        return tiltAngles.index(minTiltAngle)


    def readTiltFile(self, tltFile):
        tiltValues = []
        with open(tltFile, 'r') as tiltValuesFile:
            for tiltValue in tiltValuesFile.readlines():
                tiltValues.append(float(tiltValue))
        return tiltValues


    def getCtfFind4Results(self, diagnosticTxtFile):
        with open(diagnosticTxtFile, 'r') as diagnosticFile:
            for line in diagnosticFile.readlines():
                if line.split()[0] != "#":
                    return [float(i) for i in line.split()]

    def getGctfResults(self, diagnosticStarFile):
        mdDefocusData = MetaData(diagnosticStarFile)
        return [1.0, mdDefocusData.data_[0].rlnDefocusU, mdDefocusData.data_[0].rlnDefocusV, mdDefocusData.data_[0].rlnDefocusAngle, 0.0, mdDefocusData.data_[0].rlnCtfFigureOfMerit, mdDefocusData.data_[0].rlnFinalResolution]


    def writeCtfFind4ResultsSumary(self, ctfFind4Results, outputFile):
        with open(outputFile, 'w') as ctfFind4outFile:
            for ctfFind4Result in ctfFind4Results:
                ctfFind4outFile.write("%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n" % tuple(ctfFind4Result))


    def combineCtfMrcOutput(self, nrOfTilts, tmpDir, outFile):
        combinedData = self.readMrcData("%s/tilt_%s_ctf.mrc" % (tmpDir, 0))

        for tilt in range(nrOfTilts - 1):
            combinedData = np.append(combinedData, self.readMrcData("%s/tilt_%s_ctf.mrc" % (tmpDir, tilt + 1)), 0)

        self.writeMrcFile(combinedData, "%s/tilt_%s_ctf.mrc" % (tmpDir, 0), nrOfTilts, outFile)

    def main(self):
        total_start = time.perf_counter()

        if not os.path.exists(self.tmpDir):
            os.makedirs(self.tmpDir)

        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)

        nonSplitMrcData = self.readMrcData(self.tiltsFile);
        imageSizeX, imageSizeY, imageSizeZ, pixelSize, originX, originY, originZ = self.readMrcSizeApix(self.tiltsFile)

        imageSizeXY = imageSizeX * imageSizeY

        perTiltMrcData = nonSplitMrcData.reshape(-1, imageSizeXY)

        print("Preparing data %s for CtfFind4 fit...." % self.tiltsFile )
        for i in range(imageSizeZ):
            self.writeMrcFile(perTiltMrcData[i], self.tiltsFile, 1, "%s/tilt_%s.mrc" % (self.tmpDir, i))

        tiltValues = self.readTiltFile(self.tltFile)
        smallestTiltAngleIndex = self.findSmallestTiltFromTlt(tiltValues)

        # ctf estimate center tilt
        print("processing tilt %s/%s" % (self.tiltcounter, len(tiltValues)))
        self.tiltcounter += 1

        if not self.doGctf:
            self.runCtfFind4("%s/tilt_%s.mrc" % (self.tmpDir, smallestTiltAngleIndex),
                        "%s/tilt_%s_ctf.mrc" % (self.tmpDir, smallestTiltAngleIndex), self.apix, self.voltage, self.cs_val, self.amp_cont, self.spectr_size,
                        self.min_res, self.max_res, self.min_defoc, self.max_defoc, self.step_defoc, self.threads)
            ctfFind4ResultSmalestAngle = self.getCtfFind4Results("%s/tilt_%s_ctf.txt" % (self.tmpDir, smallestTiltAngleIndex))

        else:
            self.runGctf("%s/tilt_%s.mrc" % (self.tmpDir, smallestTiltAngleIndex),
                        "%s/tilt_%s_ctf.star" % (self.tmpDir, smallestTiltAngleIndex), self.apix, self.voltage, self.cs_val, self.amp_cont, self.spectr_size,
                        self.min_res, self.max_res, self.min_defoc, self.max_defoc, self.step_defoc)
            ctfFind4ResultSmalestAngle = self.getGctfResults("%s/tilt_%s_ctf.star" % (self.tmpDir, smallestTiltAngleIndex))
            shutil.move("%s/tilt_%s.ctf" % (self.tmpDir, smallestTiltAngleIndex), "%s/tilt_%s_ctf.mrc" % (self.tmpDir, smallestTiltAngleIndex))

#        exit()
        # set correct tilt nr
        ctfFind4ResultSmalestAngle[0] = smallestTiltAngleIndex + 1

        previousTiltAngleDefocus = (ctfFind4ResultSmalestAngle[1] + ctfFind4ResultSmalestAngle[2]) / 2

        self.runCtfFindOnSeriesOfTilts(previousTiltAngleDefocus, tiltValues, range(smallestTiltAngleIndex - 1, -1, -1))

        self.ctfFind4Results.reverse()
        self.ctfFind4Results.append(ctfFind4ResultSmalestAngle)
        previousTiltAngleDefocus = (ctfFind4ResultSmalestAngle[1] + ctfFind4ResultSmalestAngle[2]) / 2

        self.runCtfFindOnSeriesOfTilts(previousTiltAngleDefocus, tiltValues, range(smallestTiltAngleIndex + 1, len(tiltValues)))

        self.writeCtfFind4ResultsSumary(self.ctfFind4Results, self.outDir + "/" + self.outFilePreffix + "_ctf.txt")
        self.combineCtfMrcOutput(len(tiltValues), self.tmpDir, self.outDir + "/" + self.outFilePreffix + "_ctf.mrc")

        # cleanup
        try:
            shutil.rmtree(self.tmpDir)
        except OSError as e:
            print("Error: %s - %s." % (e.filename, e.strerror))

        print("%s written...." % (self.outFilePreffix + "_ctf.txt") )
        print("%s written...." % (self.outFilePreffix + "_ctf.mrc"))
        print("===>Total run time: %0.2f sec" % ((time.perf_counter() - total_start)))
        print("All done. Have fun!")

if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Performs merging of separate tilt files from PACE-tomo or dose symmetric SerialEM script into tilt-series files.")
    add = parser.add_argument
    add('--i', help="Input tilt-series MRC file.")
    add('--o', help="Output directory of results.")
    add('--tlt', help="*.tlt or *.rawtlt file of the tilt series" )
    add('--defocus_window', type=int, default=10000,
        help="Defocus window in Angstroms used for the restrained ctf search of the subsequent tilts. (Default: 10 000)")
    add('--apix', type=float, default=1.0,
        help="Pixel-size in Angstroms. (Default: 1.0)")
    add('--voltage', type=float, default=300.0,
        help="Acceleration voltage of the microscope in kV (Default: 300.0)")
    add('--cs_val', type=float, default=2.7,
        help="Cs value of the microscope. (Default: 2.7)")
    add('--amp_cont', type=float, default=0.07,
        help="Amplitude contrast used for CTF fitting. (Default: 0.07)")
    add('--spectr_size', type=int, default=512,
        help="Size of the Fourier spectrum used for CTF fitting. (Default: 512)")
    add('--min_res', type=float, default=30.0,
        help="Minimum resolution in Angstroms used for zero-tilt CTF fitting. (Default: 30.0)")
    add('--max_res', type=float, default=8.0,
        help="Maximum resolution  in Angstroms used for zero-tilt CTF fitting. (Default: 8.0)")
    add('--min_defoc', type=float, default=5000.0,
        help="Minimum defocus in Angstroms used for zero-tilt CTF fitting. (Default: 5000.0)")
    add('--max_defoc', type=float, default=75000.0,
        help="Maximum defocus in Angstroms used for zero-tilt CTF fitting. (Default: 75000.0)")
    add('--step_defoc', type=float, default=100.0,
        help="Defocus step in Angstroms used for CTF fitting. (Default: 100.0)")
    add('--threads', type=int, default=10,
        help="Number of parallel threads used for calculation. (Default: 10)")
    add('--gctf', dest='gctf', action='store_true', default=False,
        help="Use gCTF instead of CtfFind4.")
    add('--gpu', type=int, default=0,
        help="GPU card ID to be used by gCTF. (Default: 0)")
    add('--tmpDir', type=str, default="ctffind_tmp",
        help="Temp directory for storage of the intermediate files. (Default: ctffind_tmp)")
    add('--ctffind_cmd', type=str, default="ctffind",
        help="Name of the Ctffind4 command. (Default: ctffind)")
    add('--gctf_cmd', type=str, default="gctf",
        help="Name of the gCTF command. (Default: gctf)")
    add('--path', type=str, default="",
        help="Path to be included in env PATH for Ctffind4. (Default: empty)")

    add('--verb', type=int, default="0",
        help="Verbosity level (0,1). (Default: 0)")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()

    ctfFitTilts = CtfFitTiltSeriesCtFind4(args.i, args.tlt, args.o, defocus_window=args.defocus_window, apix=args.apix, voltage=args.voltage, cs_val=args.cs_val, amp_cont=args.amp_cont, spectr_size=args.spectr_size, min_res=args.min_res, max_res=args.max_res, min_defoc=args.min_defoc, max_defoc=args.max_defoc, step_defoc=args.step_defoc, threads=args.threads, tmpDir=args.tmpDir, ctffind_cmd=args.ctffind_cmd, gctf=args.gctf, gpu=args.gpu, gctf_cmd=args.gctf_cmd, path=args.path, verbosity=args.verb)

    ctfFitTilts.main()