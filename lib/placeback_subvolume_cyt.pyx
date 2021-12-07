from lib.metadata import MetaData
import numpy as np
cimport numpy as np
cimport cython
from cython.parallel import prange
from libc.math cimport sin, cos, M_PI, round
import struct
import sys
import time

np.import_array()

DTYPE = np.float32
ctypedef float DTYPE_t

def readMrcSizeApix(mrcFileName):
    cdef int imageSizeX, imageSizeY, imageSizeZ, mx
    cdef float xlen, apix, originX, originY, originZ
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

def readMrcData(mrcFileName):
    cdef int imageSizeX, imageSizeY, imageSizeZ
    cdef np.ndarray[DTYPE_t, ndim=1] mrcData

    with open(mrcFileName, "rb") as mrcFile:
        imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
        mrcData = np.fromfile(mrcFile, dtype=np.dtype(np.float32), count=(imageSizeX * imageSizeY * imageSizeZ),
                              offset=1024 - 12)
    return mrcData

def writeMrcFile(np.ndarray[DTYPE_t, ndim=1] mrcData, stencilFile, outFile):
    with open(stencilFile, "rb") as mrcStencilFile:
        mrcHeader = mrcStencilFile.read(1024)
    with open(outFile, 'wb+') as mrcFile:
        mrcFile.write(mrcHeader)
        mrcFile.seek(12, 0)
        mrcFile.write(b"\x02\x00")
        mrcFile.seek(1024, 0)
        mrcData.astype('float32').tofile(mrcFile)

def writeCmmFile(cmm_coordinates, outputCmmFile, inputMapMaxX, apix, mapOriginX, mapOriginY, mapOriginZ):
    with open(outputCmmFile, 'w') as cmmFile:
        cmmFile.write("<marker_set name=\"marker set 1\">\n")
        id = 1
        for coordinate in cmm_coordinates:
            cmmFile.write(
                "<marker id=\"%d\" x=\"%0.3f\" y=\"%0.3f\" z=\"%0.3f\" r=\"%0.3f\" g=\"%0.3f\" b=\"%0.3f\" radius=\"%f\" coordX=\"%0.3fpx\" coordY=\"%0.3fpx\" coordZ=\"%0.3fpx\"/> \n" % (
                    id, -coordinate[0] + coordinate[3] * apix + mapOriginX,
                    -coordinate[1] + coordinate[4] * apix + mapOriginY,
                    -coordinate[2] + coordinate[5] * apix + mapOriginZ, coordinate[6][0],
                    coordinate[6][1], coordinate[6][2], inputMapMaxX / 6, coordinate[3], coordinate[4], coordinate[5]))
            id += 1
        cmmFile.write("</marker_set>")

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef matrix_from_euler(DTYPE_t rot, DTYPE_t tilt, DTYPE_t psi):
    """create a rotation matrix from three Euler anges in ZYZ convention"""
    cdef np.ndarray[DTYPE_t, ndim=2] m = np.full((3, 3), 0.0, dtype=DTYPE)

    m[0, 0] = cos(psi) * cos(tilt) * cos(rot) - sin(psi) * sin(rot)
    m[0, 1] = cos(psi) * cos(tilt) * sin(rot) + sin(psi) * cos(rot)
    m[0, 2] = -cos(psi) * sin(tilt)
    m[1, 0] = -sin(psi) * cos(tilt) * cos(rot) - cos(psi) * sin(rot)
    m[1, 1] = -sin(psi) * cos(tilt) * sin(rot) + cos(psi) * cos(rot)
    m[1, 2] = sin(psi) * sin(tilt)
    m[2, 0] = sin(tilt) * cos(rot)
    m[2, 1] = sin(tilt) * sin(rot)
    m[2, 2] = cos(tilt)
    return m

cdef DTYPE_t c_radians(DTYPE_t degree):
    return degree * (M_PI / 180.0);

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def rotateVolume(np.ndarray[DTYPE_t, ndim=1] mrcData, int sizeX, int sizeY, int sizeZ, DTYPE_t rot, DTYPE_t tilt,
                 DTYPE_t psi):
    cdef int MapMaxX = sizeX, MapMaxY = sizeY, MapMaxZ = sizeZ
    cdef np.ndarray[DTYPE_t, ndim=1] rotatedMapArray = np.full((MapMaxX * MapMaxY * MapMaxZ), 0.0, dtype=DTYPE)

    cdef int pixPos = 0
    cdef int rotatedX, rotatedY, rotatedZ
    cdef DTYPE_t origXi, origYi, origZi, pixelValue

    # use inverse rotation matrix as we are looking for the origin of the rotated pixe vecto@cython.boundscheck(False)
    cdef np.ndarray[DTYPE_t, ndim=2] rotMatrixInv = np.linalg.inv(
        matrix_from_euler(c_radians(rot), c_radians(tilt), c_radians(psi)))

    cdef int rotatedMapArraySize = len(rotatedMapArray)
    cdef DTYPE_t d000, d001, d010, d011, d100, d101, d110, d111, dx00, dx01, dx10, dx11, dxy0, dxy1
    cdef int x0, y0, z0, x1, y1, z1
    cdef DTYPE_t fx, fy, fz

    # for every pixel position in the rotated map look for the unrotate original pixel position and trilinear interpolate the value of the original pixel
    for pixPos in prange(rotatedMapArraySize, nogil=True):
        #for pixPos in range(rotatedMapArraySize):
        rotatedZ = pixPos // (MapMaxX * MapMaxY)
        rotatedY = (pixPos % (MapMaxX * MapMaxY)) // MapMaxY
        rotatedX = (pixPos % (MapMaxX * MapMaxY)) % MapMaxY

        origXi = (rotMatrixInv[0, 0] * (rotatedX - MapMaxX / 2) + rotMatrixInv[0, 1] * (rotatedY - MapMaxY / 2) +
                  rotMatrixInv[0, 2] * (rotatedZ - MapMaxZ / 2))
        origYi = (rotMatrixInv[1, 0] * (rotatedX - MapMaxX / 2) + rotMatrixInv[1, 1] * (rotatedY - MapMaxY / 2) +
                  rotMatrixInv[1, 2] * (rotatedZ - MapMaxZ / 2))
        origZi = (rotMatrixInv[2, 0] * (rotatedX - MapMaxX / 2) + rotMatrixInv[2, 1] * (rotatedY - MapMaxY / 2) +
                  rotMatrixInv[2, 2] * (rotatedZ - MapMaxZ / 2))

        origXi = origXi + MapMaxX / 2
        origYi = origYi + MapMaxY / 2
        origZi = origZi + MapMaxZ / 2

        if origXi >= MapMaxX - 1 or origYi >= MapMaxY - 1 or origZi >= MapMaxZ - 1 or origXi < 0 or origYi < 0 or origZi < 0:
            pixelValue = 0
        else:
            # make trilinear interpolation of the pixel value from the original map, using interpolation between origXi, origYi, origZi and their integer parts
            if <int> origXi >= MapMaxX - 1:
                x0 = MapMaxX - 2
            else:
                x0 = <int> origXi
            fx = origXi - x0
            x1 = x0 + 1

            if <int> origYi >= MapMaxY - 1:
                y0 = MapMaxY - 2
            else:
                y0 = <int> origYi
            fy = origYi - y0
            y1 = y0 + 1

            if <int> origZi >= MapMaxZ - 1:
                z0 = MapMaxZ - 2
            else:
                z0 = <int> origZi
            fz = origZi - z0
            z1 = z0 + 1

            d000 = mrcData[x0 + y0 * MapMaxX + z0 * MapMaxX * MapMaxY]
            d001 = mrcData[x1 + y0 * MapMaxX + z0 * MapMaxX * MapMaxY]
            d010 = mrcData[x0 + y1 * MapMaxX + z0 * MapMaxX * MapMaxY]
            d011 = mrcData[x1 + y1 * MapMaxX + z0 * MapMaxX * MapMaxY]
            d100 = mrcData[x0 + y0 * MapMaxX + z1 * MapMaxX * MapMaxY]
            d101 = mrcData[x1 + y0 * MapMaxX + z1 * MapMaxX * MapMaxY]
            d110 = mrcData[x0 + y1 * MapMaxX + z1 * MapMaxX * MapMaxY]
            d111 = mrcData[x1 + y1 * MapMaxX + z1 * MapMaxX * MapMaxY]

            dx00 = ((d000) + ((d001) - (d000)) * (fx))
            dx01 = ((d100) + ((d101) - (d100)) * (fx))
            dx10 = ((d010) + ((d011) - (d010)) * (fx))
            dx11 = ((d110) + ((d111) - (d110)) * (fx))
            dxy0 = ((dx00) + ((dx10) - (dx00)) * (fy))
            dxy1 = ((dx01) + ((dx11) - (dx01)) * (fy))

            pixelValue = ((dxy0) + ((dxy1) - (dxy0)) * (fz))

        rotatedMapArray[rotatedX + rotatedY * MapMaxX + rotatedZ * MapMaxX * MapMaxY] = pixelValue

    return rotatedMapArray

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def placeSubvolumes(inputStarFile, inputVolumeToPlace, outputMapStencil, outputPrefix, outputCmm, binning,
                    placePartialVolumes, recenter, coloringLabel, outputColorMap, colorMapTreshold, colorMapExtend):
    #read in star file
    md = MetaData(inputStarFile)
    particles = []
    for particle in md:
        particles.append(particle)

    optic_groups = []
    for optic_group in md.data_optics:
        optic_groups.append(optic_group)

    #get unbinned apix from star file
    cdef float apix
    apix = float(optic_groups[0].rlnTomoTiltSeriesPixelSize)

    #get MRC size
    cdef int inputMapMaxX, inputMapMaxY, inputMapMaxZ
    cdef float inputApix, inputMapOriginX, inputMapOriginY, inputMapOriginZ
    inputMapMaxX, inputMapMaxY, inputMapMaxZ, inputApix, inputMapOriginX, inputMapOriginY, inputMapOriginZ = readMrcSizeApix(
        inputVolumeToPlace)
    if inputMapMaxX == 0 or inputMapMaxY == 0 or inputMapMaxZ == 0:
        print("Error: One of the dimensions of the input volume to place is 0!")
        exit()

    cdef int mapMaxX, mapMaxY, mapMaxZ
    cdef float mapApix, mapOriginX, mapOriginY, mapOriginZ
    mapMaxX, mapMaxY, mapMaxZ, mapApix, mapOriginX, mapOriginY, mapOriginZ = readMrcSizeApix(outputMapStencil)

    cdef float minColor, maxColor
    if outputCmm:
        # create a 512 steps red/green rainbow array for cmm - each element represents r/g/b (0 = min; 1=max)
        if coloringLabel != "":
            rainbowArray = []
            rangeColor = 512
            for i in range(int(rangeColor / 2)):  # from red -> yellow
                rainbowArray.append([1, i / (rangeColor / 2), 0])
            for i in range(int(rangeColor / 2), 0, -1):  # from yellow -> green
                rainbowArray.append([i / (rangeColor / 2), 1, 0])

            #get the min and max value (ie range) of the desired label in the star file
            minColor = getattr(particles[0], coloringLabel)
            maxColor = getattr(particles[0], coloringLabel)
            for particle in particles:
                minColor = min(minColor, getattr(particle, coloringLabel))
                maxColor = max(maxColor, getattr(particle, coloringLabel))
        else:
            rainbowArray = [[1, 1, 0]]

    print("Output map size [x, y, z]: %i x %i x %i" % (mapMaxX, mapMaxY, mapMaxZ))

    #get binnig if not defined by user
    cdef float coordinateBinningFactor
    if binning == 0.0:
        coordinateBinningFactor = mapApix / apix
        print(
                "Binnig factor not provided. Calculated from stencil header and star file => Binning = %0.2f" % coordinateBinningFactor)
    else:
        coordinateBinningFactor = binning

    if abs(inputApix - mapApix) > 0.01:
        print(
                "!!!WARNING: Stencil tomogram has different apix than subvolume to be placed (%0.3f vs %0.3f). The subvolumes placed in final volume might not be in scale." % (
            mapApix, inputApix))

    if recenter and not (
            hasattr(particles[0], 'rlnOriginXAngst') and hasattr(particles[0], 'rlnOriginYAngst') and hasattr(
            particles[0], 'rlnOriginZAngst')):
        print(
            "WARNING: Recentering was enable but no rlnOriginXAngst/rlnOriginYAngst/rlnOriginZAngst in star file. Disabling recentering!")
        recenter = False

    #initialize outputMapArray array
    cdef np.ndarray[DTYPE_t, ndim=1] outputMapArray = np.full((mapMaxX * mapMaxY * mapMaxZ), 0.0, dtype=DTYPE)

    cdef bint makeColorMap = False
    cdef bint makeCmm = False

    #convert python boolean to bint
    if outputColorMap:
        makeColorMap = True
    if outputCmm:
        makeCmm = True

    cdef np.ndarray[DTYPE_t, ndim=1] outputColorMapArray = np.full((mapMaxX * mapMaxY * mapMaxZ), 0.0, dtype=DTYPE)

    cdef int nrOfParticles = len(particles)
    currentParticleNr = 1

    total_start = time.perf_counter()
    mapToRotate = readMrcData(inputVolumeToPlace)

    cdef float colorMapTreshold_c = colorMapTreshold
    cdef int colorMapExtend_c = int(colorMapExtend)
    cdef float coloringValue
    cdef float xpos, ypos, zpos
    cdef int pixPos, origX, origY, origZ, newX, newY, newZ, x_extend, y_extend, z_extend
    cdef DTYPE_t pixel

    cdef np.ndarray[DTYPE_t, ndim=1] rotatedMap = np.full((inputMapMaxX * inputMapMaxY * inputMapMaxZ), 0.0,
                                                          dtype=DTYPE)
    cdef int rotatedMapArraySize = len(mapToRotate)

    cmm_coordinates = []

    print("Placing %i subvolumes:" % nrOfParticles)
    for particle in particles:
        if recenter:
            xpos = (particle.rlnCoordinateX - particle.rlnOriginXAngst / apix) / coordinateBinningFactor
            ypos = (particle.rlnCoordinateY - particle.rlnOriginYAngst / apix) / coordinateBinningFactor
            zpos = (particle.rlnCoordinateZ - particle.rlnOriginZAngst / apix) / coordinateBinningFactor
        else:
            xpos = particle.rlnCoordinateX / coordinateBinningFactor
            ypos = particle.rlnCoordinateY / coordinateBinningFactor
            zpos = particle.rlnCoordinateZ / coordinateBinningFactor

        if coloringLabel != "":
            coloringValue = round(float(getattr(particle, coloringLabel)) * 100) / 100

        if makeCmm:
            #color value index in raibowArray
            if coloringLabel != "":
                colorIndex = int((coloringValue - minColor) / (maxColor - minColor) * (rangeColor - 1))
            else:
                colorIndex = 0
            if recenter:
                cmm_coordinates.append([particle.rlnOriginXAngst, particle.rlnOriginYAngst, particle.rlnOriginZAngst,
                                        particle.rlnCoordinateX, particle.rlnCoordinateY, particle.rlnCoordinateZ,
                                        rainbowArray[colorIndex]])
            else:
                cmm_coordinates.append(
                    [0, 0, 0, particle.rlnCoordinateX, particle.rlnCoordinateY, particle.rlnCoordinateZ,
                     rainbowArray[colorIndex]])

        rotatedMap = rotateVolume(mapToRotate, inputMapMaxX, inputMapMaxY, inputMapMaxZ, particle.rlnAngleRot,
                                  particle.rlnAngleTilt, particle.rlnAnglePsi)

        if ((round(xpos) + inputMapMaxX / 2) < mapMaxX and (round(ypos) + inputMapMaxY / 2) < mapMaxY and (
                round(zpos) + inputMapMaxZ / 2) < mapMaxZ and (round(xpos) - inputMapMaxX / 2) >= 0 and (
                    round(ypos) - inputMapMaxY / 2) >= 0 and (
                    round(zpos) - inputMapMaxZ / 2) >= 0) or placePartialVolumes:
            # a simple progress bar
            sys.stdout.write('\r')
            progress = int(currentParticleNr / nrOfParticles * 20)
            sys.stdout.write("[%-20s] %d%%" % ('=' * progress, 5 * progress))
            sys.stdout.flush()

            pixPos = 0
            for pixPos in prange(rotatedMapArraySize, nogil=True):
                #for pixPos in range(rotatedMapArraySize):
                pixel = rotatedMap[pixPos]
                origZ = pixPos // (inputMapMaxX * inputMapMaxY)
                origY = (pixPos % (inputMapMaxX * inputMapMaxY)) // inputMapMaxY
                origX = (pixPos % (inputMapMaxX * inputMapMaxY)) % inputMapMaxY
                newZ = <int> round(origZ - inputMapMaxZ / 2 + zpos)
                newY = <int> round(origY - inputMapMaxY / 2 + ypos)
                newX = <int> round(origX - inputMapMaxX / 2 + xpos)
                # data[x][y][z] = data[x + ymaxX + zmaxX*maxY] 
                if newZ < mapMaxZ and newY < mapMaxY and newX < mapMaxX and newZ >= 0 and newY >= 0 and newX >= 0:
                    if outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] == 0:
                        outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] = pixel
                    elif outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] < pixel:
                        outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] = pixel

                    if makeColorMap:
                        if pixel >= colorMapTreshold_c:
                            for x_extend in range(newX-colorMapExtend_c, newX+colorMapExtend_c, 1):
                                for y_extend in range(newY - colorMapExtend_c, newY + colorMapExtend_c, 1):
                                    for z_extend in range(newZ - colorMapExtend_c, newZ + colorMapExtend_c, 1):
                                        if z_extend < mapMaxZ and y_extend < mapMaxY and x_extend < mapMaxX and z_extend >= 0 and y_extend >= 0 and x_extend >= 0:
                                            outputColorMapArray[x_extend + y_extend * mapMaxX + z_extend * mapMaxX * mapMaxY] = coloringValue

        currentParticleNr += 1

    #end the progress bar
    sys.stdout.write('\r\n')

    print("Writing MRC file...")
    writeMrcFile(outputMapArray, outputMapStencil, outputPrefix+".mrc")
    print("%s.mrc written." % outputPrefix)

    if makeColorMap:
        writeMrcFile(outputColorMapArray, outputMapStencil, outputPrefix+"_color.mrc")
        print("%s_color.mrc written." % outputPrefix)

    if outputCmm:
        writeCmmFile(cmm_coordinates, outputPrefix + ".cmm", inputMapMaxX * inputApix, apix, mapOriginX, mapOriginY, mapOriginZ)
        print("%s.cmm written." % outputPrefix)
    print("===>Total run time: %0.2f sec" % ((time.perf_counter() - total_start)))
    print("All done. Have fun!")
