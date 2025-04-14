from lib.metadata import MetaData
import numpy as np
from lib.euler import euler_from_matrix
import lib.matrix3 as matrix3
from lib.chimcmm import cmmData
from math import sin, cos, ceil, radians, pi
import struct
import sys
import time
import ctypes
from cuda import cuda
from lib.cuda_common import common
from lib.cuda_common.helper_cuda import checkCudaErrors, findCudaDeviceDRV

placeBackKernels = """
        extern "C"{
              __global__ void rotate(float *mrcData_gpu, float *rotatedMapArray_gpu, float *rotMatrixInv_g, float shiftX, float shiftY, float shiftZ, int MapMaxX, int MapMaxY, int MapMaxZ )
              {
                size_t pixPos = blockIdx.x * blockDim.x + threadIdx.x;
                __shared__ float rotMatrixInv_gpu[9];
                if (threadIdx.x < 9){
                 rotMatrixInv_gpu[threadIdx.x] = rotMatrixInv_g[threadIdx.x];
                }               
                __syncthreads();
                int MapMaxXY = MapMaxX * MapMaxY;


                float pixelValue = 0.0;      
                if(pixPos < (MapMaxX*MapMaxY*MapMaxZ)){

                int rotatedZ = floorf(pixPos / (MapMaxXY));
                int rotatedY = floorf((pixPos % (MapMaxXY)) / MapMaxY);
                int rotatedX = (pixPos % (MapMaxXY)) % MapMaxY;

                float origXi = ((rotMatrixInv_gpu[0] * (rotatedX - MapMaxX / 2) + rotMatrixInv_gpu[1] * (rotatedY - MapMaxY / 2) + rotMatrixInv_gpu[2] * (rotatedZ - MapMaxZ / 2))) - shiftX;
                float origYi = ((rotMatrixInv_gpu[3] * (rotatedX - MapMaxX / 2) + rotMatrixInv_gpu[4] * (rotatedY - MapMaxY / 2) + rotMatrixInv_gpu[5] * (rotatedZ - MapMaxZ / 2))) - shiftY;
                float origZi = ((rotMatrixInv_gpu[6] * (rotatedX - MapMaxX / 2) + rotMatrixInv_gpu[7] * (rotatedY - MapMaxY / 2) + rotMatrixInv_gpu[8] * (rotatedZ - MapMaxZ / 2))) - shiftZ;

                origXi = origXi + MapMaxX / 2;
                origYi = origYi + MapMaxY / 2;
                origZi = origZi + MapMaxZ / 2;

                if ((origXi > MapMaxX) || (origYi > MapMaxY) || (origZi > MapMaxZ) || (origXi < 0.0) || (origYi < 0.0) || (origZi < 0.0)){
                    pixelValue = 0.0; }
                else{
                    // make trilinear interpolation of the pixel value from the original map, using interpolation between origXi, origYi, origZi and their integer parts

                    int x0 = 0;
                    if (origXi > MapMaxX){
                        x0 = MapMaxX - 1; }
                    else{
                        x0 = floorf(origXi); }
                    float fx = origXi - x0;
                    int x1 = x0 + 1;

                    int y0 = 0;
                    if (origYi > MapMaxY){
                        y0 = MapMaxY - 1; }
                    else{
                        y0 = floorf(origYi);}
                    float fy = origYi - y0;
                    int y1 = y0 + 1;

                    int z0 = 0;
                    if (origZi > MapMaxZ){
                        z0 = MapMaxZ - 1; }
                    else{
                        z0 = floorf(origZi); }
                    float fz = origZi - z0;
                    int z1 = z0 + 1;

                    float d000 = mrcData_gpu[x0 + y0 * MapMaxX + z0 * MapMaxXY];
                    float d001 = mrcData_gpu[x1 + y0 * MapMaxX + z0 * MapMaxXY];
                    float d010 = mrcData_gpu[x0 + y1 * MapMaxX + z0 * MapMaxXY];
                    float d011 = mrcData_gpu[x1 + y1 * MapMaxX + z0 * MapMaxXY];
                    float d100 = mrcData_gpu[x0 + y0 * MapMaxX + z1 * MapMaxXY];
                    float d101 = mrcData_gpu[x1 + y0 * MapMaxX + z1 * MapMaxXY];
                    float d110 = mrcData_gpu[x0 + y1 * MapMaxX + z1 * MapMaxXY];
                    float d111 = mrcData_gpu[x1 + y1 * MapMaxX + z1 * MapMaxXY];

                    float dx00 = ((d000) + ((d001) - (d000)) * (fx));
                    float dx01 = ((d100) + ((d101) - (d100)) * (fx));
                    float dx10 = ((d010) + ((d011) - (d010)) * (fx));
                    float dx11 = ((d110) + ((d111) - (d110)) * (fx));
                    float dxy0 = ((dx00) + ((dx10) - (dx00)) * (fy));
                    float dxy1 = ((dx01) + ((dx11) - (dx01)) * (fy));

                    pixelValue = ((dxy0) + ((dxy1) - (dxy0)) * (fz));
                    }
                rotatedMapArray_gpu[pixPos] = pixelValue;
                }

              }
}
extern "C"{
__global__ void placeVolume(float *rotatedMap, float *outputMapArray, float *outputColorMapArray,
                            int rotatedMapMaxX, int rotatedMapMaxY, int rotatedMapMaxZ,
                            float shiftX, float shiftY, float shiftZ, int mapMaxX, int mapMaxY, int mapMaxZ,
                            bool outputColorMap, float colorMapTreshold, int colorMapExtend, float coloringValue, bool radial_color, float inputApix)
{

    size_t pixPos = blockIdx.x * blockDim.x + threadIdx.x;

    
    if (pixPos < (rotatedMapMaxX*rotatedMapMaxY*rotatedMapMaxZ)){
        int rotatedMapMaxXY = rotatedMapMaxX * rotatedMapMaxY;
        
        float pixel = rotatedMap[pixPos];
        
        int origZ = floorf(pixPos / (rotatedMapMaxXY));
        int origY = floorf((pixPos % (rotatedMapMaxXY)) / rotatedMapMaxY);
        int origX = (pixPos % (rotatedMapMaxXY)) % rotatedMapMaxY;

        int newZ = round(origZ - rotatedMapMaxZ / 2 + shiftZ);
        int newY = round(origY - rotatedMapMaxY / 2 + shiftY);
        int newX = round(origX - rotatedMapMaxX / 2 + shiftX);

        // data[x][y][z] = data[x + ymaxX + zmaxX*maxY]
        if ((newZ < mapMaxZ) && (newY < mapMaxY) && (newX < mapMaxX) && (newZ >= 0) && (newY >= 0) && (newX >= 0)){
            if (outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] == 0){
                outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] = pixel; }
            else if (outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] < pixel){
                outputMapArray[newX + newY * mapMaxX + newZ * mapMaxX * mapMaxY] = pixel; }
        }

        if (outputColorMap) {
            if (pixel >= colorMapTreshold) {
                for (int x_extend = (newX - colorMapExtend); x_extend < (newX + colorMapExtend); x_extend++) {
                    for (int y_extend = (newY - colorMapExtend); y_extend < (newY + colorMapExtend); y_extend++) {
                        for (int z_extend = (newZ - colorMapExtend); z_extend < (newZ + colorMapExtend); z_extend++) {
                            if ((z_extend < mapMaxZ) && (y_extend < mapMaxY) && (x_extend < mapMaxX) &&
                                (z_extend >= 0) &&
                                (y_extend >= 0) && (x_extend >= 0)) {
                                if (radial_color) {
                                   coloringValue = sqrt((x_extend-shiftX)*(x_extend-shiftX) + (y_extend-shiftY)*(y_extend-shiftY) + (z_extend-shiftZ)*(z_extend-shiftZ))*inputApix; }
                                outputColorMapArray[x_extend + y_extend * mapMaxX +
                                                    z_extend * mapMaxX * mapMaxY] = coloringValue;
                            }
                        }
                    }
                }
            }
        }
    }
}
}
                """

def readMrcSizeApix(mrcFileName):
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
    with open(mrcFileName, "rb") as mrcFile:
        imageSizeX = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeY = int(struct.unpack('i', mrcFile.read(4))[0])
        imageSizeZ = int(struct.unpack('i', mrcFile.read(4))[0])
        mrcData = np.fromfile(mrcFile, dtype=np.dtype(np.float32), count=(imageSizeX * imageSizeY * imageSizeZ),
                              offset=1024 - 12)
    return mrcData

def writeMrcFile(mrcData, stencilFile, outFile):
    with open(stencilFile, "rb") as mrcStencilFile:
        mrcHeader = mrcStencilFile.read(1024)
    with open(outFile, 'wb+') as mrcFile:
        mrcFile.write(mrcHeader)
        mrcFile.seek(12, 0)
        mrcFile.write(b"\x02\x00")
        mrcFile.seek(1024, 0)
        mrcData.astype('float32').tofile(mrcFile)

def matrix_from_euler(rot, tilt, psi):
    """create a rotation matrix from three Euler anges in ZYZ convention"""
    m = np.full((3, 3), 0.0, dtype=np.float32)

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

def placeSubvolumes_gpu(inputStarFile, inputVolumeToPlace, outputMapStencil, outputPrefix, outputCmm, tomoName, starApix,
                    binning, placePartialVolumes, recenter, coloringLabel, outputColorMap, colorMapTreshold,
                    colorMapExtend, radial_color, Xtilt, Ytilt):

    total_start = time.perf_counter()

    #read in star file
    md = MetaData(inputStarFile)
    particles = []

    for particle in md:
        particles.append(particle)

    #get unbinned apix from star file
    if starApix == 0:
        if hasattr(md, "data_optics"):
            apix = float(md.data_optics[0].rlnTomoTiltSeriesPixelSize)
            print("Apix of star got form the first optic group. Apix = %f0.2" % starApix)
        elif hasattr(md, "data_"):
            apix = float(md.data_[0].rlnDetectorPixelSize)
            print(
                "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % starApix)
        elif hasattr(md, "data_particles"):
            apix = float(md.data_particles[0].rlnDetectorPixelSize)
            print(
                "No optic groups in star file. Apix of particles star got form the first particle rlnDetectorPixelSize. Apix = %0.2f" % starApix)
        else:
            print("Could not get the apix of the particles from the star file. Define it using --star_apix parameter.")
            exit()
    else:
        apix = starApix

    if tomoName != "":
        filteredParticles = []
        for particle in particles:
            if particle.rlnTomoName == tomoName:
                filteredParticles.append(particle)
        print("%s particles included in selection from %s tomogram." % (len(filteredParticles), tomoName))
        particles = filteredParticles

    #get MRC size

    inputMapMaxX, inputMapMaxY, inputMapMaxZ, inputApix, inputMapOriginX, inputMapOriginY, inputMapOriginZ = readMrcSizeApix(
        inputVolumeToPlace)

    if inputMapMaxX == 0 or inputMapMaxY == 0 or inputMapMaxZ == 0:
        print("Error: One of the dimensions of the input volume to place is 0!")
        exit()

    mapMaxX, mapMaxY, mapMaxZ, mapApix, mapOriginX, mapOriginY, mapOriginZ = readMrcSizeApix(outputMapStencil)

    print("Output map size [x, y, z]: %i x %i x %i" % (mapMaxX, mapMaxY, mapMaxZ))

    #get binnig if not defined by user
    if binning == 0.0:
        coordinateBinningFactor = mapApix / apix
        print(
                "Binning factor not provided. Calculated from stencil header and star file => Binning = %0.2f" % coordinateBinningFactor)
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
            "WARNING: Re-centering was enable but no rlnOriginXAngst/rlnOriginYAngst/rlnOriginZAngst in star file. Disabling recentering!")
        recenter = False

    #initialize outputMapArray, outputColorMapArray, mapToRotate array
    outputMapArray = np.full((mapMaxX * mapMaxY * mapMaxZ), 0.0, dtype=np.float32)
    outputColorMapArray = np.full((mapMaxX * mapMaxY * mapMaxZ), 0.0, dtype=np.float32)
    mapToRotate = readMrcData(inputVolumeToPlace)

    nrOfParticles = len(particles)
    currentParticleNr = 1

    if tomoName != "":
        outCmmItems = cmmData(tomoName)
    else:
        outCmmItems = cmmData(inputStarFile[:-5])

    # Initialize CUDA
    checkCudaErrors(cuda.cuInit(0));
    cuDevice = findCudaDeviceDRV()
    # Create context
    cuContext = checkCudaErrors(cuda.cuCtxCreate(0, cuDevice))

    uvaSupported = checkCudaErrors(cuda.cuDeviceGetAttribute(cuda.CUdevice_attribute.CU_DEVICE_ATTRIBUTE_UNIFIED_ADDRESSING, cuDevice))
    if not uvaSupported:
        print("Accessing pageable memory directly requires UVA")
        return

    kernelHelper = common.KernelHelper(placeBackKernels, int(cuDevice))
    _rotate_kernel = kernelHelper.getFunction(b'rotate')
    _placeVolume_kernel = kernelHelper.getFunction(b'placeVolume')

    # allocate memory for mapToRotate, rotatedMap, rotMatrix, outputMapArray, outputColorMapArray
    d_mapToRotate = checkCudaErrors(cuda.cuMemAlloc(mapToRotate.nbytes))
    d_rotatedMap = checkCudaErrors(cuda.cuMemAlloc(mapToRotate.nbytes))
    d_rotMatrix = checkCudaErrors(cuda.cuMemAlloc(9*4))
    d_outputMapArray = checkCudaErrors(cuda.cuMemAlloc(outputMapArray.nbytes))
    d_outputColorMapArray = checkCudaErrors(cuda.cuMemAlloc(outputColorMapArray.nbytes))

    # Copy vectors from host memory to device memory
    checkCudaErrors(cuda.cuMemcpyHtoD(d_mapToRotate, mapToRotate, mapToRotate.nbytes))
    checkCudaErrors(cuda.cuMemcpyHtoD(d_outputMapArray, outputMapArray, outputMapArray.nbytes))
    checkCudaErrors(cuda.cuMemcpyHtoD(d_outputColorMapArray, outputColorMapArray, outputColorMapArray.nbytes))

    NUM_THREADS = 1024 # Threads per block
    NUM_BLOCKS = int(ceil(len(mapToRotate) / NUM_THREADS))  # Blocks per grid

    print("Placing %i subvolumes:" % nrOfParticles)
    for particle in particles:

        # a simple progress bar
        sys.stdout.write('\r')
        progress = int(currentParticleNr / nrOfParticles * 20)
        sys.stdout.write("[%-20s] %d%%" % ('=' * progress, 5 * progress))
        sys.stdout.flush()

        # convert Relion5 centered coordinates
        if hasattr(particle, "rlnCenteredCoordinateXAngst"):
            setattr(particle, "rlnCoordinateX",
                    (particle.rlnCenteredCoordinateXAngst + mapMaxX / 2 * mapApix) / apix)
            setattr(particle, "rlnCoordinateY",
                    (particle.rlnCenteredCoordinateYAngst + mapMaxY / 2 * mapApix) / apix)
            setattr(particle, "rlnCoordinateZ",
                    (particle.rlnCenteredCoordinateZAngst + mapMaxZ / 2 * mapApix) / apix)

        if recenter:
            xpos = (particle.rlnCoordinateX - particle.rlnOriginXAngst / apix) / coordinateBinningFactor
            ypos = (particle.rlnCoordinateY - particle.rlnOriginYAngst / apix) / coordinateBinningFactor
            zpos = (particle.rlnCoordinateZ - particle.rlnOriginZAngst / apix) / coordinateBinningFactor
        else:
            xpos = particle.rlnCoordinateX / coordinateBinningFactor
            ypos = particle.rlnCoordinateY / coordinateBinningFactor
            zpos = particle.rlnCoordinateZ / coordinateBinningFactor

        if Xtilt != 0 or Ytilt != 0:
            dx = xpos - mapMaxX / 2
            dy = ypos - mapMaxY / 2
            dz = zpos - mapMaxZ / 2

            rot_m = matrix3.matrix_from_euler_xyz(radians(Xtilt), radians(Ytilt), 0)
            xpos = mapMaxX / 2 + dx * rot_m.m[0][0] + dy * rot_m.m[0][1] + dz * rot_m.m[0][2]
            ypos = mapMaxY / 2 + dx * rot_m.m[1][0] + dy * rot_m.m[1][1] + dz * rot_m.m[1][2]
            zpos = mapMaxZ / 2 + dx * rot_m.m[2][0] + dy * rot_m.m[2][1] + dz * rot_m.m[2][2]

            particle_rot_matrix = matrix3.matrix_from_euler(radians(particle.rlnAngleRot),
                                                            radians(particle.rlnAngleTilt),
                                                            radians(particle.rlnAnglePsi))

            rotated_particle_matrix = matrix3.matrix_multiply(rot_m, particle_rot_matrix)

            angleRot, angleTilt, anglePsi = [angle_in_radian * 180.0 / pi for angle_in_radian in
                                             euler_from_matrix(rotated_particle_matrix)]
        else:
            angleRot = particle.rlnAngleRot
            angleTilt = particle.rlnAngleTilt
            anglePsi = particle.rlnAnglePsi

        ###### kernel to rotate volume
        rotMatrixInv = np.linalg.inv(matrix_from_euler(radians(angleRot), radians(angleTilt), radians(anglePsi)))
        rotMatrixInv = rotMatrixInv.astype(dtype=np.dtype(np.float32))
        checkCudaErrors(cuda.cuMemcpyHtoD(d_rotMatrix, rotMatrixInv, rotMatrixInv.nbytes))

        kernelArgs = ((d_mapToRotate, d_rotatedMap, d_rotMatrix, 0.0, 0.0, 0.0, inputMapMaxX, inputMapMaxY, inputMapMaxZ),
                      (None, None, None, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, ctypes.c_int, ctypes.c_int))

        checkCudaErrors(cuda.cuLaunchKernel(_rotate_kernel,
                                            NUM_BLOCKS, 1, 1,
                                            NUM_THREADS, 1, 1,
                                            0, 0,
                                            kernelArgs, 0))


        if ((round(xpos) + inputMapMaxX / 2) < mapMaxX and (round(ypos) + inputMapMaxY / 2) < mapMaxY and (
                round(zpos) + inputMapMaxZ / 2) < mapMaxZ and (round(xpos) - inputMapMaxX / 2) >= 0 and (
                    round(ypos) - inputMapMaxY / 2) >= 0 and (
                    round(zpos) - inputMapMaxZ / 2) >= 0) or placePartialVolumes:

            if coloringLabel != "":
                coloringValue = round(float(getattr(particle, coloringLabel)) * 100) / 100
            else:
                coloringValue = 0

            if outputCmm:
                # to keep the cmm selection possibility wee need in any case the original coords from the star
                origCoordX = particle.rlnCoordinateX
                origCoordY = particle.rlnCoordinateY
                origCoordZ = particle.rlnCoordinateZ

                if recenter:
                    particle.rlnCoordinateX -= particle.rlnOriginXAngst / apix
                    particle.rlnCoordinateY -= particle.rlnOriginYAngst / apix
                    particle.rlnCoordinateZ -= particle.rlnOriginZAngst / apix

                if Xtilt != 0 or Ytilt != 0:
                    dx = particle.rlnCoordinateX - mapMaxX * coordinateBinningFactor / 2
                    dy = particle.rlnCoordinateY - mapMaxY * coordinateBinningFactor / 2
                    dz = particle.rlnCoordinateZ - mapMaxZ * coordinateBinningFactor / 2

                    rot_m = matrix3.matrix_from_euler_xyz(radians(Xtilt), radians(Ytilt), 0)
                    particle.rlnCoordinateX = mapMaxX * coordinateBinningFactor / 2 + dx * rot_m.m[0][0] + dy * \
                                              rot_m.m[0][1] + dz * rot_m.m[0][2]
                    particle.rlnCoordinateY = mapMaxY * coordinateBinningFactor / 2 + dx * rot_m.m[1][0] + dy * \
                                              rot_m.m[1][1] + dz * rot_m.m[1][2]
                    particle.rlnCoordinateZ = mapMaxZ * coordinateBinningFactor / 2 + dx * rot_m.m[2][0] + dy * \
                                              rot_m.m[2][1] + dz * rot_m.m[2][2]

                outCmmItems.addItem(particle.rlnCoordinateX * apix, particle.rlnCoordinateY * apix,
                                    particle.rlnCoordinateZ * apix,
                                    coloringValue, origCoordX, origCoordY,
                                    origCoordZ,
                                    inputMapMaxX * inputApix / 6)

            # place volume CUDA kernel
            kernelArgs = (
            (d_rotatedMap, d_outputMapArray, d_outputColorMapArray, inputMapMaxX, inputMapMaxY, inputMapMaxZ, xpos, ypos, zpos, mapMaxX, mapMaxY, mapMaxZ, outputColorMap, colorMapTreshold, colorMapExtend, coloringValue, radial_color, inputApix),
            (None, None, None, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_bool, ctypes.c_float, ctypes.c_int, ctypes.c_float, ctypes.c_bool, ctypes.c_float,))

            checkCudaErrors(cuda.cuLaunchKernel(_placeVolume_kernel,
                                                NUM_BLOCKS, 1, 1,
                                                NUM_THREADS, 1, 1,
                                                0, 0,
                                                kernelArgs, 0))
        currentParticleNr += 1

    #end the progress bar
    sys.stdout.write('\r\n')

    # copy resulting placebacked volumes back to host memory
    checkCudaErrors(cuda.cuMemcpyDtoH(outputMapArray, d_outputMapArray, outputMapArray.nbytes))
    checkCudaErrors(cuda.cuMemcpyDtoH(outputColorMapArray, d_outputColorMapArray, outputColorMapArray.nbytes))

    # free CUDA allocated memory
    checkCudaErrors(cuda.cuMemFree(d_mapToRotate))
    checkCudaErrors(cuda.cuMemFree(d_rotatedMap))
    checkCudaErrors(cuda.cuMemFree(d_rotMatrix))
    checkCudaErrors(cuda.cuMemFree(d_outputMapArray))
    checkCudaErrors(cuda.cuMemFree(d_outputColorMapArray))

    checkCudaErrors(cuda.cuCtxDestroy(cuContext))

    print("Writing MRC file...")
    writeMrcFile(outputMapArray, outputMapStencil, outputPrefix + ".mrc")
    print("%s.mrc written." % outputPrefix)

    if outputColorMap:
        writeMrcFile(outputColorMapArray, outputMapStencil, outputPrefix + "_color.mrc")
        print("%s_color.mrc written." % outputPrefix)

    if outputCmm:
        outCmmItems.writeCmmFile(outputPrefix + ".cmm")
        print("%s.cmm written." % outputPrefix)

    print("===>Total run time: %0.2f sec" % ((time.perf_counter() - total_start)))
    print("All done. Have fun!")