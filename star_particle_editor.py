#!/usr/bin/env python3
import numpy as np
import struct
import napari
from napari.layers import Points
from magicgui import magicgui
import vispy.color
from lib.metadata import MetaData
from copy import deepcopy
import argparse

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

def main(inputStars, inputMrcs, outputStarFile, tomoName, coloringLb, binning, pointSize):

#    inputMrcs = "tomo67_1_bin8_filtered.rec,tomo67_1_bin8_convmap.mrc"
#    inputStars = "tomo67_clustered.star,tomo67_clustered_napari.star,tomo67_clustered_napari3.star"
#    outputStarFile = "tomo67_clustered_napari.star"
#    tomoName = "tomo67"
#    coloringLb = "rlnCtfFigureOfMerit"
#    binning = 0
#    pointSize = 7

    if "," in inputMrcs:
        inputMrcs = inputMrcs.replace(","," ")

    inputMrcFiles = inputMrcs.split(" ")

    if "," in inputStars:
        inputStars = inputStars.replace(",", " ")

    inputStarFiles = inputStars.split(" ")

    viewer = napari.Viewer(title='Star particle editor')

    for inputMrcFile in inputMrcFiles:
        mrcData = readMrcData(inputMrcFile)
        imageSizeX, imageSizeY, imageSizeZ, apix, originX, originY, originZ = readMrcSizeApix(inputMrcFile)
        mrcData = mrcData.reshape((-1, imageSizeY, imageSizeX))
        viewer.add_image(mrcData, name=inputMrcFile)

    # take the apix form the first mrc file (for binnig calculation)
        imageSizeX, imageSizeY, imageSizeZ, apix, originX, originY, originZ = readMrcSizeApix(inputMrcFiles[0])

    metaDatas = []
    pointLayers = []
    metaDataCounter = 0

    gredYellowGreen = "yellow"
    if coloringLb != "":
        rainbowArray = []
        rangeColor = 512
        for i in range(int(rangeColor / 2)):  # from red -> yellow
            rainbowArray.append([1, i / (rangeColor / 2), 0])
        for i in range(int(rangeColor / 2), 0, -1):  # from yellow -> green
            rainbowArray.append([i / (rangeColor / 2), 1, 0])
        gredYellowGreen = vispy.color.Colormap(rainbowArray)


    # calculate point coordinates from particles
    for inputStarFile in inputStarFiles:
        md = MetaData(inputStarFile)
        metaDatas.append(md)
        particles = []
        for particle in md:
            particles.append(particle)

        #get unbinned apix from star file
        starApix = float(md.data_optics[0].rlnTomoTiltSeriesPixelSize)

        if binning == 0: #binnig not set by user
            binning = apix / starApix

        print("For %s binning is: %0.2f" % (inputStarFile, binning))

        pointCoordinates = []
        coloringValues = []
        particleIDs = []
        particleCounter = 0
        for particle in particles:
            if tomoName == "" or tomoName == particle.rlnTomoName:
                pointCoordinates.append([particle.rlnCoordinateZ/binning, particle.rlnCoordinateY/binning, particle.rlnCoordinateX/binning])
                if coloringLb != "":
                    coloringValues.append(particle.rlnCtfFigureOfMerit)
                else:
                    coloringValues.append(0.0)
                particleIDs.append(particleCounter)
            particleCounter += 1

        if len(pointCoordinates) == 0:
            print("No particles found for tomoName %s in star file %s." % (tomoName, inputStarFile))
            viewer.text_overlay.text = "No particles found for tomoName %s in star file %s." % (tomoName, inputStarFile)
            viewer.text_overlay.visible = True
            continue

        coloringValues = np.array(coloringValues)
        particleIDs = np.array(particleIDs)
        metadataID = metaDataCounter

        point_features = {'particleID': particleIDs, 'metadataID': metadataID, 'coloringValue': coloringValues }
        pointLayers.append(viewer.add_points(pointCoordinates, features=point_features, size=pointSize, name=inputStarFile, face_color='coloringValue', face_colormap=gredYellowGreen, edge_colormap=gredYellowGreen, edge_color='coloringValue', out_of_slice_display=True))
        pointLayers[metadataID].feature_defaults['particleID'] = "NaN"
        pointLayers[metadataID].feature_defaults['metadataID'] = metadataID

        metaDataCounter += 1

    @magicgui(pointsToSave = {'label': 'Save:'}, outputFile = {'label': 'Star file name:'}, call_button='Save')
    def saveStarFile(pointsToSave: Points, outputFile=outputStarFile):
        pointLayerID = pointLayers.index(pointsToSave)

        newParticles = []
        nonEditedParticles = []
        pointCounter = 0
        for point in pointLayers[pointLayerID].data:
            if pointLayers[pointLayerID].features['particleID'][pointCounter] != "NaN":
                particleID = int(pointLayers[pointLayerID].features['particleID'][pointCounter])
            else:
                # take the first particle from the required tomoName as template
                particleID = int(pointLayers[pointLayerID].features['particleID'][0])

            metadataID = int(pointLayers[pointLayerID].features['metadataID'][pointCounter])

            newParticle = metaDatas[metadataID].data_particles[particleID]
            if pointLayers[pointLayerID].features['particleID'][pointCounter] == "NaN":
                newParticle.rlnAngleRot = 0.0
                newParticle.rlnAngleTilt = 0.0
                newParticle.rlnAnglePsi = 0.0
            # only change the coordinates in the star if the difference is more than a pixel (avoid change by rounding errors)
            if abs(newParticle.rlnCoordinateX - point[2] * binning) > 1:
                newParticle.rlnCoordinateX = point[2] * binning
            if abs(newParticle.rlnCoordinateY - point[1] * binning) > 1:
                newParticle.rlnCoordinateY = point[1] * binning
            if abs(newParticle.rlnCoordinateZ - point[0] * binning) > 1:
                newParticle.rlnCoordinateZ = point[0] * binning
            newParticles.append(deepcopy(newParticle))
            pointCounter += 1

        mdOut = metaDatas[metadataID].clone()
        dataTableName = "data_particles"

        mdOut.removeDataTable(dataTableName)
        mdOut.addDataTable(dataTableName, metaDatas[metadataID].isLoop(dataTableName))
        mdOut.addLabels(dataTableName, metaDatas[metadataID].getLabels(dataTableName))
        # add the records of non-selected tomoName
        for particle in metaDatas[metadataID].data_particles:
            if particle.rlnTomoName != tomoName:
                nonEditedParticles.append(particle)
        mdOut.addData(dataTableName, nonEditedParticles)
        mdOut.addData(dataTableName, newParticles)
        mdOut.write(outputFile)
        viewer.text_overlay.text = "Output star file %s written." % outputFile
        viewer.text_overlay.visible = True
        print("Output star file %s written." % outputFile)

    @magicgui(copyFrom={'label': 'Copy from:'}, copyTo={'label': 'Copy to:'}, call_button='Copy')
    def copyPointsBetweenLayers(copyFrom: Points, copyTo: Points):
        if copyFrom == copyTo:
            print("Cannot copy between the same layers!")
            viewer.text_overlay.text = "Cannot copy between the same layers!"
            viewer.text_overlay.visible = True
            return

        copyFromLayerID = pointLayers.index(copyFrom)
        copyToLayerID = pointLayers.index(copyTo)

        for point in pointLayers[copyFromLayerID].selected_data:
            pointLayers[copyToLayerID].feature_defaults['particleID'] =  pointLayers[copyFromLayerID].features['particleID'][point]
            pointLayers[copyToLayerID].feature_defaults['metadataID'] =  pointLayers[copyFromLayerID].features['metadataID'][point]
            pointLayers[copyToLayerID].feature_defaults['coloringValue'] =  pointLayers[copyFromLayerID].features['coloringValue'][point]
            pointLayers[copyToLayerID].add(pointLayers[copyFromLayerID].data[point])
            pointLayers[copyToLayerID].refresh_colors()

        # set back the feature default of the copyTo layer
        pointLayers[copyToLayerID].feature_defaults['particleID'] = 'NaN'
        pointLayers[copyToLayerID].feature_defaults['metadataID'] = copyToLayerID

        viewer.text_overlay.text = "%s points copied from layer %s to layer %s." % (len(pointLayers[copyFromLayerID].selected_data),pointLayers[copyFromLayerID].name,pointLayers[copyToLayerID].name)
        viewer.text_overlay.visible = True

    @magicgui( call_button='Color IMOD style')
    def imodStylePoints():
        for pointLayer in pointLayers:
            pointLayer.face_color = [0,0,0,0]
            pointLayer.edge_width = 0.05
            pointLayer.edge_width_is_relative = True
            pointLayer.edge_colormap = None
            pointLayer.edge_color = 'lime'
            pointLayer.refresh_colors()

    viewer.window.add_dock_widget(copyPointsBetweenLayers, name="Copy Points Between Layers")
    viewer.window.add_dock_widget(saveStarFile, name="Save Star File")
    viewer.window.add_dock_widget(imodStylePoints, name="Color IMOD style")
    napari.run()


#    if napari_view:
#        def next_on_click(layer, event):
#            """Mouse click binding to advance the label when a point is added"""
#            if layer.mode == 'add':
#                print(pointsLayer.features['particle_data'])

#        green = vispy.color.Colormap([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Visual editor to add, remove, combine particles from tomo STAR files.")
    add = parser.add_argument
    add('--i', help="Input star file. Possible to open multiple star files delimited by comma.")
    add('--itomo', help="Input mrc file. Possible to open multiple mrc files delimited by comma.")
    add('--o', help="Output star file name. Possible to change later in the gui.")
    add('--tomo_name', type=str, default="",
        help="Use only particles from tomogram equal in rlnTomoName. OPTIONAL: If not set all particles will be visualized.")
    add('--color_lb', type=str, default="",
        help="Label from the star file that will be used for rainbow coloring of the markers.")
    add('--bin', type=float, default="0.0",
        help="User provided binnig of the --itomo. If not set the apix of the first MRC file in --itomo is taken to calcualte the binning.")
    add('--point_size', type=int, default="7",
        help="Size of the points in pixels.")

    args = parser.parse_args()

    #if args.i == "" or args.o == "" or len(sys.argv) == 1:
    #    parser.print_help()
    #    exit()
    main(args.i, args.itomo, args.o, args.tomo_name, args.color_lb, args.bin, args.point_size)
    #main()