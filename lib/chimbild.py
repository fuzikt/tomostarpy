class bildArrowItem:
    def __init__(self, x=0, y=0, z=0, value=0, shiftX=0, shiftY=0, shiftZ=0, thickness=5, r=1, g=1, b=0):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.g = g
        self.b = b
        self.thickness = thickness
        self.shiftX = shiftX
        self.shiftY = shiftY
        self.shiftZ = shiftZ
        self.value = value


class bildArrowData:

    def __init__(self):
        self.bildArrowItems = []

    def addItem(self, x, y, z, value=0, shiftX=0, shiftY=0, shiftZ=0, thickness=5, r=1, g=1, b=0):
        self.bildArrowItems.append(bildArrowItem(x, y, z, value, shiftX, shiftY, shiftZ, thickness, r, g, b))

    def createColorGradient(self):
        rainbowArray = []
        rangeColor = 512
        for i in range(int(rangeColor / 2)):  # from red -> yellow
            rainbowArray.append([1, i / (rangeColor / 2), 0])
        for i in range(int(rangeColor / 2), 0, -1):  # from yellow -> green
            rainbowArray.append([i / (rangeColor / 2), 1, 0])

        # get the min and max value (ie range) of the desired label in the star file
        minColor = self.bildArrowItems[0].value
        maxColor = self.bildArrowItems[0].value
        for item in self.bildArrowItems:
            minColor = min(minColor, item.value)
            maxColor = max(maxColor, item.value)

        for item in self.bildArrowItems:
            if minColor != maxColor:
                coloringArrayElement = rainbowArray[
                    int((item.value - minColor) / (maxColor - minColor) * (rangeColor - 1))]
            else:
                coloringArrayElement = [1, 1, 0]
            self.bildArrowItems[self.bildArrowItems.index(item)].r = coloringArrayElement[0]
            self.bildArrowItems[self.bildArrowItems.index(item)].g = coloringArrayElement[1]
            self.bildArrowItems[self.bildArrowItems.index(item)].b = coloringArrayElement[2]

    def writeBuilsArrowFile(self, outputBuildArrowFile):
        self.createColorGradient()

        with open(outputBuildArrowFile, 'w') as bildArrowFile:
            for item in self.bildArrowItems:
                bildArrowFile.write(".color %0.3f %0.3f %0.3f\n" % (item.r, item.g, item.b))
                bildArrowFile.write(".arrow %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n" % (
                item.shiftX, item.shiftY, item.shiftZ, item.x + item.shiftX, item.y + item.shiftY, item.z + item.shiftZ,
                item.thickness))
