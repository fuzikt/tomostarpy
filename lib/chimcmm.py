class cmmItem:
    def __init__(self, x=0, y=0, z=0, value=0, coordX=0, coordY=0, coordZ=0, radius=5, r=1, g=1, b=0):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.g = g
        self.b = b
        self.radius = radius
        self.coordX = coordX
        self.coordY = coordY
        self.coordZ = coordZ
        self.value = value


class cmmData:

    def __init__(self, setName=None):
        if setName == None:
            self.setName = "marker set 1"
        else:
            self.setName = setName
        self.cmmItems = []

    def addItem(self, x, y, z, value=0, coordX=0, coordY=0, coordZ=0, radius=5, r=1, g=1, b=0):
        self.cmmItems.append(cmmItem(x, y, z, value, coordX, coordY, coordZ, radius, r, g, b))

    def createColorGradient(self):
        rainbowArray = []
        rangeColor = 512
        for i in range(int(rangeColor / 2)):  # from red -> yellow
            rainbowArray.append([1, i / (rangeColor / 2), 0])
        for i in range(int(rangeColor / 2), 0, -1):  # from yellow -> green
            rainbowArray.append([i / (rangeColor / 2), 1, 0])

        # get the min and max value (ie range) of the desired label in the star file
        minColor = self.cmmItems[0].value
        maxColor = self.cmmItems[0].value
        for item in self.cmmItems:
            minColor = min(minColor, item.value)
            maxColor = max(maxColor, item.value)

        for item in self.cmmItems:
            if minColor != maxColor:
                coloringArrayElement = rainbowArray[
                    int((item.value - minColor) / (maxColor - minColor) * (rangeColor - 1))]
            else:
                coloringArrayElement = [1, 1, 0]
            self.cmmItems[self.cmmItems.index(item)].r = coloringArrayElement[0]
            self.cmmItems[self.cmmItems.index(item)].g = coloringArrayElement[1]
            self.cmmItems[self.cmmItems.index(item)].b = coloringArrayElement[2]

    def writeCmmFile(self, outputCmmFile):
        self.createColorGradient()

        with open(outputCmmFile, 'w') as cmmFile:
            cmmFile.write("<marker_set name=\"%s\">\n" % self.setName)
            cmm_id = 1
            for item in self.cmmItems:
                cmmFile.write(
                    "<marker id=\"%d\" x=\"%0.3f\" y=\"%0.3f\" z=\"%0.3f\" r=\"%0.3f\" g=\"%0.3f\" b=\"%0.3f\" radius=\"%f\" coordX=\"%0.3fpx\" coordY=\"%0.3fpx\" coordZ=\"%0.3fpx\"/> \n" % (
                        cmm_id, item.x, item.y, item.z, item.r, item.g, item.b, item.radius, item.coordX, item.coordY,
                        item.coordZ
                    ))
                cmm_id += 1
            cmmFile.write("</marker_set>")
