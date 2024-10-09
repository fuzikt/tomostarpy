#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import argparse
import sys
from lib.newstlist import *

def removeElementsByIndex(parent, indices):
    """Remove lines from the parent node's text by line indices."""
    # Split the parent text by lines
    lines = parent.text.strip().splitlines()

    # Filter out lines that are not in the indices_to_remove
    filteredLines = [line for idx, line in enumerate(lines) if idx not in indices]

    # Update the parent text by joining the remaining lines
    parent.text = "\n" + "\n".join(filteredLines) + "\n"

def removeAndRenumberTags(parent, tag, indices):
    """Remove elements based on indices and renumber remaining elements."""
    elements_to_keep = []
    # Iterate through child elements with the specified tag and remove the ones at specified indices
    for idx, element in enumerate(parent.findall(f"./{tag}")):
        if int(element.attrib['ID']) in indices:
            parent.remove(element)
    # Now renumber the kept elements starting from 0
    for idx, element in enumerate(parent.findall(f"./{tag}")):
        element.attrib['ID'] = str(idx)


def removeNodesByZvalue(parent, tag, indices):
    """Remove <Node> elements whose 'Z' attribute matches the indices and renumber the remaining 'Z' values."""
    # Find all <Node> elements inside the given tag
    preservedZvalues = []
    for node in parent.findall(f".//{tag}/Node"):
        zValue = int(node.attrib.get("Z", -1))  # Get the Z attribute, default to -1 if missing
        if zValue in indices:
            # tree.find(f"./{tag}").remove(node)
            parent.find(f"./{tag}").remove(node)
        else:
            preservedZvalues.append(zValue)

    # Now renumber the remaining <Node> elements so that Z starts from 0
    uniquePreservedZvalues = list(set(preservedZvalues))
    for node in parent.findall(f".//{tag}/Node"):
        node.attrib['Z'] = str(uniquePreservedZvalues.index(int(node.attrib.get("Z", -1))))

def excludeTiltsFromTomoStar(inputFile, outputFile, indicesToRemove):
    with open(inputFile, 'r') as tomoStarFile:
        tomoStarLines =  tomoStarFile.readlines()
    tiltCounter = 0
    filteredTomoStarLines = []

    for line in tomoStarLines:
        if line.startswith("_wrp") or line.startswith('\n') or line.startswith('\r\n') or line.startswith('data_') or line.startswith('loop_'):
            filteredTomoStarLines.append(line)
        else:
            if tiltCounter not in indicesToRemove:
                filteredTomoStarLines.append(line)
            tiltCounter +=1

    with open(outputFile, 'w') as outTiltListFile:
        for line in filteredTomoStarLines:
            outTiltListFile.write("%s" % line)

def main(inputFile, outputPrefix, tomostarDir, excludeFile):
    with open(excludeFile, 'r') as tiltListFile:
        rangeList = tiltListFile.read()

    if len(rangeList.strip()) == 0:
        print("No tilts listed to be excluded...")
        exit()

    indicesToRemove = extractMembersOfRangeList(rangeList)
    print("Tilts that are going to be excluded: %s" % indicesToRemove)
    # Convert 1-based indices to 0-based for Python
    indicesToRemove = [i - 1 for i in indicesToRemove]

    # Load the XML file
    tree = ET.parse(inputFile)
    root = tree.getroot()

    # Define the tags to process for removing text members
    tagsToProcess = ["Angles", "Dose", "UseTilt", "AxisAngle", "AxisOffsetX", "AxisOffsetY", "MoviePath"]

    # Define the tags to process for removing entire elements and renumbering the remaining elements
    tagsToRenumber = ["TiltPS1D", "TiltSimulatedScale"]

    # Define the tags containing <Node> elements where elements are removed according to Z attribute value and nodes are renumbered
    tagsWithNodes = [
        "GridCTF", "GridCTFDefocusDelta", "GridCTFDefocusAngle", "GridCTFPhase",
        "GridMovementX", "GridMovementY", "GridAngleX", "GridAngleY", "GridAngleZ",
        "GridDoseBfacs", "GridDoseBfacsDelta", "GridDoseBfacsAngle", "GridDoseWeights"
    ]

    for tag in tagsToProcess:
        # Find the element for the current tag
        element = root.find(f".//{tag}")

        if element is not None and element.text is not None:
            # Remove the elements at the specified indices
            removeElementsByIndex(element, indicesToRemove)

    for tag in tagsToRenumber:
        parent = root.find(".")  # Root element to modify the list of <TiltPS1D> or <TiltSimulatedScale> elements
        if parent is not None:
            removeAndRenumberTags(parent, tag, indicesToRemove)

    for tag in tagsWithNodes:
        parent = root.find(".")  # Root element
        if parent is not None:
            removeNodesByZvalue(parent, tag, indicesToRemove)

    # Save the modified XML to a new file
    tree.write(outputPrefix+".xml", encoding='utf-8', xml_declaration=True)

    print("XML has been modified and saved to %s" % outputPrefix)

    tomoStarFileName=tomostarDir+"/"+args.i.replace("xml", "tomostar")
    outTomoStarFileName=tomostarDir+"/"+args.o+".tomostar"
    excludeTiltsFromTomoStar(tomoStarFileName, outTomoStarFileName, indicesToRemove)

    print("Tomostar has been modified and saved to %s" % outTomoStarFileName)

    print("All done! Have fun!")


if __name__ == "__main__":
    # if called directly from command line
    parser = argparse.ArgumentParser(
        description="Removes tilts defined in --exclude_file from Warp2 tilt-series XML file and tomostar file.")
    add = parser.add_argument
    add('--i', help="Input Warp2 tilt-series XML file.")
    add('--o', help="Output prefix for XML and tomostar file.")
    add('--tomostar_dir', type=str, default="./",
        help="Directory where tomostar files are located. (Default ./)")
    add('--exclude_file', type=str, default="",
        help="Textfile with comma separated list of excluded views (or range of views as in newstack).")

    args = parser.parse_args()

    if args.i == "" or args.o == "" or len(sys.argv) == 1:
        parser.print_help()
        exit()

    main(args.i, args.o, args.tomostar_dir, args.exclude_file)
