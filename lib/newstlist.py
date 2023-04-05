# converts newstacks compatible list of included/excluded views to a python list of separate views
# e.g. "1-6,10,44-46"  -> [1, 2, 3, 4, 5, 6, 10, 44, 45, 46]
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

# creates a newstacks compatible list of included/excluded views from a python list
# e.g. [1, 2, 3, 4, 5, 6, 10, 44, 45, 46] -> "1-6,10,44-46"
def rangeListCreate(excludeList):
    rangeList = ""
    rangeStart = False

    # sort and unique
    excludeList = sorted(set(excludeList))

    for itemNr in range(len(excludeList)):
        if itemNr == 0:
            rangeList += str(excludeList[itemNr])
            continue
        if excludeList[itemNr] - 1 == excludeList[itemNr - 1]:
            rangeStart = True
        else:
            if rangeStart:
                rangeStart = False
                rangeList += "-" + str(excludeList[itemNr - 1]) + "," + str(excludeList[itemNr])
            else:
                rangeList += "," + str(excludeList[itemNr])

    if rangeStart:
        rangeList += "-" + str(excludeList[itemNr])

    return rangeList