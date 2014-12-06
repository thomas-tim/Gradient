#!/usr/bin/env python
'''This program can calculate if you should use your aerobars and how the of you route is distributed.'''

from math import radians, cos, sin, asin, sqrt
import collections
import sys
import getopt
import urllib
import json
import gpxpy
import numpy
import matplotlib.pyplot as plt
import scipy.stats

__author__ = "Thomas Schaper"
__copyright__ = "Copyright 2014"
__credits__ = ["Thomas Schaper"]
__license__ = "GPL"
__version__ = "3.14.15"
__maintainer__ = "Thomas Schaper"
__status__ = "Production"

def encode_coords(coords):
    result = []
    prev_lat = 0
    prev_lng = 0

    for y, x in coords:
        lat, lng = int(y * 1e5), int(x * 1e5)

        d_lat = _encode_value(lat - prev_lat)
        d_lng = _encode_value(lng - prev_lng)

        prev_lat, prev_lng = lat, lng

        result.append(d_lat)
        result.append(d_lng)

    return ''.join(c for r in result for c in r)

def _split_into_chunks(value):
    while value >= 32:
        yield (value & 31) | 0x20
        value >>= 5
    yield value

def _encode_value(value):
    # Step 2 & 4
    value = ~(value << 1) if value < 0 else (value << 1)

    # Step 5 - 8
    chunks = _split_into_chunks(value)

    # Step 9-10
    return (chr(chunk + 63) for chunk in chunks)

def decode(point_str):
    # sone coordinate offset is represented by 4 to 5 binary chunks
    coord_chunks = [[]]
    for char in point_str:

        # convert each character to decimal from ascii
        value = ord(char) - 63

        # values that have a chunk following have an extra 1 on the left
        split_after = not (value & 0x20)
        value &= 0x1F

        coord_chunks[-1].append(value)

        if split_after:
                coord_chunks.append([])

    del coord_chunks[-1]

    coords = []

    for coord_chunk in coord_chunks:
        coord = 0

        for i, chunk in enumerate(coord_chunk):
            coord |= chunk << (i * 5)

        #there is a 1 on the right if the coord is negative
        if coord & 0x1:
            coord = ~coord #invert
        coord >>= 1
        coord /= 100000.0

        coords.append(coord)

    # convert the 1 dimensional list to a 2 dimensional list and offsets to
    # actual values
    points = []
    prev_x = 0
    prev_y = 0
    for i in xrange(0, len(coords) - 1, 2):
        if coords[i] == 0 and coords[i + 1] == 0:
            continue

        prev_x += coords[i + 1]
        prev_y += coords[i]
        # a round to 6 digits ensures that the floats are the same as when
        # they were encoded
        points.append((round(prev_x, 6), round(prev_y, 6)))

    return points

def distancePoints(lon1, lat1, lon2, lat2):
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))

    # 6367 km is the radius of the Earth
    m = 6367000 * c
    return int(m)

def timeExtra(grade, distance, weightTotal, weightAdded, dragReduced, power):
    if grade > 3:
        dragReduced = 0
    if grade <0:
        weightAdded = 0
    if grade <-4:
        return (1,0)
    grade = grade/100
    speedNormalFormula = numpy.poly1d([(0.3*(1) * 1.226),0,(9.8 * grade * \
                            (weightTotal)),((0.005 * weightTotal * 2.78)-power)])
    speedNormal = speedNormalFormula.r.real.max()
    speedReducedDragFormula = numpy.poly1d([(0.3*(1-dragReduced)* 1.226),0,(9.8 * grade * \
                                (weightTotal+weightAdded)),((0.005 * weightTotal+weightAdded * 2.78)-power)])
    speedReducedDrag = speedReducedDragFormula.r.real.max()
    timeNeededReducedDrag =(distance/(speedReducedDrag))
    timeNormal = ((distance)/(speedNormal))
    timeDelta = timeNormal-timeNeededReducedDrag
    return (timeNormal,timeDelta)

def addElevation(coordinatesDict):
    print("downloading elevation data from the Google Maps API")
    zzz=300
    iteration=0
    for i in range(1,len(coordinatesList)/zzz+2):
        if coordinatesList[(i-1)*zzz:i*zzz] == []:
            break
        else:
            encoded = encode_coords(coordinatesList[(i-1)*zzz:i*zzz])
            url = "http://maps.googleapis.com/maps/api/elevation/json?locations=enc:"+encoded
            data = urllib.urlopen(url).read()
            jsonData = json.loads(data)
            for result in jsonData['results']:
                coordinatesDict[iteration].append(result['elevation'])
                iteration +=1

def getGradient(list1, list2):
    distance = distancePoints(list1[0],list1[1],list2[0],list2[1])
    deltaElevation = list2[2]-list1[2]
    gradient = (deltaElevation/float(distance))*100
    return gradient

def addGradient(coordinatesList, coordinatesDict, rounding=1):
    gradeDict = {}
    for i in range(1,len(coordinatesList)):
        distance = distancePoints(coordinatesDict[i][0],coordinatesDict[i][1],coordinatesDict[i-1][0],coordinatesDict[i-1][1])
        grade = round(getGradient(coordinatesDict[i-1],coordinatesDict[i]),rounding)
        if abs(grade) <=19:
            try:
                gradeDict[grade]+=distance
            except KeyError:
                gradeDict[grade] = distance
    return gradeDict

def constructCoordinates(gpxfile):
    gpx = gpxpy.parse(open(gpxfile))
    print("file loaded")
    totalDistance=0
    coordinatesDict = {}
    coordinatesList = []
    first = True
    z =0
    elevationWasIncluded = False
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                if point.elevation != None:
                    elevationWasIncluded = True
                    if first:
                        oldLat = point.latitude
                        oldLon = point.longitude
                        oldElv = point.elevation
                        coordinatesDict[z] = [oldLon,oldLat,oldElv]
                        coordinatesList.append((oldLat,oldLon))
                        first = False
                        z+=1
                    else:
                        if distancePoints(point.longitude,point.latitude,oldLon,oldLat) >0:
                            oldLat = point.latitude
                            oldLon = point.longitude
                            oldElv = point.elevation
                            coordinatesDict[z]=[point.longitude,point.latitude,point.elevation]
                            coordinatesList.append((point.latitude,point.longitude))
                            z+=1
                else:
                    if first:
                        oldLat = point.latitude
                        oldLon = point.longitude
                        coordinatesDict[z] = [oldLon,oldLat]
                        coordinatesList.append((oldLat,oldLon))
                        first = False
                        z+=1
                    else:
                        distance = distancePoints(point.longitude,point.latitude,oldLon,oldLat)
                        if distance >= 1000 or point == segment.points[-1]:
                            totalDistance+=distance
                            oldLat = point.latitude
                            oldLon = point.longitude
                            coordinatesDict[z]=[point.longitude,point.latitude]
                            coordinatesList.append((point.latitude,point.longitude))
                            z+=1
    if elevationWasIncluded == False:
        addElevation()
    return (coordinatesDict, coordinatesList)

def parseArg(argv):
    inputFile = None
    weightTotal = None
    weightAdded = None
    dragReduced = None
    power = None
    function = None
    printString = False
    rounding = 0
    try:
        opts, args = getopt.getopt(argv,"hi:w:a:d:p:f:r:",["ifile=","ofile="])
    except getopt.GetoptError:
        print '-i <inputfile> -w <total weight including bike, excluding aerobars> -a <weight of aerobars> -d <drag reduced by aerobars> -p <power you want to use for the calculation> -f <the function, distribution/graph/timesaved> -r <rounding of the gradient, default 0> '
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print "-i <inputfile> -w <total weight including bike, excluding aerobars> -a <weight of aerobars> -d <drag reduced by aerobars> -p <power you want to use for the calculation> -f <the function, distribution/graph/timesaved> -r <rounding of the gradient, default 0>\nDon't know some data? Just run the program"
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputFile = arg
        elif opt in ("-w"):
            try:
                weightTotal = float(arg)
            except ValueError:
                print("weight has to be a number")
        elif opt in ("-a"):
            try:
                weightAdded = float(arg)
            except ValueError:
                print("weight of aerobars has to be a number")
        elif opt in ("-d"):
            try:
                dragReduced = float(arg)
            except ValueError:
                print("the ruduced drag has to be a number")
        elif opt in ("-p"):
            try:
                power = float(arg)
            except ValueError:
                print("power has to be a number")
        elif opt in ("-f"):
            function = arg
        elif opt in ("-r"):
            try:
                rounding = int(arg)
            except ValueError:
                print("rounding has to be a number")
    if None in (inputFile, weightTotal, weightAdded, dragReduced, power, function, rounding):
        printString = True
    while inputFile == None:
        try:
            inputFile = raw_input("What is the exact location of your gpx file?\n")
            gpxpy.parse(open(inputFile))
        except ValueError:
            pass
        except IOError:
            print("The file does not exist!")
            inputFile = None
    while weightTotal == None:
        try:
            weightTotal = float(raw_input("What is the weight of you and your bike and including all the things on your bike but excluding your aerobars? It has to be a number!\n"))
        except ValueError:
            pass
    while weightAdded == None:
        try:
            weightAdded = float(raw_input("What is the weight of your aerobars? It has to be a number!\n"))
        except ValueError:
            pass
    while dragReduced == None:
        dragReduced = 0.08
        try:
            dragReduced = float(raw_input("How much does it reduce your drag with aerobars? You can test this by riding a period of time with a fixed power first without aerobars and then with aerobars on a flat road. Then input the answer from this formula: 1-(speed_without_aerobars/speed_with_aerobars)^3.\nAccording to Specialized it is approximately 7%. So it defaults to 0.07 if you don't enter anything!' It has to be a number!\n"))
        except ValueError:
            pass
    while power == None:
        power = 185
        try:
            power = float(raw_input("How much power are you going to average over the segemnt you entered? (Will default to 185 Watt) It has to be a number!\n"))
        except ValueError:
            pass
    while function not in ("distribution", "graph", "timesaved"):
        function = raw_input("What function do you want to use? Please enter distribution, graph or timesaved\n")

    if printString:
        theString = sys.argv[0] + " -i " + inputFile + " -w " + str(weightTotal) + " -a" + str(weightAdded) + " -d " + str(dragReduced) \
        + " -p " + str(power) + " -f " + function + " -r " + str(rounding)
        print("Next time you can use this command: \"" + theString + "\"")
    return (inputFile, weightTotal, weightAdded, dragReduced, power, function, rounding)

def makeDistribution(gradeDict, rounding):
    print rounding
    gradeList = []
    for key in gradeDict.keys():
        for grade in range(gradeDict[key]):
            gradeList.append(key)
    gradeList.sort()
    rangeWanted =  numpy.linspace(gradeList[0],gradeList[-1], (rounding+1)*(gradeList[-1]-gradeList[0]))
    plt.hist(gradeList, rangeWanted)
    plt.show()

def makeGraph(gradeDict, weightTotal, weightAdded, dragReduced, rounding):
    startPowerRange, endPowerRange = 1,0
    try:
        startPowerRange = int(raw_input("At what power do you want to start the graph?\n"))
        endPowerRange = int(raw_input("At what power do you want to end the graph?\n"))
    except ValueError:
        print("they have to be numbers!")
    while startPowerRange > endPowerRange:
        print("the start has to be smaller than the end")
        try:
            startPowerRange = int(raw_input("At what power do you want to start the graph?\n"))
            endPowerRange = int(raw_input("At what power do you want to end the graph?\n"))
        except ValueError:
            print("they have to be numbers!")
    timeDict = {}
    plt.xlabel("Power (W)")
    plt.ylabel("Time saved (%)")
    plt.grid(True)
    for i in numpy.linspace(startPowerRange,endPowerRange,1000/((rounding*10)+1)):
        timeTotal, timeSavedTotal = 0,0
        for (grade,distance) in gradeDict.items():
            (timeAll, timeSaved) = timeExtra(grade, distance, weightTotal, weightAdded, dragReduced, i)
            timeSavedTotal+=timeSaved
            timeTotal += timeAll
        timeDict[i]=timeSavedTotal/timeTotal*100
    timeDict = collections.OrderedDict(sorted(timeDict.items()))
    plt.plot(timeDict.keys(), timeDict.values())
    plt.show()

def timeSaved(gradeDict, weightTotal, weightAdded, dragReduced, power, rounding):
    timeSavedTotal = 0
    totalDistance = 0
    timeSavedDict = {}
    totalTime = 0
    timeTotal = 0
    k=0
    gradeDict = collections.OrderedDict(sorted(gradeDict.items()))
    print("calculating time now!")
    for grade in gradeDict.keys():
        timeExtraOutput = timeExtra(grade, gradeDict[grade], weightTotal, weightAdded, dragReduced, power)
        timeSavedPercentage = float(timeExtraOutput[-1]/timeExtraOutput[0])*100
        totalDistance += gradeDict[grade]
        if timeSavedPercentage > 0:
            print("The time saved with aerobars on a grade of " +str(grade) + "% is: " +str(round(timeSavedPercentage,1)) + "% percent")
        elif timeSavedPercentage<0:
            print("The time EXTRA SPENDED with aerobars on a grade of " +str(grade) + "% is: " +str(abs(round(timeSavedPercentage,1))) + "% percent")
        elif timeSavedPercentage==0:
            print("No difference with aerobars on a grade of " +str(grade) + "%!")
        timeSavedTotal+=timeExtraOutput[-1]
        totalTime += timeExtraOutput[0]
        k += grade*gradeDict[grade]
    print("total distance:" + str(round(totalDistance/1000,3)) + "km, this is less precize than RideWithGPS!")
    print("average gradient is: " +str(round(k/totalDistance,rounding)) + "%")
    if timeSavedTotal >0:
        print("the total saved time is: " + str(round(100*(timeSavedTotal/totalTime),2)) + "% percent")
    elif timeSavedTotal <0:
        print("the total EXTRA SPENDED time is: " + str(abs(round(100*(timeSavedTotal/totalTime),2))) + "% percent")
    elif timeSavedTotal ==0:
        print("No time saved or extra SPEND. Better not do it!")

def main(argv):
    (inputFile, weightTotal, weightAdded, dragReduced, power, function, rounding) = parseArg(argv)
    (coordinatesDict, coordinatesList) = constructCoordinates(inputFile)

    gradeDict = addGradient(coordinatesList, coordinatesDict, rounding)

    if function == "distribution":
        makeDistribution(gradeDict, rounding)
    elif function == "timesaved":
        timeSaved(gradeDict,weightTotal,weightAdded,dragReduced,power,rounding)
    elif function == "graph":
        makeGraph(gradeDict, weightTotal, weightAdded, dragReduced, rounding)
    exit()

if __name__ =="__main__":
    main(sys.argv[1:])
