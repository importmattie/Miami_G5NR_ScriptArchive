import datetime
import netCDF4
import numpy
import sys
import time

#DO NOT TOUCH - USER INPUT AND STATIC ENTRY
lonGC = int(sys.argv[1])
latGC = int(sys.argv[2])
lonBinLen = lonGC/10
latBinLen = latGC/9
lonRes = 360/lonGC
latRes = 180/latGC
hourRes = int(sys.argv[3])
lonStart = 0+lonRes/2.
latStart = -90+latRes/2.
var2Bin = int(sys.argv[4])
var3Bin = int(sys.argv[5])
var4Bin = int(sys.argv[6])
var5Bin = int(sys.argv[7])
var6Bin = int(sys.argv[8])
var7Bin = int(sys.argv[9])
var8Bin = int(sys.argv[10])
var9Bin = int(sys.argv[11])
outputForBundleScript = int(sys.argv[12])

fileTag = 'r'+str(lonGC)+'x'+str(latGC)+'_'+str(hourRes)

startHourInS = hourRes*(3600./2.)
startDate = datetime.datetime(2005,06,01)+datetime.timedelta(0,startHourInS)

lonBin = var2Bin
latBin = var3Bin
lonFirst = lonBin*lonBinLen
latFirst = latBin*latBinLen
lonLast = lonFirst+lonBinLen-1
latLast = latFirst+latBinLen-1

cdfIn = netCDF4.Dataset('/data2/niznik/histOutput_'+fileTag+'_V2.nc4','r')
precBinned = cdfIn.variables['PRECBINNED'][:,latFirst:latLast+1,lonFirst:lonLast+1]
w500Binned = cdfIn.variables['W500BINNED'][:,latFirst:latLast+1,lonFirst:lonLast+1]
wPuPBinned = cdfIn.variables['WPUPBINNED'][:,latFirst:latLast+1,lonFirst:lonLast+1]
TEFBinned = cdfIn.variables['TEFBINNED'][:,latFirst:latLast+1,lonFirst:lonLast+1]
KEDotBinned = cdfIn.variables['KEDOTBINNED'][:,latFirst:latLast+1,lonFirst:lonLast+1]
HMVBinned = cdfIn.variables['HMVBINNED'][:,latFirst:latLast+1,lonFirst:lonLast+1]

matchesAllBins = numpy.zeros(precBinned.shape)

if(0 <= var4Bin <= 11):
    matchesAllBins += numpy.where(precBinned == var4Bin,0,-1)
if(0 <= var5Bin <= 11):
    matchesAllBins += numpy.where(w500Binned == var5Bin,0,-1)
if(0 <= var6Bin <= 11):
    matchesAllBins += numpy.where(wPuPBinned == var6Bin,0,-1)
if(0 <= var7Bin <= 11):
    matchesAllBins += numpy.where(TEFBinned == var7Bin,0,-1)
if(0 <= var8Bin <= 11):
    matchesAllBins += numpy.where(KEDotBinned == var8Bin,0,-1)
if(0 <= var9Bin <= 11):
    matchesAllBins += numpy.where(HMVBinned == var9Bin,0,-1)

results = numpy.where(matchesAllBins == 0)

if(outputForBundleScript != 1):
    print 'Found '+str(len(results[0]))+ ' results!'

    probableDuplicates = 0
    for ii in range(0,len(results[0])):
        if(ii != 0 and abs(results[0][ii]-results[0][ii-1]) < 2 and
           abs(results[1][ii]-results[1][ii-1]) < 2 and
           abs(results[2][ii]-results[2][ii-1]) < 2):
            probableDuplicates += 1
        else:
            currentResultTime = startDate+datetime.timedelta(0,(hourRes*3600.)*results[0][ii])
            currentLon = lonStart+lonRes*(lonFirst+results[2][ii])
            if(currentLon > 180):
                currentLon = currentLon-360
            currentLat = latStart+latRes*(latFirst+results[1][ii])
            print (str("%3i" % ii)+' '+str(currentResultTime.year)+'-'+str(currentResultTime.month)+
                   '-'+str(currentResultTime.day)+' '+str(currentResultTime.hour)+':'+
                   str(currentResultTime.minute)+' '+str(currentLon)+' '+str(currentLat))

    print '(Excluded '+str(probableDuplicates)+' probable duplicates)'
else:
    probableDuplicates = 0
    for ii in range(0,len(results[0])):
        if(ii != 0 and abs(results[0][ii]-results[0][ii-1]) < 2 and
           abs(results[1][ii]-results[1][ii-1]) < 2 and
           abs(results[2][ii]-results[2][ii-1]) < 2):
            probableDuplicates += 1
        else:
            currentResultTime = time.mktime((startDate+datetime.timedelta(0,(hourRes*3600.)*results[0][ii])).timetuple())
            currentLon = lonStart+lonRes*(lonFirst+results[2][ii])
            if(currentLon > 180):
                currentLon = currentLon-360
            currentLat = latStart+latRes*(latFirst+results[1][ii])
            print (str(int(currentResultTime))+' '+str(currentLon)+' '+str(currentLat))
            
