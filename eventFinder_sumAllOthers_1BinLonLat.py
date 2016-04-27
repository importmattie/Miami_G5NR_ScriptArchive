import datetime
import netCDF4
import numpy
import sys

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
var1Bin = int(sys.argv[4])
var2Bin = int(sys.argv[5])

fileTag = 'r'+str(lonGC)+'x'+str(latGC)+'_'+str(hourRes)

startHourInS = hourRes*(3600./2.)
startDate = datetime.datetime(2005,06,01)+datetime.timedelta(0,startHourInS)

print fileTag

cdfIn = netCDF4.Dataset('/data2/niznik/histOutput_'+fileTag+'.nc4','r')
precBinned = cdfIn.variables['PRECBINNED'][:]
w500Binned = cdfIn.variables['W500BINNED'][:]
wpupBinned = cdfIn.variables['WPUPBINNED'][:]

#EDIT THESE LINES
lonBin = 6-1
latBin = 6-1
#END EDIT AREA

lonFirst = lonBin*lonBinLen
latFirst = latBin*latBinLen
lonLast = lonFirst+lonBinLen-1
latLast = latFirst+latBinLen-1

#EDIT THESE LINES TOO
var1Binned = precBinned[:,latFirst:latLast+1,lonFirst:lonLast+1]
var2Binned = wpupBinned[:,latFirst:latLast+1,lonFirst:lonLast+1]
#var2Binned = kedotBinned[:,latFirst:latLast+1,lonFirst:lonLast+1]
#var1Bin = 0
#var2Bin = 9
#END EDIT AREA

rightVar1Bin = numpy.where(var1Binned == var1Bin,0,-1)
rightVar2Bin = numpy.where(var2Binned == var2Bin,0,-1)

results = numpy.where(rightVar1Bin+rightVar2Bin == 0)

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
