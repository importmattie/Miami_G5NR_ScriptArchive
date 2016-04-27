import datetime
import netCDF4
import numpy

cdfIn = netCDF4.Dataset('histOutput_r90x45_3.nc4','r')

startHourInS = hourRes*(3600./2.)
startDate = datetime.datetime(2005,06,01)+datetime.timedelta(0,startHourInS)

hist9D = cdfIn.variables['HIST9D'][:]
precBinned = cdfIn.variables['PRECBINNED'][:]
wpupBinned = cdfIn.variables['WPUPBINNED'][:]

lonBin = 8-1
latBin = 7-1

var1Dim = 4
var2Dim = 6
var1Bin = 8
var2Bin = 10
var1Binned = precBinned
var2Binned = wpupBinned

sumVars = range(0,8,1)
sumVars.remove(var1Dim)
sumVars.remove(var2Dim)

rightVar1Bin = numpy.where(var1Binned == var1Bin,0,-1)
rightVar2Bin = numpy.where(var2Binned == var2Bin,0,-1)

results = numpy.where(rightVar1Bin+rightVar2Bin == 0)

for ii in range(0,len(results[0])):
    currentResultTime = startDate+datetime.timedelta(0,(hourRes*3600.)*results[0][ii]
    print (str("%3i" % ii)+' '+str(currentResultTime.year)+'-'+str(currentResultTime.month)+
           '-'+str(currentResultTime.day)+' '+str(currentResultTime.hour)+':'+
           str(currentResultTime.minute)+' '+
