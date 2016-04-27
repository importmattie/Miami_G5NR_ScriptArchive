import netCDF4
import numpy
import calendar
import datetime
import sys

Nx = int(sys.argv[1])
Ny = int(sys.argv[2])
hourAv = int(sys.argv[3])

inputTag = 'r'+str(Nx)+'x'+str(Ny)+'_'+str(hourAv)

inputFile = '/home/niznik/niznik-data2/allVars_'+inputTag+'.nc4'
cdfIn = netCDF4.Dataset(inputFile,'r')
outputFile = '/home/niznik/niznik-data2/histOutput_'+inputTag+'_V2.nc4'

inFileFirstDay = datetime.date(2005,05,16)
inFileLastDay = datetime.date(2007,06,15)
outFileFirstDay = datetime.date(2005,06,01)
outFileLastDay = datetime.date(2007,05,31)

totalDays = (outFileLastDay-outFileFirstDay).days+1
totalHours = totalDays*24
totalHoursAv = totalHours/hourAv
hrPerDay = 24/hourAv

firstHourZ = datetime.datetime(2000,1,1)+datetime.timedelta(0,3600*(0.5+(hourAv-1)/2.))
startTimeStringZ = (str("%02i" % firstHourZ.hour)+':'+
                    str("%02i" % firstHourZ.minute)+':'+
                    str("%02i" % firstHourZ.second))

numOfLonBins = 10
numOfLatBins = 9
numOfValueBins = 12
numOfMonthBins = 12

spaceBinIntervalX = Nx/numOfLonBins
spaceBinIntervalY = Ny/numOfLatBins

#___Ranges___
#
smallInc = 1e-8
#Precip: 0- mm/day
precMin = -9.
precMax = 111.
precBEs = numpy.arange(precMin,precMax+smallInc,(precMax-precMin)/numOfValueBins)
precBEs[0] = 0.
precBEs[numOfValueBins] = 10000.
#w500: 0-2 m/s (?)
w500Min = -0.165
w500Max = 0.195
w500BEs = numpy.arange(w500Min,w500Max+smallInc,(w500Max-w500Min)/numOfValueBins)
w500BEs[0] = -10.
w500BEs[numOfValueBins] = 10.
#w'u': (?)
wPuPMin = -0.22
wPuPMax = 0.26
wPuPBEs = numpy.arange(wPuPMin,wPuPMax+smallInc,(wPuPMax-wPuPMin)/numOfValueBins)
wPuPBEs[0] = -10.
wPuPBEs[numOfValueBins] = 10.
#Cp*[w'T']+L*[w'q'] (Total Enthalpy Flux or TEF)
TEFMin = -20.
TEFMax = 460.
TEFBEs = numpy.arange(TEFMin,TEFMax+smallInc,(TEFMax-TEFMin)/numOfValueBins)
TEFBEs[0] = -1000.
TEFBEs[numOfValueBins] = 10000.
#KEDot
KEDotMin = -520.
KEDotMax = 440.
KEDotBEs = numpy.arange(KEDotMin,KEDotMax+smallInc,(KEDotMax-KEDotMin)/numOfValueBins)
KEDotBEs[0] = -10000.
KEDotBEs[numOfValueBins] = 10000.
#HMV
HMVMin = 0.
HMVMax = 48.
HMVBEs = numpy.arange(HMVMin,HMVMax+smallInc,(HMVMax-HMVMin)/numOfValueBins)
HMVBEs[numOfValueBins] = 10000.

def main():

    cdfOut = netCDF4.Dataset(outputFile,'w',format='NETCDF4')

    #lon/lat bins are static, so make them easy to find
    lonBins = numpy.zeros(Nx)
    for xx in range(0,Nx):
        lonBins[xx] = xx/spaceBinIntervalX
    latBins = numpy.zeros(Ny)
    for yy in range(0,Ny):
        latBins[yy] = yy/spaceBinIntervalY

    #Dimensions are:
    #0 - Month (static)
    #1 - Longitude (static)
    #2 - Latitude (static)
    #3 - Rain Rate (search)
    #4 - w500 (search)
    #5 - w'u' (search)
    #6 - Cp*[w'T']+L*[w'q'] (search)
    #7 - KEDot (search)
    #8 - HMV (search)

    histogram = numpy.zeros((numOfMonthBins,numOfLonBins,numOfLatBins,numOfValueBins,numOfValueBins,numOfValueBins,numOfValueBins,numOfValueBins,numOfValueBins),dtype='i4')
    precBinTracker = numpy.zeros((totalHoursAv,Ny,Nx),dtype='i1')-1
    w500BinTracker = numpy.zeros((totalHoursAv,Ny,Nx),dtype='i1')-1
    wPuPBinTracker = numpy.zeros((totalHoursAv,Ny,Nx),dtype='i1')-1
    TEFBinTracker = numpy.zeros((totalHoursAv,Ny,Nx),dtype='i1')-1
    KEDotBinTracker = numpy.zeros((totalHoursAv,Ny,Nx),dtype='i1')-1
    HMVBinTracker = numpy.zeros((totalHoursAv,Ny,Nx),dtype='i1')-1

    thisTimeStart = datetime.datetime.now()
    lastTimeStart = thisTimeStart

#For Testing...
#    for year_int in range(2005,2006):
#        for month_int in range(6,7):
#            for day_int in range(1,2):
    for year_int in range(2005,2007+1):
        for month_int in range(1,12+1):
           for day_int in range(1,calendar.monthrange(year_int,month_int)[1]+1):
                currentDay = datetime.date(year_int,month_int,day_int)
                if(currentDay <= outFileLastDay and currentDay >= outFileFirstDay):

                    daysSinceInFirst = (currentDay-inFileFirstDay).days
                    daysSinceOutFirst = (currentDay-outFileFirstDay).days

                    loadChunkSize = hrPerDay

                    precipRate = cdfIn.variables['PREC'][(daysSinceInFirst*hrPerDay):(daysSinceInFirst*hrPerDay+hrPerDay),:,:]*86400
                    w500 = cdfIn.variables['W'][(daysSinceInFirst*hrPerDay):(daysSinceInFirst*hrPerDay+hrPerDay),:,:]
                    wPuP = cdfIn.variables['WPUP'][(daysSinceInFirst*hrPerDay):(daysSinceInFirst*hrPerDay+hrPerDay),:,:]
                    TEF = cdfIn.variables['TEF'][(daysSinceInFirst*hrPerDay):(daysSinceInFirst*hrPerDay+hrPerDay),:,:]
                    KEDot = cdfIn.variables['KEDOT'][(daysSinceInFirst*hrPerDay):(daysSinceInFirst*hrPerDay+hrPerDay),:,:]
                    HMV = cdfIn.variables['HMV'][(daysSinceInFirst*hrPerDay):(daysSinceInFirst*hrPerDay+hrPerDay),:,:]

                    for hour_int in range(0,hrPerDay):
                        for xx in range(0,Nx):
                            for yy in range(0,Ny):
                                curPrecBin = findBin(precipRate[hour_int,yy,xx],precBEs)
                                curw500Bin = findBin(w500[hour_int,yy,xx],w500BEs)
                                curwPuPBin = findBin(wPuP[hour_int,yy,xx],wPuPBEs)
                                curTEFBin = findBin(TEF[hour_int,yy,xx],TEFBEs)
                                curKEDotBin = findBin(KEDot[hour_int,yy,xx],KEDotBEs)
                                curHMVBin = findBin(HMV[hour_int,yy,xx],HMVBEs)
                                if(curPrecBin > -1 and curw500Bin > -1 and curwPuPBin > -1 and curTEFBin > -1 and curKEDotBin > -1):
                                    histogram[month_int-1,lonBins[xx],latBins[yy],curPrecBin,
                                                curw500Bin,curwPuPBin,curTEFBin,curKEDotBin,curHMVBin] += 1
                                    precBinTracker[daysSinceOutFirst*(hrPerDay)+hour_int,yy,xx] = curPrecBin
                                    w500BinTracker[daysSinceOutFirst*(hrPerDay)+hour_int,yy,xx] = curw500Bin
                                    wPuPBinTracker[daysSinceOutFirst*(hrPerDay)+hour_int,yy,xx] = curwPuPBin
                                    TEFBinTracker[daysSinceOutFirst*(hrPerDay)+hour_int,yy,xx] = curTEFBin
                                    KEDotBinTracker[daysSinceOutFirst*(hrPerDay)+hour_int,yy,xx] = curKEDotBin
                                    HMVBinTracker[daysSinceOutFirst*(hrPerDay)+hour_int,yy,xx] = curHMVBin

                    lastTimeStart = thisTimeStart
                    thisTimeStart = datetime.datetime.now()
                    print (str(currentDay.year)+'-'+str(currentDay.month)+'-'+str(currentDay.day)+' ('+
                           str((thisTimeStart-lastTimeStart).total_seconds())+')')

    monthBinDim = cdfOut.createDimension('MONTHBINS',numOfMonthBins)
    lonBinDim = cdfOut.createDimension('LONBINS',numOfLonBins)
    latBinDim = cdfOut.createDimension('LATBINS',numOfLatBins)
    precBinDim = cdfOut.createDimension('PRECBINS',numOfValueBins)
    w500BinDim = cdfOut.createDimension('W500BINS',numOfValueBins)
    wPuPBinDim = cdfOut.createDimension('WPUPBINS',numOfValueBins)
    TEFBinDim = cdfOut.createDimension('TEFBINS',numOfValueBins)
    KEDotBinDim = cdfOut.createDimension('KEDOTBINS',numOfValueBins)
    HMVBinDim = cdfOut.createDimension('HMVBINS',numOfValueBins)
    valueBinDim = cdfOut.createDimension('VALUEBINS',numOfValueBins)

    lonDim = cdfOut.createDimension('LON',Nx)
    latDim = cdfOut.createDimension('LAT',Ny)
    timeDim = cdfOut.createDimension('TIME',totalHoursAv)

    monthBinVar = cdfOut.createVariable('MONTHBINS','i1',('MONTHBINS',))
    lonBinVar = cdfOut.createVariable('LONBINS','i1',('LONBINS',))
    latBinVar = cdfOut.createVariable('LATBINS','i1',('LATBINS',))
    precBinVar = cdfOut.createVariable('PRECBINS','i1',('PRECBINS',))
    w500BinVar = cdfOut.createVariable('W500BINS','i1',('W500BINS',))
    wPuPBinVar = cdfOut.createVariable('WPUPBINS','i1',('WPUPBINS',))
    TEFBinVar = cdfOut.createVariable('TEFBINS','i1',('TEFBINS',))
    KEDotBinVar = cdfOut.createVariable('KEDOTBINS','i1',('KEDOTBINS',))
    HMVBinVar = cdfOut.createVariable('HMVBINS','i1',('HMVBINS',))

    lonVar = cdfOut.createVariable('LON','f4',('LON',))
    latVar = cdfOut.createVariable('LAT','f4',('LAT',))
    timeVar = cdfOut.createVariable('TIME','i4',('TIME',))

    precBinMinsVar = cdfOut.createVariable('PRECBINMINS','f4',('VALUEBINS',))
    precBinMaxsVar = cdfOut.createVariable('PRECBINMAXS','f4',('VALUEBINS',))
    w500BinMinsVar = cdfOut.createVariable('W500BINMINS','f4',('VALUEBINS',))
    w500BinMaxsVar = cdfOut.createVariable('W500BINMAXS','f4',('VALUEBINS',))
    wPuPBinMinsVar = cdfOut.createVariable('WPUPBINMINS','f4',('VALUEBINS',))
    wPuPBinMaxsVar = cdfOut.createVariable('WPUPBINMAXS','f4',('VALUEBINS',))
    TEFBinMinsVar = cdfOut.createVariable('TEFBINMINS','f4',('VALUEBINS',))
    TEFBinMaxsVar = cdfOut.createVariable('TEFBINMAXS','f4',('VALUEBINS',))
    KEDotBinMinsVar = cdfOut.createVariable('KEDOTBINMINS','f4',('VALUEBINS',))
    KEDotBinMaxsVar = cdfOut.createVariable('KEDOTBINMAXS','f4',('VALUEBINS',))
    HMVBinMinsVar = cdfOut.createVariable('HMVBINMINS','f4',('VALUEBINS',))
    HMVBinMaxsVar = cdfOut.createVariable('HMVBINMAXS','f4',('VALUEBINS',))

    histogramVar = cdfOut.createVariable('HIST','i4',('MONTHBINS','LONBINS','LATBINS',
                                                      'PRECBINS','W500BINS','WPUPBINS',
                                                      'TEFBINS','KEDOTBINS','HMVBINS'))
    precBinTrackerVar = cdfOut.createVariable('PRECBINNED','i1',('TIME','LAT','LON',))
    w500BinTrackerVar = cdfOut.createVariable('W500BINNED','i1',('TIME','LAT','LON',))
    wPuPBinTrackerVar = cdfOut.createVariable('WPUPBINNED','i1',('TIME','LAT','LON',))
    TEFBinTrackerVar = cdfOut.createVariable('TEFBINNED','i1',('TIME','LAT','LON',))
    KEDotBinTrackerVar = cdfOut.createVariable('KEDOTBINNED','i1',('TIME','LAT','LON',))
    HMVBinTrackerVar = cdfOut.createVariable('HMVBINNED','i1',('TIME','LAT','LON',))

    monthBinVar[:] = numpy.arange(1,numOfMonthBins+1,1)
    lonBinVar[:] = numpy.arange(1,numOfLonBins+1,1)
    latBinVar[:] = numpy.arange(1,numOfLatBins+1,1)
    precBinVar[:] = numpy.arange(1,numOfValueBins+1,1)
    w500BinVar[:] = numpy.arange(1,numOfValueBins+1,1)
    wPuPBinVar[:] = numpy.arange(1,numOfValueBins+1,1)
    TEFBinVar[:] = numpy.arange(1,numOfValueBins+1,1)
    KEDotBinVar[:] = numpy.arange(1,numOfValueBins+1,1)
    HMVBinVar[:] = numpy.arange(1,numOfValueBins+1,1)

    precBinMinsVar[:] = precBEs[0:numOfValueBins]
    precBinMaxsVar[:] = precBEs[1:numOfValueBins+1]
    w500BinMinsVar[:] = w500BEs[0:numOfValueBins]
    w500BinMaxsVar[:] = w500BEs[1:numOfValueBins+1]
    wPuPBinMinsVar[:] = wPuPBEs[0:numOfValueBins]
    wPuPBinMaxsVar[:] = wPuPBEs[1:numOfValueBins+1]
    TEFBinMinsVar[:] = TEFBEs[0:numOfValueBins]
    TEFBinMaxsVar[:] = TEFBEs[1:numOfValueBins+1]
    KEDotBinMinsVar[:] = KEDotBEs[0:numOfValueBins]
    KEDotBinMaxsVar[:] = KEDotBEs[1:numOfValueBins+1]
    HMVBinMinsVar[:] = HMVBEs[0:numOfValueBins]
    HMVBinMaxsVar[:] = HMVBEs[1:numOfValueBins+1]

    lonVar[:] = cdfIn.variables['lon'][:]
    latVar[:] = cdfIn.variables['lat'][:]
    timeVar[:] = numpy.arange(0,totalHours*60,hourAv*60)

    histogramVar[:] = histogram[:]
    precBinTrackerVar[:] = precBinTracker[:]
    w500BinTrackerVar[:] = w500BinTracker[:]
    wPuPBinTrackerVar[:] = wPuPBinTracker[:]
    TEFBinTrackerVar[:] = TEFBinTracker[:]
    KEDotBinTrackerVar[:] = KEDotBinTracker[:]
    HMVBinTrackerVar[:] = HMVBinTracker[:]

    setattr(lonVar,'long_name','longitude')
    setattr(lonVar,'units','degrees_east')

    setattr(latVar,'long_name','latitude')
    setattr(latVar,'units','degrees_north')

    setattr(timeVar,'long_name','time')
    setattr(timeVar,'units','minutes since 2005-06-01 '+startTimeStringZ)

    cdfIn.close()
    cdfOut.close()

def findBin(value,binEdges):
    position = numpy.searchsorted(binEdges,value)-1
    if(position >= 0 and position <= numOfValueBins-1):
        return position
    return -1

if __name__ == "__main__":
    main()
