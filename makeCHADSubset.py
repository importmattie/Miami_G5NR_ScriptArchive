import netCDF4
import numpy as np
import sys

#--- Edit These ---#
var1EdgeName = 'HMVBIN'
var2EdgeName = 'KEDOTBIN'
var1BinnedName = 'HMVBINNED'
var2BinnedName = 'KEDOTBINNED'
var1ValueName = 'HMV'
var2ValueName = 'KEDOT'
var1LongName = ''
var2LongName = ''
var1Units = ''
var2Units = ''
lonBin = 0
latBin = 5
Nx = 90
Ny = 45
hourInterval = 3
#--- End Edit

cellHourTag = 'r'+str(Nx)+'x'+str(Ny)+'_'+str(hourInterval)
cdfInHist = netCDF4.Dataset('/home/niznik/niznik-data2/histOutput_'+cellHourTag+'_V2.nc4','r')
cdfInValues = netCDF4.Dataset('/home/niznik/niznik-data2/allVars_'+cellHourTag+'.nc4','r')

varNamePos = ['HOUR','PREC','W500','WPUP','TEEF','KEDOT','HMV']
varLongNamePos = ['hour','precipitation','w at 500 hPa','w\'u\' at 500 hPa','TEEF','KEDot','HMV']
varUnits = ['hr','mm/day','m s-1','m2 s-2','J m kg-1 s-1','m3 s-3','m2 s-2']
var1Pos = varNamePos.index(var1ValueName)
var2Pos = varNamePos.index(var2ValueName)
print var1Pos,var2Pos
dimsToSumTemp = [0,1,2,3,4,5,6]
dimsToSumTemp.remove(var1Pos)
dimsToSumTemp.remove(var2Pos)
dimsToSum = tuple(dimsToSumTemp)
var1LongName = varLongNamePos[var1Pos]
var2LongName = varLongNamePos[var2Pos]
var1Units = varUnits[var1Pos]
var2Units = varUnits[var2Pos]

histName = 'HIST'
lonValueName = 'lon'
latValueName = 'lat'
timeValueName = 'time'
numLonBins = 10
numLatBins = 9
firstTimeValues = 8*16
lastTimeValues = -1*8*15

hist_2d = np.sum(cdfInHist.variables[histName][:,lonBin,latBin,:,:,:,:,:,:],dimsToSum)
var1Edges = np.append(cdfInHist.variables[var1EdgeName+'MINS'][0],cdfInHist.variables[var1EdgeName+'MAXS'][:])
var2Edges = np.append(cdfInHist.variables[var2EdgeName+'MINS'][0],cdfInHist.variables[var2EdgeName+'MAXS'][:])

lonAVMin = lonBin*(Nx/numLonBins)
lonAVMax = lonAVMin+(Nx/numLonBins)
latAVMin = latBin*(Ny/numLatBins)
latAVMax = latAVMin+(Ny/numLatBins)

var1Binned = cdfInHist.variables[var1BinnedName][:,latAVMin:latAVMax,lonAVMin:lonAVMax]
var2Binned = cdfInHist.variables[var2BinnedName][:,latAVMin:latAVMax,lonAVMin:lonAVMax]
cdfInHist.close()

var1Values = cdfInValues.variables[var1ValueName][firstTimeValues:lastTimeValues,latAVMin:latAVMax,lonAVMin:lonAVMax]
var2Values = cdfInValues.variables[var2ValueName][firstTimeValues:lastTimeValues,latAVMin:latAVMax,lonAVMin:lonAVMax]

lonValues = cdfInValues.variables[lonValueName][lonAVMin:lonAVMax]
latValues = cdfInValues.variables[latValueName][latAVMin:latAVMax]
timeValues = cdfInValues.variables[timeValueName][firstTimeValues:lastTimeValues]
cdfInValues.close()

cdfOut = netCDF4.Dataset('/home/niznik/niznik-data2/CHAD_Data/HV_'+cellHourTag+
                         '_0E5N_'+str(var1ValueName)+'_'+str(var2ValueName)+'.nc4','w')

timeDim = cdfOut.createDimension('time',len(timeValues))
latDim = cdfOut.createDimension('lat',len(latValues))
lonDim = cdfOut.createDimension('lon',len(lonValues))

var1BinDim = cdfOut.createDimension(var1EdgeName+'S',hist_2d.shape[0])
var2BinDim = cdfOut.createDimension(var2EdgeName+'S',hist_2d.shape[1])
var1BinEdgesDim = cdfOut.createDimension(var1EdgeName+'EDGES',len(var1Edges))
var2BinEdgesDim = cdfOut.createDimension(var2EdgeName+'EDGES',len(var2Edges))

timeValuesVar = cdfOut.createVariable('time','f4',('time',))
latValuesVar = cdfOut.createVariable('lat','f4',('lat',))
lonValuesVar = cdfOut.createVariable('lon','f4',('lon',))

histVar = cdfOut.createVariable('HIST','f4',(var1EdgeName+'S',var2EdgeName+'S',))
var1BinEdgesVar = cdfOut.createVariable(var1EdgeName+'EDGES','f4',(var1EdgeName+'EDGES',))
var2BinEdgesVar = cdfOut.createVariable(var2EdgeName+'EDGES','f4',(var2EdgeName+'EDGES',))
var1BinnedVar = cdfOut.createVariable(var1EdgeName+'NED','f4',('time','lat','lon',))
var2BinnedVar = cdfOut.createVariable(var2EdgeName+'NED','f4',('time','lat','lon',))

var1ValuesVar = cdfOut.createVariable(var1ValueName,'f4',('time','lat','lon',))
var2ValuesVar = cdfOut.createVariable(var2ValueName,'f4',('time','lat','lon',))

timeValuesVar[:] = timeValues[:]
latValuesVar[:] = latValues[:]
lonValuesVar[:] = lonValues[:]

histVar[:] = hist_2d[:]
var1BinEdgesVar[:] = var1Edges[:]
var2BinEdgesVar[:] = var2Edges[:]
var1BinnedVar[:] = var1Binned[:]
var2BinnedVar[:] = var2Binned[:]

var1ValuesVar[:] = var1Values[:]
var2ValuesVar[:] = var2Values[:]

setattr(var1ValuesVar,'long_name',var1LongName)
setattr(var1ValuesVar,'units',var1Units)
setattr(var2ValuesVar,'long_name',var2LongName)
setattr(var2ValuesVar,'units',var2Units)
cdfOut.close()
