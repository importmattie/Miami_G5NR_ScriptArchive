import netCDF4
import numpy

inputFilepath = '/home/niznik/niznik-data2/dlVars_r720x360_1.nc4'
cdfIn = netCDF4.Dataset(inputFilepath,'r')

lonInput = cdfIn.variables['lon'][:]
latInput = cdfIn.variables['lat'][:]
timeInput = cdfIn.variables['time'][:]

outputFilepath = '/home/niznik/niznik-data2/allVars_r720x360_1.nc4'
cdfOut = netCDF4.Dataset(outputFilepath,'w')

#set up all the info we need for the output file,
#then just fill it gradually

lonDim = cdfOut.createDimension('lon',len(lonInput))
latDim = cdfOut.createDimension('lat',len(latInput))
timeDim = cdfOut.createDimension('time',len(timeInput))

lonVar = cdfOut.createVariable('lon','f4',('lon',))
latVar = cdfOut.createVariable('lat','f4',('lat',))
timeVar = cdfOut.createVariable('time','i4',('time',))

lonVar[:] = lonInput[:]
latVar[:] = latInput[:]
timeVar[:] = timeInput[:]

precVar = cdfOut.createVariable('PREC','f4',('time','lat','lon',))
QVar = cdfOut.createVariable('Q','f4',('time','lat','lon',))
TVar = cdfOut.createVariable('T','f4',('time','lat','lon',))
UVar = cdfOut.createVariable('U','f4',('time','lat','lon',))
VVar = cdfOut.createVariable('V','f4',('time','lat','lon',))
WVar = cdfOut.createVariable('W','f4',('time','lat','lon',))

WQVar = cdfOut.createVariable('WQ','f4',('time','lat','lon',))
WTVar = cdfOut.createVariable('WT','f4',('time','lat','lon',))
WUVar = cdfOut.createVariable('WU','f4',('time','lat','lon',))
WVVar = cdfOut.createVariable('WV','f4',('time','lat','lon',))

wPtPVar = cdfOut.createVariable('WPTP','f4',('time','lat','lon',))
wPqPVar = cdfOut.createVariable('WPQP','f4',('time','lat','lon',))
wPuPVar = cdfOut.createVariable('WPUP','f4',('time','lat','lon',))
TEFVar = cdfOut.createVariable('TEF','f4',('time','lat','lon',))
UUVar = cdfOut.createVariable('UU','f4',('time','lat','lon',))
VVVar = cdfOut.createVariable('VV','f4',('time','lat','lon',))

numOfHours = 18288
numOfChunks = 48
div = numOfHours/numOfChunks

#Need to calculate:
#[w'T'],[w'u'],[w'q']
#cp*[w'T']+L*[w'q'] (J*m*kg-1*s-1)+(J*m*kg-1*s-1)

cp = 1005. #J*kg-1*K-1
L = 2.5e6  #J*kg-1

for tt in range(0,numOfChunks):
    firstHour = tt*div
    lastHour = tt*div+div-1
    print 'Processing hours '+str(firstHour)+' through '+str(lastHour)+'...'
    prec = cdfIn.variables['PRECTOT'][firstHour:lastHour+1,:,:]
    Q = cdfIn.variables['QV'][firstHour:lastHour+1,0,:,:]
    T = cdfIn.variables['T'][firstHour:lastHour+1,0,:,:]
    U = cdfIn.variables['U'][firstHour:lastHour+1,0,:,:]
    V = cdfIn.variables['V'][firstHour:lastHour+1,0,:,:]
    W = cdfIn.variables['W'][firstHour:lastHour+1,0,:,:]
    WQ = cdfIn.variables['WQ'][firstHour:lastHour+1,0,:,:]
    WT = cdfIn.variables['WT'][firstHour:lastHour+1,0,:,:]
    WU = cdfIn.variables['WU'][firstHour:lastHour+1,0,:,:]
    WV = cdfIn.variables['WV'][firstHour:lastHour+1,0,:,:]

    wPqP = WQ-(W*Q)
    wPtP = WT-(W*T)
    wPuP = WU-(W*U)
    TEF = (cp*wPtP)+(L*wPqP)
    UU = U*U
    VV = V*V

    precVar[firstHour:lastHour+1,:,:] = prec[:]
    QVar[firstHour:lastHour+1,:,:] = Q[:]
    TVar[firstHour:lastHour+1,:,:] = T[:]
    UVar[firstHour:lastHour+1,:,:] = U[:]
    VVar[firstHour:lastHour+1,:,:] = V[:]
    WVar[firstHour:lastHour+1,:,:] = W[:]

    WQVar[firstHour:lastHour+1,:,:] = WQ[:]
    WTVar[firstHour:lastHour+1,:,:] = WT[:]
    WUVar[firstHour:lastHour+1,:,:] = WU[:]
    WVVar[firstHour:lastHour+1,:,:] = WV[:]

    wPtPVar[firstHour:lastHour+1,:,:] = wPtP[:]
    wPqPVar[firstHour:lastHour+1,:,:] = wPqP[:]
    wPuPVar[firstHour:lastHour+1,:,:] = wPuP[:]
    TEFVar[firstHour:lastHour+1,:,:] = TEF[:]
    UUVar[firstHour:lastHour+1,:,:] = UU[:]
    VVVar[firstHour:lastHour+1,:,:] = VV[:]

setattr(lonVar,'long_name','longitude')
setattr(lonVar,'units','degrees_east')

setattr(latVar,'long_name','latitude')
setattr(latVar,'units','degrees_north')

setattr(timeVar,'long_name','time')
setattr(timeVar,'units','minutes since 2005-05-16 00:30:00')

setattr(precVar,'long_name','total_precipitation')
setattr(precVar,'units','kg m-2 s-1')

setattr(QVar,'long_name','specific_humidity')
setattr(QVar,'units','kg kg-1')

setattr(TVar,'long_name','air_temperature')
setattr(TVar,'units','K')

setattr(UVar,'long_name','eastward_wind')
setattr(UVar,'units','m s-1')

setattr(VVar,'long_name','northward_wind')
setattr(VVar,'units','m s-1')

setattr(WVar,'long_name','vertical_velocity')
setattr(WVar,'units','m s-1')

setattr(WQVar,'long_name','mean vertical moisture flux')
setattr(WQVar,'units','kg m kg-1 s-1')

setattr(WTVar,'long_name','mean vertical heat flux')
setattr(WTVar,'units','K m s-1')

setattr(WUVar,'long_name','mean vertical zonal momentum flux')
setattr(WUVar,'units','m2 s-2')

setattr(WVVar,'long_name','mean vertical meridional momentum flux')
setattr(WVVar,'units','m2 s-2')

setattr(wPtPVar,'long_name','eddy vertical heat flux (w\'T\')')
setattr(wPtPVar,'units','K m s-1')

setattr(wPqPVar,'long_name','eddy vertical moisture flux (w\'q\')')
setattr(wPqPVar,'units','m s-1')

setattr(wPuPVar,'long_name','eddy vertical zonal momentum flux (w\'u\')')
setattr(wPuPVar,'units','m2 s-2')

setattr(TEFVar,'long_name','eddy vertical enthalpy flux (cp*w\'T\'+L*w\'q\')')
setattr(TEFVar,'units','J m kg-1 s-1')

setattr(UUVar,'long_name','square of eastward_wind')
setattr(UUVar,'units','m2 s-2')

setattr(VVVar,'long_name','square of northward_wind')
setattr(VVVar,'units','m2 s-2')

cdfIn.close()
cdfOut.close()
