import netCDF4
import numpy as np
import sys

inputPath = '/home/niznik/niznik-data2/dlVars_3D/'
resTag = ''

uInputFile = inputPath+'U'+resTag+'.nc4'
vInputFile = inputPath+'V'+resTag+'.nc4'
wInputFile = inputPath+'W'+resTag+'.nc4'
wuInputFile = inputPath+'WU'+resTag+'.nc4'
wvInputFile = inputPath+'WV'+resTag+'.nc4'

uCDFIn = netCDF4.Dataset(uInputFile,'r')
vCDFIn = netCDF4.Dataset(vInputFile,'r')
wCDFIn = netCDF4.Dataset(wInputFile,'r')
wuCDFIn = netCDF4.Dataset(wuInputFile,'r')
wvCDFIn = netCDF4.Dataset(wvInputFile,'r')

lons = uCDFIn.variables['lon'][:]
lats = uCDFIn.variables['lat'][:]
levs = uCDFIn.variables['lev'][:]*100.
z100 = np.where(levs == 10000.)[0][0]
z500 = np.where(levs == 50000.)[0][0]
z1000 = np.where(levs == 100000.)[0][0]
times = uCDFIn.variables['time'][:]

outputFilepath = '/home/niznik/niznik-data2/calcVars_3D/SKEDot'+resTag+'.nc4'
cdfOut = netCDF4.Dataset(outputFilepath,'w')

lonDim = cdfOut.createDimension('lon',len(lons))
latDim = cdfOut.createDimension('lat',len(lats))
levDim = cdfOut.createDimension('lev',len(levs))
timeDim = cdfOut.createDimension('time',len(times))

lonVar = cdfOut.createVariable('lon','f4',('lon',))
latVar = cdfOut.createVariable('lat','f4',('lat',))
levVar = cdfOut.createVariable('lev','f4',('lev',))
timeVar = cdfOut.createVariable('time','i4',('time',))

lonVar[:] = lons[:]
latVar[:] = lats[:]
levVar[:] = levs[:]
timeVar[:] = times[:]

SKEDot_dpVar = cdfOut.createVariable('SKEDOT_DP','f4',('time','lev','lat','lon',))
SKEDotVar = cdfOut.createVariable('SKEDOT','f4',('time','lat','lon',))
SKEDot_ZonVar = cdfOut.createVariable('SKEDOT_ZON','f4',('time','lat','lon',))
SKEDot_LAVar = cdfOut.createVariable('SKEDOT_LA','f4',('time','lat','lon',))
SKEDot_UAVar = cdfOut.createVariable('SKEDOT_UA','f4',('time','lat','lon',))

#Based on a table found at:
#usatoday30.usatoday.com/weather/wstdatmo.htm
rho = -4e-11*levs**2+2e-5*levs+0.0256

for tt in range(0,len(times)):

    print 'Started time '+str(tt)+' of 18624...'

    #Load the chunked variables
    u = uCDFIn.variables['U'][tt,:,:,:]
    v = vCDFIn.variables['V'][tt,:,:,:]
    w = wCDFIn.variables['W'][tt,:,:,:]
    wu = wuCDFIn.variables['WU'][tt,:,:,:]
    wv = wvCDFIn.variables['WV'][tt,:,:,:]

    SKEDot_dp_Zon = np.zeros(u.shape)
    SKEDot_dp_Mer = np.zeros(u.shape)
    SKEDot_dp = np.zeros(u.shape)
    for zz in range(1,len(levs)-1):
        
        #SKEDot = VertInt( -d/dp([uw]-[u][w])*u - d/dp([vw]-[v][w])*v)
        lev1 = zz+1
        lev2 = zz-1
        pDiff = levs[lev1]-levs[lev2]
        #SKEDot_dp_Zon = (-1.*((wu[lev1,:,:]-u[lev1,:,:]*w[lev1,:,:])-(wu[lev2,:,:]-u[lev2,:,:]*w[lev2,:,:]))/pDiff*u[zz,:,:])
        #SKEDot_dp_Mer = (-1.*((wv[lev1,:,:]-v[lev1,:,:]*w[lev1,:,:])-(wv[lev2,:,:]-v[lev2,:,:]*w[lev2,:,:]))/pDiff*v[zz,:,:])
        SKEDot_dp_Zon[zz,:,:] = (-1.*((wu[lev1,:,:]-u[lev1,:,:]*w[lev1,:,:])-(wu[lev2,:,:]-u[lev2,:,:]*w[lev2,:,:]))/pDiff*(-1.*u[zz,:,:]*rho[zz]*9.81))
        SKEDot_dp_Mer[zz,:,:] = (-1.*((wv[lev1,:,:]-v[lev1,:,:]*w[lev1,:,:])-(wv[lev2,:,:]-v[lev2,:,:]*w[lev2,:,:]))/pDiff*(-1.*v[zz,:,:]*rho[zz]*9.81))
        SKEDot_dp[zz,:,:] = SKEDot_dp_Zon[zz,:,:]+SKEDot_dp_Mer[zz,:,:]

    SKEDot_dp_Zon[:] = np.where(np.abs(SKEDot_dp_Zon) > 1e6,0,SKEDot_dp_Zon)
    SKEDot_dp[:] = np.where(np.abs(SKEDot_dp) > 1e6,0,SKEDot_dp)

    #calculate SKEDot/SKEDot(1000-500)/SKEDot(500-100)
    SKEDot = np.zeros(SKEDot_dp[0,:,:].shape)
    SKEDot_Zon = np.zeros(SKEDot_dp_Zon[0,:,:].shape)
    SKEDot_Lower = np.zeros(SKEDot_dp[0,:,:].shape)
    SKEDot_Upper = np.zeros(SKEDot_dp[0,:,:].shape)
    for zz in range(z100,z500):
        pWidth = (levs[zz+1]-levs[zz-1])/2.
        SKEDot[:] += SKEDot_dp[zz,:,:]*pWidth
        SKEDot_Zon[:] += SKEDot_dp_Zon[zz,:,:]*pWidth
        SKEDot_Upper[:] += SKEDot_dp[zz,:,:]*pWidth
    for zz in range(z500,z1000):
        pWidth = (levs[zz+1]-levs[zz-1])/2.
        SKEDot[:] += SKEDot_dp[zz,:,:]*pWidth
        SKEDot_Zon[:] += SKEDot_dp_Zon[zz,:,:]*pWidth
        SKEDot_Lower[:] += SKEDot_dp[zz,:,:]*pWidth

    #Make it an average instead of an integral
    SKEDot = SKEDot/((1000.*100.)-(100.*100.))
    SKEDot_Zon = SKEDot_Zon/((1000.*100.)-(100.*100.))
    SKEDot_Upper = SKEDot_Upper/((500.*100.)-(100.*100.))
    SKEDot_Lower = SKEDot_Lower/((1000.*100.)-(500.*100.))

    #save this time of SKEDot
    SKEDot_dpVar[tt,:,:,:] = SKEDot_dp[:,:,:]
    SKEDotVar[tt,:,:] = SKEDot[:,:]
    SKEDot_ZonVar[tt,:,:] = SKEDot_Zon[:,:]
    SKEDot_LAVar[tt,:,:] = SKEDot_Lower[:,:]
    SKEDot_UAVar[tt,:,:] = SKEDot_Upper[:,:]

setattr(lonVar,'long_name','longitude')
setattr(lonVar,'units','degrees_east')

setattr(latVar,'long_name','latitude')
setattr(latVar,'units','degrees_north')

setattr(levVar,'long_name','pressure')
setattr(levVar,'units','hPa')

setattr(timeVar,'long_name','time')
setattr(timeVar,'units','minutes since 2005-05-16 00:30:00')

setattr(SKEDotVar,'long_name','Integral(-d/dp([uw]-[u][w])*u - d/dp([vw]-[v][w])*v)')
setattr(SKEDotVar,'units','m2 s-3')

setattr(SKEDot_ZonVar,'long_name','Integral(-d/dp([uw]-[u][w])*u)')
setattr(SKEDot_ZonVar,'units','m2 s-3')

setattr(SKEDot_LAVar,'long_name','SKEDot, from 1000 hPa - 500 hPa')
setattr(SKEDot_LAVar,'units','m2 s-3')

setattr(SKEDot_UAVar,'long_name','SKEDot, from 500 hPa - 100 hPa')
setattr(SKEDot_UAVar,'units','m2 s-3')

setattr(SKEDot_dpVar,'long_name','dSKEDot/dp')
setattr(SKEDot_dpVar,'units','m2 s-3 Pa-1')

uCDFIn.close()
vCDFIn.close()
wCDFIn.close()
wuCDFIn.close()
wvCDFIn.close()
cdfOut.close()
