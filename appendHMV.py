import netCDF4
import numpy
from subprocess import call
import sys

res = str(sys.argv[1])

call('cp /home/niznik/niznik-data2/shared_repository/allVars_'+res+'_1.nc4 /home/niznik/niznik-data2/shared_repository/allVars_'+res+'_1.nc4.bckp',shell=True)

inputFilepath = '/home/niznik/niznik-data2/shared_repository/allVars_'+res+'_1.nc4'
cdfIn = netCDF4.Dataset(inputFilepath,'a')

U = cdfIn.variables['U'][:]
V = cdfIn.variables['V'][:]
UU = cdfIn.variables['UU'][:]
VV = cdfIn.variables['VV'][:]

HMV = (UU-(U*U))+(VV-(V*V))

HMVVar = cdfIn.createVariable('HMV','f4',('time','lat','lon',))
HMVVar[:] = HMV[:]

setattr(HMVVar,'long_name','Horizontal Mesoscale Variability')
setattr(HMVVar,'units','m2 s-2')

cdfIn.close()
