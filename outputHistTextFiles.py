import os
import subprocess
from subprocess import call

outputDir = '/home/niznik/niznik-data2/HistText/prec_w_TEF/'
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

lonRes = 90
latRes = 45
hourRes = 3

resFolder = str(lonRes)+'_'+str(latRes)+'_'+str(hourRes)+'/'

if not os.path.exists(outputDir+resFolder):
    os.makedirs(outputDir+resFolder)

lenTolerance = 10
lonBin = 6
latBin = 3

binFormatString = "%02i"

for precB in range(0,12):
    for wB in range(0,12):
        for TEFB in range(0,12):
            datOutputName = str(binFormatString%precB)+'_'+str(binFormatString%wB)+'_'+str(binFormatString%TEFB)+'.dat'
            fullPathToDatOutput = outputDir+resFolder+datOutputName
            call('python /home/niznik/bin/GEOS5/eventFinderAdapt.py '+str(lonRes)+' '+str(latRes)+' '+str(hourRes)+' '+
                 str(lonBin)+' '+str(latBin)+' '+str(precB)+' '+str(wB)+' -1 '+str(TEFB)+' -1 -1 1 > '+
                 fullPathToDatOutput,shell=True)
            print fullPathToDatOutput
            length = int(subprocess.check_output('wc -l '+fullPathToDatOutput,shell=True).split()[0])
            size = os.stat(fullPathToDatOutput).st_size
            if(length > lenTolerance or size <= 1):
                call('rm '+fullPathToDatOutput,shell=True)
                print datOutputName+' omitted'
            else:
                print datOutputName+' created'
                
