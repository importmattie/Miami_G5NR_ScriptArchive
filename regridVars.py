from subprocess import call
import sys

res = str(sys.argv[1])

call(['cdo','gencon,'+res,'/home/niznik/niznik-data2/calcVars.nc4','/home/niznik/niznik-data2/weights_'+res+'.nc4'])
call(['cdo','remap,'+res+',/home/niznik/niznik-data2/weights_'+res+'.nc4','/home/niznik/niznik-data2/calcVars.nc4','/home/niznik/niznik-data2/allVars_'+res+'.nc4'])
