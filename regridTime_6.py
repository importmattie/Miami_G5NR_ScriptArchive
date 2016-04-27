from subprocess import call
import sys

res = str(sys.argv[1])

for tt in range(0,18264,6):
    t1 = tt
    t2 = t1+2
    print t1
    call('ncra -O -d time,'+str(t1)+','+str(t2)+
         ' /home/niznik/niznik-data2/shared_repository/allVars_'+res+'_1.nc4'+
         ' /home/niznik/niznik-data2/TempCat/'+str("%05d" % t1)+'_'+res+'.nc4',shell=True)

call('ncrcat -O /home/niznik/niznik-data2/TempCat/?????_'+res+'.nc4'+
     ' /home/niznik/niznik-data2/shared_repository/allVars_'+res+'_6.nc4',shell=True)
