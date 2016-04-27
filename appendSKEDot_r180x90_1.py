from subprocess import call
import sys

call('cp /home/niznik/niznik-data2/shared_repository/allVars_r180x90_1.nc4 /home/niznik/niznik-data2/shared_repository/allVars_r180x90_1_old.nc4',shell=True)
call('ncks -A -v SKEDOT_ZON /home/niznik/niznik-data2/calcVars_3D/SKEDot_r180x90.nc4 /home/niznik/niznik-data2/shared_repository/allVars_r180x90_1.nc4',shell=True)
call('ncatted -a long_name,SKEDOT_ZON,o,c,\'eddy vertical momentum flux\' /home/niznik/niznik-data2/shared_repository/allVars_r180x90_1.nc4',shell=True)
call('ncatted -a units,SKEDOT_ZON,o,c,\'m2 s-3\' /home/niznik/niznik-data2/shared_repository/allVars_r180x90_1.nc4',shell=True)
