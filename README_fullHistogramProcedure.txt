Here is the order of tasks that must be completed in order to successfully get to the current histogram:

1) Starting with the dlVars_r720x360_1.nc4 file, calculate the appropriate metrics and save the file to 'allVars_r720x360_1.nc4' using the script 'calcAuxVariables.py'
2) Using CDO, regrid 'allVars_r720x360_1.nc4' to the 2 degree grid (filename 'allVars_r180x90_1.nc4')
3) Append KEDot to 'allVars_r180x90_1.nc4' using the script 'appendKEDot_r180x90_1.py'
4) Using CDO, regrid 'allVars_r180x90_1.nc4' to the 4 degree grid (filename 'allVars_r90x45_1.nc4')
5) Append HMV to 'allVars_r180x90_1.nc4' and 'allVars_r90x45_1.nc4' using the script 'appendHMV.py' (needs resolution argument)
6) Regrid 'allVars_r*x??_1.nc4' to 3-hourly using the script 'regridTime_3.py' (need to run it once for each res.)
7) Run 'makeHistogram.py' for each of the spatial and temporal resolutions

Done!
