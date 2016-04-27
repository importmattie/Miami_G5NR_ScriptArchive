import numpy
import netCDF4
import sys
from matplotlib import pylab
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages

Nx = 90
Ny = 45
hourAv = 3
totalHoursInValues = 18264/hourAv
valueTimeOffsetBegin = 16*(24/hourAv)
valueTimeOffsetEnd = totalHoursInValues-(15*(24/hourAv))

#Nx = int(sys.argv[1])
#Ny = int(sys.argv[2])
resString = 'r'+str(Nx)+'x'+str(Ny)
#hourAv = int(sys.argv[3])

cdfInHist = netCDF4.Dataset('/home/niznik/niznik-data2/histOutput_'+resString+'_'+str(hourAv)+'_V2.nc4','r')

hist = cdfInHist.variables['HIST'][:]

monthEdges = numpy.arange(0.5,12.5+1.,1.)

lonEdges = numpy.arange(0.,360.+1.,360./10.)
latEdges = numpy.arange(-90.,90.+1.,180./9.)

precEdges = numpy.append(cdfInHist.variables['PRECBINMINS'][0],cdfInHist.variables['PRECBINMAXS'][:])
w500Edges = numpy.append(cdfInHist.variables['W500BINMINS'][0],cdfInHist.variables['W500BINMAXS'][:])
wPuPEdges = numpy.append(cdfInHist.variables['WPUPBINMINS'][0],cdfInHist.variables['WPUPBINMAXS'][:])
TEEFEdges = numpy.append(cdfInHist.variables['TEFBINMINS'][0],cdfInHist.variables['TEFBINMAXS'][:])
KEDotEdges = numpy.append(cdfInHist.variables['KEDOTBINMINS'][0],cdfInHist.variables['KEDOTBINMAXS'][:])
HMVEdges = numpy.append(cdfInHist.variables['HMVBINMINS'][0],cdfInHist.variables['HMVBINMAXS'][:])

cdfInVars = netCDF4.Dataset('/home/niznik/niznik-data2/allVars_'+resString+'_'+str(hourAv)+'.nc4','r')

precValues = cdfInVars.variables['PREC'][:]*86400.
TEEFValues = cdfInVars.variables['TEF'][:]

dimEdges = [monthEdges,lonEdges,latEdges,precEdges,w500Edges,wPuPEdges,TEEFEdges,KEDotEdges,HMVEdges]
dimLongNames = ['Month','Longitude (degrees East)','Latitude (degrees North)',
                'Precipitation (mm/day)','Vertical Velocity (m/s)',
                'Eddy Zonal Momentum Flux (m2 s-2)','Total Eddy Enthalpy Flux (J m kg-1 s-1)',
                'Integrated Wind Shear Metric (m3 s-3)','Horizonal Mesoscale Variability (m2 s-2)']
dimShortNames = ['month','lon','lat','prec','w500','wPuP','TEEF','KEDot','HMV']

def main():

    with PdfPages('/home/niznik/niznik-data2/PLOTS/LatLonScatter_1-2D_'+resString+'_'+
                  str(hourAv)+'_pr_v_TEEF.pdf') as pdfSave:
        for lat in range(5,6):
            for lon in range(0,1):
                v1 = 3
                v2 = 6
                lonAVMin = lon*(Nx/10)
                lonAVMax = lonAVMin+(Nx/10)
                latAVMin = lat*(Ny/9)
                latAVMax = latAVMin+(Ny/9)
                 hist2D = numpy.sum(hist[:,lon,lat,:,:,:,:,:,:],(0,2,3,5,6))
                plot1D_2DScatter(precValues[valueTimeOffsetBegin:valueTimeOffsetEnd,
                                            latAVMin:latAVMax,lonAVMin:lonAVMax],
                                 TEEFValues[valueTimeOffsetBegin:valueTimeOffsetEnd,
                                            latAVMin:latAVMax,lonAVMin:lonAVMax],
                                 hist2D,v1,v2,lon,lon+1,lat,lat+1,pdfSave)

def plot1D_2DScatter(values1,values2,hist,var1Pos,var2Pos,
                     lonBinMin,lonBinMax,latBinMin,latBinMax,
                     pdfFile):

    fig,axarr = plt.subplots(1,1,figsize=(7,6),dpi=300)
    fig.tight_layout(h_pad=1.0,w_pad=1.0,rect=(0.22,0.22,0.99,0.99))

    xBinEdges = dimEdges[var1Pos]
    yBinEdges = dimEdges[var2Pos]
    xBinEdgeLen = xBinEdges.size
    yBinEdgeLen = yBinEdges.size
    unlabeledX = numpy.arange(0,xBinEdgeLen,1)
    unlabeledY = numpy.arange(0,yBinEdgeLen,1)

    minPower = -7
    maxPower = 0
    levTicks = numpy.arange(minPower,maxPower+1,1)

    histLog = numpy.log10(1.0*hist/(numpy.sum(hist)))
    histLog = numpy.where(histLog == -1.*numpy.inf,minPower,histLog)

    hist1D_0 = numpy.sum(hist,1)
    hist1D_1 = numpy.sum(hist,0)

    histX = []
    histY = []

    if(var1Pos < var2Pos):
        histX = hist1D_0
        histY = hist1D_1
    else:
        histX = hist1D_1
        histY = hist1D_0

    histXLog = numpy.log10(1.0*histX/(numpy.sum(histX)))
    histXLog = numpy.where(histXLog == -1*numpy.inf,minPower,histXLog)
    histYLog = numpy.log10(1.0*histY/(numpy.sum(histY)))
    histYLog = numpy.where(histYLog == -1*numpy.inf,minPower,histYLog)

    p = None
    if(var1Pos > var2Pos):
        pointValuesX = values1.flatten()
        pointValuesY = values2.flatten()
    else:
        pointValuesX = values1.flatten()
        pointValuesY = values2.flatten()

    xBinEdgeScat,yBinEdgeScat = numpy.zeros(xBinEdgeLen),numpy.zeros(yBinEdgeLen)

    xBinEdgeScat[:] = xBinEdges[:]
    xBinEdgeScatInt = xBinEdgeScat[2]-xBinEdgeScat[1]
    xBinEdgeScat[0] = xBinEdgeScat[1]-xBinEdgeScatInt
    xBinEdgeScat[-1] = xBinEdgeScat[-2]+xBinEdgeScatInt

    yBinEdgeScat[:] = yBinEdges[:]
    yBinEdgeScatInt = yBinEdgeScat[2]-yBinEdgeScat[1]
    yBinEdgeScat[0] = yBinEdgeScat[1]-yBinEdgeScatInt
    yBinEdgeScat[-1] = yBinEdgeScat[-2]+yBinEdgeScatInt

    plottedHist = numpy.zeros(histLog.shape)

    p = axarr.pcolor(unlabeledX,unlabeledY,histLog,vmin=minPower,vmax=maxPower,cmap='Spectral_r')
    axarr.cla()

    for aa in range(0,len(pointValuesX)):
    #for aa in range(0,100):
        locPointValueX = pointValuesX[aa]
        locPointValueY = pointValuesY[aa]
        xAxisBin = numpy.searchsorted(xBinEdges,locPointValueX)-1
        yAxisBin = numpy.searchsorted(yBinEdges,locPointValueY)-1
        if(plottedHist[xAxisBin,yAxisBin] < 100):
            plottedHist[xAxisBin,yAxisBin] += 1
            logCount = 7
            if(0 <= xAxisBin <= 11 and 0 <= yAxisBin <= 11):
                logCount = histLog[xAxisBin,yAxisBin]

            if(xAxisBin == 0 and locPointValueX < xBinEdgeScat[0]):
                locPointValueX = xBinEdgeScat[0]
            if(xAxisBin == xBinEdgeLen-2 and locPointValueX > xBinEdgeScat[-1]):
                locPointValueX = xBinEdgeScat[-1]
            if(yAxisBin == 0 and locPointValueY < yBinEdgeScat[0]):
                locPointValueY = yBinEdgeScat[0]
            if(yAxisBin == yBinEdgeLen-2 and locPointValueY > yBinEdgeScat[-1]):
                locPointValueY = yBinEdgeScat[-1]

            cbPercent = (logCount-minPower)*(1.0/abs(minPower))
            pointColor = pylab.cm.Spectral_r(cbPercent)
            axarr.scatter(locPointValueX,locPointValueY,c=pointColor)

    #axarr.set_xticks(numpy.arange(0,xBinEdgeLen))
    axarr.set_xlim(xBinEdgeScat[0],xBinEdgeScat[-1])
    axarr.set_xticks(xBinEdgeScat)
    axarr.set_xticklabels(xBinEdges[0:xBinEdgeLen],rotation=270)
    axarr.tick_params(axis='x',labelsize=8)
    axarr.set_xlabel(dimLongNames[var1Pos])

    #axarr.set_yticks(numpy.arange(0,yBinEdgeLen))
    axarr.set_ylim(yBinEdgeScat[0],yBinEdgeScat[-1])
    axarr.set_yticks(yBinEdgeScat)
    axarr.set_yticklabels(yBinEdges[0:yBinEdgeLen])
    axarr.tick_params(axis='y',labelsize=8)
    axarr.set_ylabel(dimLongNames[var2Pos],rotation=90)

    for yy in range(0,len(yBinEdgeScat)):
        axarr.axhline(y=yBinEdgeScat[yy],color='#000000')
    for xx in range(0,len(xBinEdgeScat)):
        axarr.axvline(x=xBinEdgeScat[xx],color='#000000')

    axarr.set_title(str("%3i"%lonEdges[lonBinMin])+' - '+str("%3i"%lonEdges[lonBinMax])+' E | '+
                    str("%3i"%latEdges[latBinMin])+' - '+str("%3i"%latEdges[latBinMax])+' N',fontsize=8)

    cbar = plt.colorbar(p,ticks=numpy.arange(minPower,maxPower+1,1))

    #XHist                                                                                                           
    xAxHist = fig.add_axes([0.28,0.05,0.54,0.1])
    xAxHist.bar(unlabeledX[:-1],histXLog-minPower,width=1.0)
    barColorsX = []
    for ii in range(0,len(histXLog)):
        cbPercent = (histXLog[ii]-minPower)*(1.0/abs(minPower))
        barColorsX.append(pylab.cm.Spectral_r(cbPercent))
        xAxHist.bar(unlabeledX[ii],histXLog[ii]-minPower,width=1.0,color=barColorsX[ii])
    xAxHist.xaxis.set_visible(False)
    xAxHist.yaxis.set_visible(False)
    #YHist                                                                                                           
    yAxHist = fig.add_axes([0.05,0.28,0.1,0.67])
    yAxHist.barh(unlabeledY[:-1],histYLog-minPower,height=1.0)
    barColorsY = []
    for ii in range(0,len(histYLog)):
        cbPercent = (histYLog[ii]-minPower)*(1.0/abs(minPower))
        barColorsY.append(pylab.cm.Spectral_r(cbPercent))
        yAxHist.barh(unlabeledY[ii],histYLog[ii]-minPower,height=1.0,color=barColorsY[ii])
    yAxHist.xaxis.set_visible(False)
    yAxHist.yaxis.set_visible(False)

    pdfFile.savefig()

    plt.close()

if (__name__ == "__main__"):
    main()
