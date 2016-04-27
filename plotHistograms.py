import sys
import netCDF4
import numpy
from mpl_toolkits.basemap import Basemap
from matplotlib import pylab
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages

Nx = int(sys.argv[1])
Ny = int(sys.argv[2])
resString = 'r'+str(Nx)+'x'+str(Ny)
hourAv = int(sys.argv[3])

cdfIn = netCDF4.Dataset('/home/niznik/niznik-data2/histOutput_'+resString+'_'+str(hourAv)+'_V2.nc4','r')

hist = cdfIn.variables['HIST'][:]

monthEdges = numpy.arange(0.5,12.5+.1,1.)

lonEdges = numpy.arange(0.,360.+.1,360./10.)
latEdges = numpy.arange(-90.,90.+.1,180./9.)

precEdges = numpy.append(cdfIn.variables['PRECBINMINS'][0],cdfIn.variables['PRECBINMAXS'][:])
w500Edges = numpy.append(cdfIn.variables['W500BINMINS'][0],cdfIn.variables['W500BINMAXS'][:])
wPuPEdges = numpy.append(cdfIn.variables['WPUPBINMINS'][0],cdfIn.variables['WPUPBINMAXS'][:])
TEFEdges = numpy.append(cdfIn.variables['TEFBINMINS'][0],cdfIn.variables['TEFBINMAXS'][:])
KEDotEdges = numpy.append(cdfIn.variables['KEDOTBINMINS'][0],cdfIn.variables['KEDOTBINMAXS'][:])
HMVEdges = numpy.append(cdfIn.variables['HMVBINMINS'][0],cdfIn.variables['HMVBINMAXS'][:])

dimEdges = [monthEdges,lonEdges,latEdges,precEdges,w500Edges,wPuPEdges,TEFEdges,KEDotEdges,HMVEdges]
dimLongNames = ['Month','Longitude (degrees East)','Latitude (degrees North)',
                'Precipitation (mm/day)','Vertical Velocity (m/s)',
                'Eddy Zonal Momentum Flux (m2 s-2)','Total Eddy Enthalpy Flux (J m kg-1 s-1)',
                'Integrated Wind Shear Metric (m3 s-3)','Horizonal Mesoscale Variability (m2 s-2)']
dimShortNames = ['month','lon','lat','prec','w500','wPuP','TEEF','KEDot','HMV']

def main():

    #Simple Global Histograms
    #with PdfPages('/home/niznik/niznik-data2/PLOTS/sumAllDims_1D_'+resString+'_'+str(hourAv)+'.pdf') as pdfSave:
    #    for vv in range(0,9):
    #        print 'Plotting '+str(vv)+' in 1D...'
    #        dimsToSum = range(8,-1,-1)
    #        dimsToSum.remove(vv)
    #        hist1D = hist
    #        for dd in dimsToSum:
    #            hist1D = numpy.sum(hist1D,dd)
    #        plot1DHistogram(hist1D,vv,pdfSave)

    #Simple 2D Global Histograms
    #with PdfPages('/home/niznik/niznik-data2/PLOTS/sumAllDims_2D_'+resString+'_'+str(hourAv)+'.pdf') as pdfSave:
    #    for v1 in range(0,9):
    #        for v2 in range(0,9):
    #            if(v2 != v1):
    #                print 'Plotting '+str(v1)+' and '+str(v2)+' in 2D...'
    #                dimsToSum = range(8,-1,-1)
    #                dimsToSum.remove(v1)
    #                dimsToSum.remove(v2)
    #                hist2D = hist
    #                for dd in dimsToSum:
    #                    hist2D = numpy.sum(hist2D,dd)
    #                plot2DHistogram(hist2D,v1,v2,pdfSave)

    #2D Histograms in Lat/Lon Grid
    #with PdfPages('/home/niznik/niznik-data2/PLOTS/LatLonGrid_2D_'+resString+'_'+
    #              str(hourAv)+'.pdf') as pdfSave:
    #    for v1 in range(0,9):
    #        for v2 in range(0,9):
    #            if(v2 != v1 and v1 != dimShortNames.index("lon") and v1 != dimShortNames.index("lat") 
    #               and v2 != dimShortNames.index("lon") and v2 != dimShortNames.index("lat")):
    #                print 'Plotting '+str(v1)+' vs '+str(v2)
    #                dimsToSum = range(8,-1,-1)
    #                dimsToSum.remove(v1)
    #                dimsToSum.remove(v2)
    #                dimsToSum.remove(dimShortNames.index("lon"))
    #                dimsToSum.remove(dimShortNames.index("lat"))
    #                hist2D = hist
    #                for dd in dimsToSum:
    #                    hist2D = numpy.sum(hist2D,dd)
    #                plot2DHistogramLatLonMap(hist2D,v1,v2,pdfSave)

    #2D Histograms, condition: precip = 0
    #with PdfPages('/home/niznik/niznik-data2/PLOTS/LatLonGrid_2D_'+resString+'_'+
    #              str(hourAv)+'_prec0.pdf') as pdfSave:
    #    for v1 in range(0,9):
    #        for v2 in range(0,9):
    #            if(v1 != v2 and (v1 not in set([2,3,4])) and (v2 not in set([2,3,4])) ):
    #                print 'Plotting '+str(v1)+' vs '+str(v2)
    #                #make a list of dimensions for summing
    #                dimsToSum = range(8,-1,-1)
    #                dimsToSum.remove(v1)
    #                dimsToSum.remove(v2)
    #                dimsToSum.remove(2)
    #                dimsToSum.remove(3)
    #                dimsToSum.remove(4)
    #                #subset precip and shift dims accordingly
    #                hist2D = hist[:,:,:,:,0,:,:,:,:]
    #                for ii in range(0,len(dimsToSum)):
    #                    if(dimsToSum[ii] > 4):
    #                        dimsToSum[ii] -= 1
    #                #sum and plot
    #                for dd in dimsToSum:
    #                    hist2D = numpy.sum(hist2D,dd)
    #                plot2DHistogramLatLonMap(hist2D,v1,v2,pdfSave)

    #1/2D Histograms, all regions, only precip vs w500 (dim 1 and 2 when lon/lat defined)
    #with PdfPages('/home/niznik/niznik-data2/PLOTS/LatLonGrid_1-2D_'+resString+'_'+
    #              str(hourAv)+'_pr_v_w500.pdf') as pdfSave:
    #    for lat in range(0,9):
    #        for lon in range(0,10):
    #            v1 = 3
    #            v2 = 4
    #            hist2D = numpy.sum(hist[:,lon,lat,:,:,:,:,:,:],(0,3,4,5,6))
    #            plot1D_2DHistogram(hist2D,v1,v2,lon,lon+1,lat,lat+1,pdfSave)

    #1/2D Histograms, all regions, only precip vs TEEF (dim 1 and 2 when lon/lat defined)
    with PdfPages('/home/niznik/niznik-data2/PLOTS/LatLonGrid_1-2D_'+resString+'_'+
                  str(hourAv)+'_pr_v_TEEF.pdf') as pdfSave:
        for lat in range(0,9):
            for lon in range(0,10):
                v1 = 3
                v2 = 6
                hist2D = numpy.sum(hist[:,lon,lat,:,:,:,:,:,:],(0,2,3,5,6))
                plot1D_2DHistogram(hist2D,v1,v2,lon,lon+1,lat,lat+1,pdfSave)

    return

def plot1DHistogram(histToPlot,varPos,pdfFile):

    fig,axarr = plt.subplots(1,1,figsize=(6,6),dpi=300)

    binEdges = dimEdges[varPos]
    unlabeledX = numpy.arange(0,histToPlot.size,1)

    histLog = numpy.log10(1.0*histToPlot/numpy.sum(histToPlot))

    axarr.plot(unlabeledX,histLog,'-',lw=5,color='#000000')

    axarr.set_xticks(unlabeledX[0:histToPlot.size:1])
    axarr.set_xticklabels(binEdges[0:histToPlot.size:1])
    axarr.set_xlabel(dimLongNames[varPos])
    axarr.set_ylim(-6,0,1)
    axarr.set_ylabel('Count (Log Base 10)')
    axarr.tick_params(axis='both',which='major',labelsize=6)

    pdfFile.savefig()

    plt.close()

    #plt.savefig('/home/niznik/niznik-data2/PLOTS/'+resString+'_1D_'+plotTag+'_'+dimShortNames[varPos]+'.pdf')

def plot2DHistogram(histToPlot,varPos1,varPos2,pdfFile):

    fig,axarr = plt.subplots(1,1,figsize=(6,6),dpi=300)
    fig.tight_layout(h_pad=1.0,w_pad=1.0,rect=(0.1,0.1,0.98,0.98))

    xBinEdges = dimEdges[varPos1]
    yBinEdges = dimEdges[varPos2]
    xBinEdgeLen = xBinEdges.size
    yBinEdgeLen = yBinEdges.size
    unlabeledX = numpy.arange(0,xBinEdgeLen,1)
    unlabeledY = numpy.arange(0,yBinEdgeLen,1)

    minPower = -7
    maxPower = 0
    levTicks = numpy.arange(minPower,maxPower+1,1)

    histLog = numpy.log10(1.0*histToPlot/(numpy.sum(histToPlot)))
    histLog = numpy.where(histLog == -1.*numpy.inf,minPower,histLog)

    p = None
    if(varPos1 > varPos2):
        p = axarr.pcolor(unlabeledX,unlabeledY,histLog,vmin=minPower,vmax=maxPower,cmap='Spectral_r')
    else:
        p = axarr.pcolor(unlabeledX,unlabeledY,numpy.transpose(histLog),vmin=minPower,vmax=maxPower,cmap='Spectral_r')

    axarr.axhline(y=unlabeledY[1],color='#000000')
    axarr.axhline(y=unlabeledY[-2],color='#000000')
    axarr.axvline(x=unlabeledX[1],color='#000000')
    axarr.axvline(x=unlabeledX[-2],color='#000000')

    axarr.set_xticks(numpy.arange(0,xBinEdgeLen))
    axarr.set_xticklabels(xBinEdges[0:xBinEdgeLen],rotation=270)
    axarr.tick_params(axis='x',labelsize=8)
    axarr.set_xlabel(dimLongNames[varPos1])

    axarr.set_yticks(numpy.arange(0,yBinEdgeLen))
    axarr.set_yticklabels(yBinEdges[0:yBinEdgeLen])
    axarr.tick_params(axis='y',labelsize=8)
    axarr.set_ylabel(dimLongNames[varPos2],rotation=90)

    cbar = plt.colorbar(p,ticks=numpy.arange(minPower,maxPower+1,1))

    pdfFile.savefig()

    plt.close()

def plot1D_2DHistogram(hist,varPos1,varPos2,
                       lonBinMin,lonBinMax,latBinMin,latBinMax,
                       pdfFile):

    fig,axarr = plt.subplots(1,1,figsize=(7,6),dpi=300)
    fig.tight_layout(h_pad=1.0,w_pad=1.0,rect=(0.22,0.22,0.99,0.99))

    xBinEdges = dimEdges[varPos1]
    yBinEdges = dimEdges[varPos2]
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

    if(varPos1 < varPos2):
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
    if(varPos1 > varPos2):
        p = axarr.pcolor(unlabeledX,unlabeledY,histLog,vmin=minPower,vmax=maxPower,cmap='Spectral_r')
    else:
        p = axarr.pcolor(unlabeledX,unlabeledY,numpy.transpose(histLog),vmin=minPower,vmax=maxPower,cmap='Spectral_r')

    axarr.axhline(y=unlabeledY[1],color='#000000')
    axarr.axhline(y=unlabeledY[-2],color='#000000')
    axarr.axvline(x=unlabeledX[1],color='#000000')
    axarr.axvline(x=unlabeledX[-2],color='#000000')

    axarr.set_xticks(numpy.arange(0,xBinEdgeLen))
    axarr.set_xticklabels(xBinEdges[0:xBinEdgeLen],rotation=270)
    axarr.tick_params(axis='x',labelsize=8)
    axarr.set_xlabel(dimLongNames[varPos1])

    axarr.set_yticks(numpy.arange(0,yBinEdgeLen))
    axarr.set_yticklabels(yBinEdges[0:yBinEdgeLen])
    axarr.tick_params(axis='y',labelsize=8)
    axarr.set_ylabel(dimLongNames[varPos2],rotation=90)

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

def plot2DHistogramLatLonMap(histToPlot,varPos1,varPos2,pdfFile):

    fig,axarr = plt.subplots(9,10,figsize=(24,12),dpi=300)
    fig.tight_layout(h_pad=-1,w_pad=-1,rect=(-0.005,0.023,0.92,0.968))

    A = fig.add_axes([0.01,0.01,0.9,0.98],zorder=-1)

    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=0,urcrnrlon=360,ax=A)
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    m.drawmapboundary(fill_color='aqua')
    m.drawparallels(numpy.arange(-70,90,20),labels=[False,False,False,False])
    m.drawmeridians(numpy.arange(36,360,36),labels=[False,False,False,False])

    xBinEdges = dimEdges[varPos1]
    yBinEdges = dimEdges[varPos2]
    xBinEdgeLen = xBinEdges.size
    yBinEdgeLen = yBinEdges.size
    unlabeledX = numpy.arange(0,xBinEdgeLen,1)
    unlabeledY = numpy.arange(0,yBinEdgeLen,1)

    minPower = -7
    maxPower = 0
    levTicks = numpy.arange(minPower,maxPower+1,1)

    for yy in range(0,9):
        for xx in range(0,10):

            histAtLatLon = histToPlot
            if(varPos1 < 2 and varPos2 < 2):
                histAtLatLon = histToPlot[:,:,xx,yy]
            elif(varPos1 < 2 or varPos2 < 2):
                histAtLatLon = histToPlot[:,xx,yy,:]
            else:
                histAtLatLon = histToPlot[xx,yy,:,:]

            histLog = numpy.log10(1.0*histAtLatLon/(numpy.sum(histAtLatLon)))
            histLog = numpy.where(histLog == -1.*numpy.inf,minPower,histLog)

            p = None
            if(varPos1 > varPos2):
                p = axarr[9-1-yy,xx].pcolor(unlabeledX,unlabeledY,histLog,vmin=minPower,vmax=maxPower,cmap='Spectral_r')
            else:
                p = axarr[9-1-yy,xx].pcolor(unlabeledX,unlabeledY,numpy.transpose(histLog),vmin=minPower,vmax=maxPower,cmap='Spectral_r')

            axarr[9-1-yy,xx].axhline(y=unlabeledY[1],color='#000000')
            axarr[9-1-yy,xx].axhline(y=unlabeledY[-2],color='#000000')
            axarr[9-1-yy,xx].axvline(x=unlabeledX[1],color='#000000')
            axarr[9-1-yy,xx].axvline(x=unlabeledX[-2],color='#000000')

            axarr[9-1-yy,xx].set_xticks([])

            axarr[9-1-yy,xx].set_yticks([])

    cax = fig.add_axes([0.94,0.05,0.02,0.9])
    cbar = fig.colorbar(p,cax=cax,ticks=numpy.arange(minPower,maxPower+1,1))

    fig.suptitle(dimLongNames[varPos2]+' (Y) versus '+dimLongNames[varPos1]+' (X)')

    pdfFile.savefig()

    plt.close()

if (__name__ == "__main__"):
    main()
