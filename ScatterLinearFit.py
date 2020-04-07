#from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))
#%matplotlib inline

from os import getcwd
import numpy as np
import math

from statistics import stdev, mean

from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator

figureOutput='PlotTemplate'

SMALL_SIZE = 10
MEDIUM_SIZE = 10
BIGGER_SIZE = 10

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
  

rcParams['axes.linewidth']   = 2
rcParams['figure.figsize']   = 4,4
rcParams['figure.dpi']       = 600

rcParams['xtick.major.size']  = 7
rcParams['xtick.major.width'] = 1.5
rcParams['xtick.minor.size']  = 5
rcParams['xtick.minor.width'] = 1
rcParams['ytick.major.size']  = 7
rcParams['ytick.major.width'] = 1.5
rcParams['ytick.minor.size']  = 5
rcParams['ytick.minor.width'] = 1


rcParams['text.usetex']      = 'false'
rcParams['font.family']      = 'arial' #'arial'
#rcParams['font.serif']       = 'Times New Roman'

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm']      = 'serif'
rcParams['mathtext.it']      = 'serif:italic'
rcParams['mathtext.bf']      = 'serif:bold'
rcParams['mathtext.default'] = 'it'

print(getcwd())

markers=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
colors=['b','g','r','c','m','y','k']
lines=['solid',(0, (1,1)),(0, (5,1)),(0, (3, 2, 1, 2)),(0, (3, 2, 1, 2, 1, 2))]

##legend texts
elements=['Ba','Be','Ca','Co','Cu','Fe','Mg','Mn','Ni','Pb','Sr','Zn']
d=0.07
dx=[-d*1.5, -d*0.5, -d,      d*0.5,  d*0.3, -0.5*d, -1.8*d,  0.5*d, -0.3*d, 0.5*d, -0.4*d,-d*1.2]
dy=[-d*0.2,  d*0.6,  d*0.5, -d*0.2, -d*0.7,  0.6*d, -0.3*d, -0.3*d, -1.3*d, 0,    -1.2*d, d*0.5]
##data
ionicRadiusShannon=[1.3500,0.4500,1.0000,0.7450,0.7300,0.7800,0.7200,0.8300,0.6900,1.1900,1.1800,0.7400]
dftBondLength=[2.7648,3.6158497513771559/2,2.3824,2.0991,2.0650,2.1755,2.1071,2.2078,2.0510,2.5722,2.5666,2.1364]
ionicRadiusShannon=np.array(ionicRadiusShannon)
dftBondLength=np.array(dftBondLength)

##data fitting
cmin, cmax = min(ionicRadiusShannon), max(ionicRadiusShannon)
pfit, stats = np.polynomial.Polynomial.fit(ionicRadiusShannon, dftBondLength, 1, full=True, window=(cmin, cmax),domain=(cmin, cmax))

##plotter
for i in range(len(elements)):
    if (i==3 or i==4 or i==6 or i==8 or i==11):
        tempcolor  = colors[i%len(colors)]
        templine   = lines[i%len(lines)]
        tempmarker = markers[i%len(markers)]
        plt.scatter(ionicRadiusShannon[i],dftBondLength[i],marker=tempmarker,label=elements[i],s=30,c=tempcolor)
        plt.text(ionicRadiusShannon[i]+dx[i],dftBondLength[i]+dy[i],elements[i],color=tempcolor)
plt.plot(np.array([cmin,cmax]),pfit(np.array([cmin,cmax])),linestyle=lines[2],color='k')

ax=plt.axes()
ax.tick_params(axis='both',which='both',direction='in',colors='k',\
               bottom=True,top=True,left=True,right=True,\
               #labelbottom=True, labeltop=True, labelleft=True, labelright=True,\
               labelrotation=0)
#ax.tick_params(axis='y',which='major',labelrotation=90,direction='inout',length=12,width=6)
#ax.tick_params(axis='y',which='minor',direction='in',length=8,width=2)
plt.xlabel('Shannon ionic radius [$\mathrm{\AA}$]')
plt.xlim(0.4,1.42)
plt.xticks(np.arange(0.4, 1.42*1.001, step=0.2))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))

plt.ylabel(r'$\mathit{A-}$O bond length [$\mathrm{\AA}$]')
plt.ylim(1.78,2.8)
plt.yticks(np.arange(1.8, 2.8*1.001, step=0.2))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))


#plt.legend(frameon=0)
plt.savefig(figureOutput+'.pdf',transparent=True, bbox_inches='tight', pad_inches=0.01)
plt.savefig(figureOutput+'.png',transparent=True, bbox_inches='tight', pad_inches=0.01)
#plt.show()
