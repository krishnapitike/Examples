from matplotlib import pyplot as plt
%matplotlib inline
import numpy as np
import random
from os import getcwd

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


def gaussianSmoothing2D(x_vals,y_vals,sigma):
    y_smth = np.zeros(y_vals.shape) 
    for it in range(0,len(x_vals)):
        x_position      = x_vals[it]
        gaussian_kernel = np.exp(-(x_vals - x_position) ** 2 / (2 * sigma ** 2))
        gaussian_kernel = gaussian_kernel / np.sum(gaussian_kernel)
        y_smth[it]      = np.sum(y_vals * gaussian_kernel)
    return(y_smth)

#some example data
x= np.linspace(0, 2*np.pi, 500)
y = np.sin(x) 
#some random confidence interval
ci = 0.3*np.random.rand(len(x))
#smooothed random confidence interval
cis=gaussianSmoothing2D(x,ci,0.05)

#plot
fig, ax = plt.subplots()
ax.plot(x,y, color='r',label='y')
ax.fill_between(x, (y-ci), (y+ci), color='b', alpha=.2,label='confidence interval')
ax.fill_between(x, (y-cis), (y+cis), color='g', alpha=.4,label='smoothed')
plt.xlabel('x')
plt.ylabel('y')
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.tick_params(axis='both',which='both',direction='in',colors='k',\
               bottom=True,top=True,left=True,right=True,\
               #labelbottom=True, labeltop=True, labelleft=True, labelright=True,\
               labelrotation=0)
plt.legend(frameon=0)
plt.show()
