#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Still need to improve the code to print all plots in pdf
#https://matplotlib.org/examples/color/named_colors.html for python colors
'''
The DOSCAR file has the following organization:
1) First there is a header.
2) The important line from this is the 6th line which contains
   a) energy range used for the density of states,
   b) the number of bins, and
   c) the Fermi energy.
3) In the case of non-spin-polarized calculation, you will see 3 columns:
   a) E,
   b) TDOS, and
   c) integrated TDOS.
4) In the case of a spin-polarized calculation, you will see 5 columns:
   a) E,
   b) TDOS up,
   c) TDOS down,
   d) integrated TDOS up, and
   e) integrated TDOS down.
5) if LORBIT = 11 then after this, at line NEDOS + 8, you'll have
the beginning of the partial density of states.
   a) Every NEDOS+1 lines, there will be another entry.
   b) Each of these blocks corresponds to each of the atoms listed in your POSCAR file.
   c) Between each the line of the header with energy ranges, etc. is repeated.
   d) The columns are ordered similarly to the PROCAR file labels.
   e) So the first column is the energy followed by: E, s, p_y, p_z, p_x, d_xy, d_yz, d_z2, d_xz, d_x2.
   f) Again, if it is spin-polarized, then there will be an up and down column for each of these.
   g) This would give 10 total columns or 19 if spin-polarized.
If you have included elements with POTCARs that specifically include f-orbitals,
then there will also be another 7 columns that deal with each of these (14 for spin polarization).
Again the order is similar to the PROCAR, except that now they are labeled with the lm number
(-3, -2, ..., 3).
Now, since each atom is listed separately, you may need to sum up atoms to see the full
contribution from that species (for example to see the O 2p band in an oxide you would have to
sum all three p orbitals for all O atoms). Or if you want to look at a layer in a material, you'll
need to manually sum up each atoms contribution in the layer. If you are curious about crystal
field effects on a given atom, the individual d-orbitals are given with respect to the Cartesian
coordinates (not the lattice constants). So if you have, say an octahedrally coordinated transition
metal, you may need to rotate the cell such that the octahedron is oriented along the Cartesian axes
(or close to it if it is distorted).
'''


# In[2]:


#from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))

#from IPython.display import display, Math, Latex
import os
import re
import csv
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.backends.backend_pdf import PdfPages

import ase.io
import ase.neighborlist
from   ase import atom

print(os.getcwd())

DOSCARfileName='./DOSCAR'    #dos data is present here
CONTCARfileName='./CONTCAR'  #required for counting and naming species etc.
OUTDIR='./ldos'
if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)
pdosPlots='pdosPlots.pdf'                     #output file

normTotal=10000
#min=-113.879
#Emax=12.031
#DeltaE=0.01
#Fermi=12.0361
#totElectronsCheck=912

system      = ase.io.read(CONTCARfileName,format="vasp")

#get_ipython().run_line_magic('matplotlib', 'inline')
#%matplotlib auto
#font = {'family' : 'DejaVu Sans',
#        'size'   : 18}

#matplotlib.rc('font', **font)
#matplotlib.rc('xtick', labelsize=18)
#matplotlib.rc('ytick', labelsize=18)

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=MEDIUM_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#global Elist
#global pdostot
#global temppdos


# In[3]:


DOSCARfile=open(DOSCARfileName,'r')
i=0
sumUP=0
sumDOWN=0
for line in DOSCARfile :
    if (i==0):
        line=line.split()
        nAtm=int(line[0])
    if (i==5):
        #print(i-5,line)
        line=line.split()
        eMax=float(line[0])
        eMin=float(line[1])
        nBins=int(line[2])
        eFermi=float(line[3])
    if (i==6):
        line=line.split()
        if(len(line)==5):
            print("Spin polarized dos")
        else:
            print("Not a spin polarized dos")
        sumUP   += float(line[1])
        sumDOWN += float(line[3])
    if (i>6 and i<=6+300):
        line=line.split()
        #print(line[0],line[1],line[3])
        if (float(line[0])<=eFermi):
            sumUP   += float(line[1])
            sumDOWN += float(line[2])
            #print(line[0],"\t",sumUP,"\t",sumDOWN,"\t",line[3],"\t",line[4])
    i=i+1
if(nAtm+1-(i-5)/(nBins+1)==0):
    print("Everything is allright")
else:
    print("check your DOSCAR")
print(sumUP, sumDOWN, nBins)
DOSCARfile.close()


# In[4]:


##creating pdf
with PdfPages(pdosPlots) as pdf:
    for i in range(0,nAtm):                 ##i is atom index 1 to 6
        print("plotting atom",i)
        DOSCARfile=open(DOSCARfileName,'r')
        it=0
        atm=i+1
        dos=np.zeros((33,nBins-1))            ##33 columns by so many rows
        totaldos=np.zeros((5,nBins))
        for line in DOSCARfile :
            if ( it>=7+(0+0)*(nBins+1)   and it<=4+(0+1)*(nBins+1) ):
                j=it-7
                line=line.split()
                for k in range(0,5):
                    totaldos[k,j]=float(line[k])
            if ( it>=7+(atm+0)*(nBins+1) and it<=4+(atm+1)*(nBins+1) ):
                j=it-7-atm*(nBins+1)
                line=line.split()
                for k in range(0,19):
                    dos[k,j]=float(line[k])
            it=it+1
        DOSCARfile.close()
        #print(dos)
        atomLegend=str(system.get_chemical_symbols()[atm-1])+str(atm)
        plt.rcParams['figure.figsize'] = [11, 8.5]
        lineSUp, =plt.plot(dos[0], dos[1] ,linestyle='-', color='b', linewidth=1.6)
        lineSDn, =plt.plot(dos[0],-dos[2] ,linestyle='-', color='b', linewidth=1.6)

        lineTUp, =plt.plot(totaldos[0], totaldos[1]/normTotal ,linestyle=':', color='black', linewidth=1.6)
        lineTDn, =plt.plot(totaldos[0],-totaldos[2]/normTotal ,linestyle=':', color='black', linewidth=1.6)

        plt.legend([lineSUp,lineTUp],[atomLegend+'$-s$','$Total$/'+str(normTotal)])
        #plt.suptitle('atom number='+str(atm), fontsize=12)

        matplotlib.pyplot.axvline(x=eFermi,linestyle=':',color='brown',linewidth=0.8)
        plt.ylabel('LDOS [arb. units]')
        #plt.ylim((-6, 6))
        plt.xlabel('$E-E_f$ [eV]')
        plt.xlim((dos[0,0], dos[0,-1]))
        pdf.savefig()
        #plt.show()
        plt.close()

        linePxUp, =plt.plot(dos[0], dos[3] ,linestyle='-', color='g', linewidth=1.6)
        linePxDn, =plt.plot(dos[0],-dos[4] ,linestyle='-', color='g', linewidth=1.6)

        linePyUp, =plt.plot(dos[0], dos[5] ,linestyle='-', color='r', linewidth=1.6)
        linePyDn, =plt.plot(dos[0],-dos[6] ,linestyle='-', color='r', linewidth=1.6)

        linePzUp, =plt.plot(dos[0], dos[7] ,linestyle='-', color='m', linewidth=1.6)
        linePzDn, =plt.plot(dos[0],-dos[8] ,linestyle='-', color='m', linewidth=1.6)

        linePTUp, =plt.plot(dos[0], dos[3]+dos[5]+dos[7], linestyle='-.', color='black', linewidth=1.6)
        linePTDn, =plt.plot(dos[0],-dos[4]-dos[6]-dos[8], linestyle='-.', color='black', linewidth=1.6)

        lineTUp,  =plt.plot(totaldos[0], totaldos[1]/normTotal ,linestyle=':', color='black', linewidth=1.6)
        lineTDn,  =plt.plot(totaldos[0],-totaldos[2]/normTotal ,linestyle=':', color='black', linewidth=1.6)

        plt.legend([linePxUp,linePyUp,linePzUp,linePTUp,lineTUp],
                   [atomLegend+'$-p_x$',atomLegend+'$-p_y$',atomLegend+'$-p_z$',atomLegend+'$-p_{Total}$','$Total$/'+str(normTotal)])
        #plt.suptitle('atom number='+str(atm), fontsize=12)

        matplotlib.pyplot.axvline(x=eFermi,linestyle=':',color='brown',linewidth=0.8)
        plt.ylabel('LDOS [arb. units]')
        #plt.ylim((-6, 6))
        plt.xlabel('$E-E_f$ [eV]')
        plt.xlim((dos[0,0], dos[0,-1]))
        pdf.savefig()
        #plt.show()
        plt.close()

        lineDxyUp,  =plt.plot(dos[0], dos[9] ,linestyle='-', color='y', linewidth=1.6)
        lineDxyDn,  =plt.plot(dos[0],-dos[10],linestyle='-', color='y', linewidth=1.6)

        lineDyzUp,  =plt.plot(dos[0], dos[11],linestyle='-', color='k', linewidth=1.6)
        lineDyzDn,  =plt.plot(dos[0],-dos[12],linestyle='-', color='k', linewidth=1.6)

        lineDzxUp,  =plt.plot(dos[0], dos[13],linestyle='-', color='c', linewidth=1.6)
        lineDzxDn,  =plt.plot(dos[0],-dos[14],linestyle='-', color='c', linewidth=1.6)

        lineDx2y2Up,=plt.plot(dos[0], dos[15],linestyle='-', color='cadetblue', linewidth=1.6)
        lineDx2y2Dn,=plt.plot(dos[0],-dos[16],linestyle='-', color='cadetblue', linewidth=1.6)

        lineDz2Up,  =plt.plot(dos[0], dos[17],linestyle='-', color='firebrick', linewidth=1.6)
        lineDz2Dn,  =plt.plot(dos[0],-dos[18],linestyle='-', color='firebrick', linewidth=1.6)

        lineDTUp,   =plt.plot(dos[0], dos[ 9]+dos[11]+dos[13]+dos[15]+dos[17],linestyle='-.', color='black', linewidth=1.6)
        lineDTDn,   =plt.plot(dos[0],-dos[10]-dos[12]-dos[14]-dos[16]-dos[18],linestyle='-.', color='black', linewidth=1.6)

        lineTUp,    =plt.plot(totaldos[0], totaldos[1]/normTotal ,linestyle=':', color='black', linewidth=1.6)
        lineTDn,    =plt.plot(totaldos[0],-totaldos[2]/normTotal ,linestyle=':', color='black', linewidth=1.6)

        plt.legend([lineDxyUp,lineDyzUp,lineDzxUp,lineDx2y2Up,lineDz2Up,lineDTUp,lineTUp],
                   [atomLegend+'$-d_{xy}$',atomLegend+'$-d_{yz}$',atomLegend+'$-d_{zx}$',
                    atomLegend+'$-d_{x^2-y^2}$',atomLegend+'$-d_{z^2}$',atomLegend+'$-d_{Total}$','$Total$/'+str(normTotal)])
        #plt.suptitle('atom number='+str(atm), fontsize=12)

        matplotlib.pyplot.axvline(x=eFermi,linestyle=':',color='brown',linewidth=0.8)
        plt.ylabel('LDOS [arb. units]')
        #plt.ylim((-6, 6))
        plt.xlabel('$E-E_f$ [eV]')
        plt.xlim((dos[0,0], dos[0,-1]))
        pdf.savefig()
        #plt.show()
        plt.close()


        lineF_3Up,  =plt.plot(dos[0], dos[19] ,linestyle='-', color='y', linewidth=1.6)
        lineF_3Dn,  =plt.plot(dos[0],-dos[20] ,linestyle='-', color='y', linewidth=1.6)

        lineF_2Up,  =plt.plot(dos[0], dos[21] ,linestyle='-', color='k', linewidth=1.6)
        lineF_2Dn,  =plt.plot(dos[0],-dos[22] ,linestyle='-', color='k', linewidth=1.6)

        lineF_1Up,  =plt.plot(dos[0], dos[23] ,linestyle='-', color='c', linewidth=1.6)
        lineF_1Dn,  =plt.plot(dos[0],-dos[24] ,linestyle='-', color='c', linewidth=1.6)

        lineF0Up,  =plt.plot(dos[0], dos[25] ,linestyle='-', color='cadetblue', linewidth=1.6)
        lineF0Dn,  =plt.plot(dos[0],-dos[26] ,linestyle='-', color='cadetblue', linewidth=1.6)

        lineF1Up,  =plt.plot(dos[0], dos[27] ,linestyle='-', color='firebrick', linewidth=1.6)
        lineF1Dn,  =plt.plot(dos[0],-dos[28] ,linestyle='-', color='firebrick', linewidth=1.6)

        lineF2Up,  =plt.plot(dos[0], dos[29] ,linestyle='-', color='darkorange', linewidth=1.6)
        lineF2Dn,  =plt.plot(dos[0],-dos[30] ,linestyle='-', color='darkorange', linewidth=1.6)

        lineF3Up,  =plt.plot(dos[0], dos[31] ,linestyle='-', color='forestgreen', linewidth=1.6)
        lineF3Dn,  =plt.plot(dos[0],-dos[32] ,linestyle='-', color='forestgreen', linewidth=1.6)

        lineFTUp,   =plt.plot(dos[0], dos[19]+dos[21]+dos[23]+dos[25]+dos[27]+dos[29]+dos[31],linestyle='-.', color='black', linewidth=1.6)
        lineFTDn,   =plt.plot(dos[0],-dos[20]-dos[22]-dos[24]-dos[26]-dos[28]-dos[30]-dos[32],linestyle='-.', color='black', linewidth=1.6)

        lineTUp,    =plt.plot(totaldos[0], totaldos[1]/normTotal ,linestyle=':', color='black', linewidth=1.6)
        lineTDn,    =plt.plot(totaldos[0],-totaldos[2]/normTotal ,linestyle=':', color='black', linewidth=1.6)

        plt.legend([lineF_3Up,lineF_2Up,lineF_1Up,lineF0Up,lineF1Up,lineF2Up,lineF3Up,lineFTUp,lineTUp],
                   [atomLegend+'$-f_{-3}$',atomLegend+'$-f_{-2}$',atomLegend+'$-f_{-1}$',
                    atomLegend+'$-f_{0}$',atomLegend+'$-f_{1}$',atomLegend+'$-f_{2}$',atomLegend+'$-f_{3}$',
                    atomLegend+'$-f_{Total}$','$Total$/'+str(normTotal)])
        #plt.suptitle('atom number='+str(atm), fontsize=12)

        matplotlib.pyplot.axvline(x=eFermi,linestyle=':',color='brown',linewidth=0.8)
        plt.ylabel('LDOS [arb. units]')
        #plt.ylim((-6, 6))
        plt.xlabel('$E-E_f$ [eV]')
        plt.xlim((dos[0,0], dos[0,-1]))
        pdf.savefig()
        #plt.show()
        plt.close()


# In[5]:


###Writing data to files

def writeTDOS(fName,totaldos,nBins):
    fileName=fName+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e} {2:12.4e}'.format(totaldos[0,i]-eFermi,totaldos[1,i],-totaldos[2,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".dos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(totaldos[0,i]-eFermi,totaldos[1,i]+totaldos[2,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()

def writeDOS(fName,totaldos,nBins):

    #sssssssssssssssssssssssssssssssss
    fileName=fName+".01-s.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[1,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".02-s.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[2,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".33-s.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[1,i]+dos[2,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()

    #pxpxpxpxpxpxpxpxpxpxpxpxpxpxpxpxpxpxpx
    fileName=fName+".03-px.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[3,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".04-px.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[4,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".34-px.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[3,i]+dos[4,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #pypypypypypypypypypypypypypypypypypypy
    fileName=fName+".05-py.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[5,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".06-py.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[6,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".35-py.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[5,i]+dos[6,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #pzpzpzpzpzpzpzpzpzpzpzpzpzpzpzpzpzpzpz
    fileName=fName+".07-pz.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[7,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".08-pz.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[8,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".36-pz.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[7,i]+dos[8,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #ppppppppppppppppppppppppppppppppppppppp
    fileName=fName+".37-p.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[3,i]+dos[5,i]+dos[7,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".38-p.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[4,i]-dos[6,i]-dos[8,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".39-p.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[3,i]+dos[4,i]+dos[5,i]+dos[6,i]+dos[7,i]+dos[8,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()

    #dxydxydxydxydxydxydxydxydxydxydxydxydxy
    fileName=fName+".09-dxy.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[9,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".10-dxy.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[10,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".40-dxy.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[9,i]+dos[10,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #dyzdyzdyzdyzdyzdyzdyzdyzdyzdyzdyzdyzdyz
    fileName=fName+".11-dyz.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[11,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".12-dyz.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[12,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".41-dyz.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[11,i]+dos[12,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #dzxdzxdzxdzxdzxdzxdzxdzxdzxdzxdzxdzxdzx
    fileName=fName+".13-dzx.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[13,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".14-dzx.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[14,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".42-dzx.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[13,i]+dos[14,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #dx2y2dx2y2dx2y2dx2y2dx2y2dx2y2dx2y2dx2y2
    fileName=fName+".15-dx2y2.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[15,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".16-dx2y2.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[16,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".43-dx2y2.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[15,i]+dos[16,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #dz2dz2dz2dz2dz2dz2dz2dz2dz2dz2dz2dz2dz2dz2
    fileName=fName+".17-dz2.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[17,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".18-dz2.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".44-dz2.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[17,i]+dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #dddddddddddddddddddddddddddddddddddddddddddd
    fileName=fName+".45-d.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[ 9,i]+dos[11,i]+dos[13,i]+dos[15,i]+dos[17,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".46-d.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[10,i]-dos[12,i]-dos[14,i]-dos[16,i]-dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".47-d.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[ 9,i]+dos[10,i]+dos[11,i]+dos[12,i]+dos[13,i]+dos[14,i]+dos[15,i]+dos[16,i]+dos[17,i]+dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()


    #f-3
    fileName=fName+".19-f-3.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[19,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".20-f-3.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[20,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".48-f-3.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[9,i]+dos[10,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #f-2
    fileName=fName+".21-f-2.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[21,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".22-f-2.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[22,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".49-f-2.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[11,i]+dos[12,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #f-1
    fileName=fName+".23-f-1.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[23,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".24-f-1.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[24,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".50-f-1.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[13,i]+dos[14,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #f0
    fileName=fName+".25-f0.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[25,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".26-f0.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[26,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".51-f0.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[15,i]+dos[16,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #f1
    fileName=fName+".27-f1.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[27,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".28-f1.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[28,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".52-f1.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[17,i]+dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #f2
    fileName=fName+".29-f2.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[29,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".30-f2.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[30,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".53-f2.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[17,i]+dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #f3
    fileName=fName+".31-f3.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[31,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".32-f3.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[32,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".54-f3.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[17,i]+dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    #f-total
    fileName=fName+".55-f.up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[ 9,i]+dos[11,i]+dos[13,i]+dos[15,i]+dos[17,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".56-f.dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[10,i]-dos[12,i]-dos[14,i]-dos[16,i]-dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".57-f.tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[ 9,i]+dos[10,i]+dos[11,i]+dos[12,i]+dos[13,i]+dos[14,i]+dos[15,i]+dos[16,i]+dos[17,i]+dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()





    #atomatomatomatomatomatomatomatomatomatomatom
    fileName=fName+".58-up"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[1,i]+dos[3,i]+dos[5,i]+dos[7,i]+dos[ 9,i]+dos[11,i]+dos[13,i]+dos[15,i]+dos[17,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".59-dn"+".ldos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,-dos[2,i]-dos[4,i]-dos[6,i]-dos[8,i]-dos[10,i]-dos[12,i]-dos[14,i]-dos[16,i]-dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()
    fileName=fName+".60-tot"+".pdos"+".dat"
    fWrite = open(fileName,'w')
    for i in range(nBins):
        writeLine='{0:12.4e} {1:12.4e}'.format(dos[0,i]-eFermi,dos[1,i]+dos[2,i]+dos[3,i]+dos[4,i]+dos[ 5,i]+dos[6,i]+dos[7,i]+dos[8,i]+dos[9,i]+dos[10,i]+dos[11,i]+dos[12,i]+dos[13,i]+dos[14,i]+dos[15,i]+dos[16,i]+dos[17,i]+dos[18,i])+"\n"
        fWrite.write(writeLine)
    fWrite.close()

for i in range(nAtm):
    print("writing data of atom",i)
    if (i==0) :
        fName=OUTDIR+"/"+"000.total"
        writeTDOS(fName,totaldos,nBins-1)
    fName=OUTDIR+"/"+str('{0:03d}'.format(i+1))+str(system.get_chemical_symbols()[i])
    writeDOS(fName,dos,nBins-1)


# In[6]:


print("done")


# In[ ]:
