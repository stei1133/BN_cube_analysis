import numpy as np
import matplotlib.pyplot as plt

def readFile(filename):
   data=[]
   headercount=0   
   global NX,NY,NZ
   with open(filename) as f:
        for line in f:
             line = line.split() # to deal with blank 
             if line:            # lines (ie skip them)
                 try:
                       line = [float(i) for i in line]
                       if len(line)==6:
                            for i in line:
                                data.append(i)
                       elif len(line)==4:
                            if headercount==0:
                                 numofatoms=int(line[0])
                            elif headercount==1:
                                 NX=int(line[0])
                            elif headercount==2:
                                 NY=int(line[0])
                            elif headercount==3:
                                 NZ=int(line[0])
                            headercount+=1                            
                 except ValueError:
                       total=0
   return data

def yplanaraverage(imported_data):
     averagedpotential=np.zeros((NZ,NX))
     total=0
     for x in range(NX):
        for z in range(NZ):
          for i in np.linspace(x*NY*NZ+z,x*NZ*NY+z+NZ*NY,num=NY,endpoint=False):
              total+=imported_data[int(i)]
          avg=total/NY
          averagedpotential[z][x]=avg
          total=0     
     return averagedpotential
   
def makeLineProfile(data,z_index):
     linedata=change_in_potential[0:NX-1][int(NZ/2+.5)]
     plt.plot(linedata)
     plt.xlabel('X-position (Angstroms)')
     plt.ylabel('Electrostatic Potential (eV)')
     tmparrayx=[0, NZ/4,NZ/2,3*NZ/4, NZ]
     labelx=[round(tmparrayx[i]*zstep,1) for i in range(len(tmparrayx))]
     plt.xticks(tmparrayx,labelx)
     plt.show()

ydimension=4.33
zdimension=35
xdimension=15.12
count=0
imported_overall_data=[]
imported_overall_data=readFile('BN3up3down_potpp.cube')
overallpotential=yplanaraverage(imported_overall_data)
#subtract HF
imported_HF_data=readFile('BN3up3down_HF.cube')
HF_potential=yplanaraverage(imported_HF_data)
graph_sub_HFpotential=overallpotential-HF_potential
#subtract BN
imported_graphene_data=readFile('BN3up3down_Cell.cube')
Cell_potential=yplanaraverage(imported_graphene_data)
vacuumdiff=graph_sub_HFpotential[0][0]-Cell_potential[0][0]
change_in_potential=(graph_sub_HFpotential-Cell_potential)*13.6056980659
#make the plots
zstep=zdimension/NZ
ystep=ydimension/NY
xstep=xdimension/NX
#makeLineProfile(change_in_potential,NZ/2+.5)
levels=np.linspace(-1,1,num=21)
CS=plt.contourf(change_in_potential,
                levels,
                cmap='bwr',origin='lower',extend='both')
plt.xlabel('X-position (Angstroms)')
plt.ylabel('Z-position (Angstroms from the surface)')
ytickarray=[0, NZ/4,NZ/2,3*NZ/4,NZ]
yticklabel=[-1*round(NZ/2*zstep,2), -1*round(NZ/4*zstep,2), 0,
       round(NZ/4*zstep,2), round(NZ/2*zstep,2)]

plt.yticks(ytickarray, yticklabel)
xtickarray=[0, NX/4,NX/2,3*NX/4, NX]
xticklabel=[round(xtickarray[i]*xstep) for i in range(len(xtickarray))]
plt.xticks(xtickarray,xticklabel)
cbar=plt.colorbar(CS)
cbar.set_label('Electrostatic Potential (eV)')
plt.show()
