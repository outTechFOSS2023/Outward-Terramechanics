# -*- coding: utf-8 -*-
#*********************************************************************************************************
#*********************************************************************************************************
#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#*********************************************************************************************************
#*********************************************************************************************************
# This Python script is a typical adaptive reduced order modeling (ROM) of a cuttter blade (CB) test in GRC-3 regolith at a  dense state (n=0.315).
# NOTE: ALL EXTERNAL FILES (DEM GEOMETRY), AND OUTPUT FOLDERES ARE PROVIDED. THE CORRECT ADDRESSES
# MUST BE PROVIDED BY THE USER FROM THEIR OWN LOCAL MACHINES. 
# This simulation has the following characteristics:
# 1- The granular assembly is purely made of spherical grains to allow for high degrees of compaction at stress free condition.
# 2- DEM assembly size is 6cm(X-direction)*2.5cm(Y-direction)x7.68cm(Z-Direction). The simulation is considered as a full-width model (X direction)
# where the blade moves along the centerline of the GRC-3 assembly.
# 3- Wall boundaries and the cutter blade has been generated using YADE built-in geometry generation functions for boxes and therefore are not deformable. 
# 4- The penetration staep of the typical cutter blade test has been skipped to accelerate the simulation. Here, the cutter blade starts moving from outside the 
# boundary box.
# 5- The cutting proces takes place in +Z direction using the YADE built-in TranslationEngine. 
# 6- In order to account for the shape effects of real grain geometries, a sophisticated rolling/twisting resistance model (Hertz-Mindlin with rolling resistance--see YADE manual for activation of the rolling resistance) has been implemented.
# and used for both the regolith and box components. However, to avoid numerical instabilities, rolling resistance parameters have been disabled for the boundaries
# and the blade .
#**********************************************************************************************************
#**********************************************************************************************************
from yade import pack, export
from yade import pack, plot 
from yade import ymport
import sys
import os
import os.path
import numpy
import math

######################################################
###   DEFINING VARIABLES AND MATERIAL PARAMETERS   ###
######################################################
#### UNITS
## LENGTH: m
## VELOCITY: m/s
## TIME: s
## PRESSURE: pa
## DENSITY: kg/m3
deposFricDegree = 21.42 # CONTACT FRICTION 
normalDamp=0.7 # NORMAL VISCOUS DAMPING
shearDamp=0.7 # SHEAR VISCOUS DAMPING
youngSoil=0.75e8 # CONTACT STIFFNESS FOR REGOLITH
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.44 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densSoil=1801.11858943976 # BULK DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.4 # NMERICAL DAMPING
# NOTE: THE NUMERICAL STABILITY AND RESPNSE ACCURACY IS HIGHLY DEPENDENT UPON THE SIZE OF THE ACTIVE ZONE (THE LARGER, THE BETTER) AND THE BIN SIZE
# (THE SMALLER, THE BETTER). THEREFORE, THESE TWO PARAMETERS ARE SUBJECT TO SENSITIVITY ANALYSES.  
binSize = 1/3*0.0192 # BIN SIZE ADDED IN FRONT OF THE ACTIVE ZONE OF THE DEM MODEL
activeDistance = 2*0.0192 # ACTIVE DEM ZONE IN THE MODEL.
IniDistanceBladefromBoundary = 0.0 #INITIAL DISTANCE OF THE RIGHT BOUNDARY OF THE CUTTER BLADE FROM THE LEFT BOUNDARY OF THE COMPACTED BED (-Z). 
DifferenceWCenterActiveEnd = activeDistance - IniDistanceBladefromBoundary
HeightBlade=0.07 #CUTTER BLADE HEIGHT (Y DIRECTION)
WidthBlade=0.02 # CUTTER BLADE WIDTH (X DIRECTION)
ThicknessBlade=0.000020 #CUTTER BLADE THICKNESS (Z DIRECTION)
InitialPeneterationofBlade = 0.010 # PENETRATION DEPTH OF THE BLADE. NOTE: HERE, THE BLADE STARTS CUTTING FROM THE OUTSIDE OF THE BOUNDARY (-Z)
# AND THIS PARAMETER IS ONLY USED FOR ITS INITIAL PLACEMENT.
InitialDistanceofBladefromTopSoil = 0 # INITIAL DISTANCE OF THE BOTTOM OF THE CUTTER BLADE FROM THE SURFACE OF THE COMPACTED BED. SINCE THE ASSEMBLY
# IS COMPLETELY SMOOTH, THIS PARAMETER CAN BE SET TO ZERO.
HorizentalvelocityofBlade = 0.025 #CUTTING VELOCITY. IN ORDER TO ENSURE THE NUMERICAL STABILITY OF THE MODEL, THE VELOCITY MUST BE 
# LOW ENOUGH TO AVOID DYNAMIC EFFECTS. THE ACCURACY OF THE FORCE-DISPLACEMENT RESPONSE IS HIGHLY DEPENDENT
# ON THE CUTTER BLADE VELOCITY. HERE , THE CUTTING VELOCITY IS SET TO A HIGH VALUE FOR DEMONESTRATION PURPOSES ONLY. 

############################################
###   DEFINING VARIABLES AND MATERIALS   ###
############################################
SoilId=FrictMat(young=youngSoil,poisson=poissonSoil,frictionAngle=radians(deposFricDegree),density=densSoil)
O.materials.append(SoilId)


ContainerId=FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=radians(0),density=densContainer)
O.materials.append(ContainerId)


###################################
#####   CREATING GEOMETRIES   #####
###################################
## IMPORTING THE ALREADY COMPACTED BIN (USE YOUR OWN LOCAL MACHINE ADDRESS)
ymport.textClumps("../Adaptive Reduced Order Modeling/GRC-3/Spherical Grains/Importable/CSP-GRC3-B1-A-SphericalGrains-ClumpFormat.txt",scale=10,orientation=Quaternion((0, 1, 0),-pi/2),shift=Vector3(0,0,0),color=Vector3(0.0,1.0,0.0),material=SoilId)
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
for b in O.bodies:
  if isinstance(b.shape,Sphere):
    O.bodies.erase(b.id)

ymport.textClumps("../Adaptive Reduced Order Modeling/GRC-3/Spherical Grains/Importable/CSP-GRC3-B1-A-SphericalGrains-ClumpFormat.txt",scale=10,orientation=Quaternion((0, 1, 0),-pi/2),shift=Vector3(-minX,0,0),color=Vector3(0.0,1.0,0.0),material=SoilId)
minZC=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZC=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
NBin = int((maxZC-minZC)/(binSize))
for b in O.bodies:
  if isinstance(b.shape,Sphere):
    if ((NBin*binSize)<(b.state.pos[2]+b.shape.radius)):
      O.bodies.erase(b.id)

minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"maxX:",maxX,"minY:",minY,"maxY:",maxY,"minZ:",minZ,"maxZ:",maxZ)
## EXPORTING THE PROCESSED ASSEMBLY (USE YOUR OWN LOCAL MACHINE ADDRESS)
export.textClumps("../Adaptive Reduced Order Modeling/GRC-3/Spherical Grains/Importable/CSP-GRC3-B1-A-SphericalGrains-ClumpFormat_Sectioned.txt")
for b in O.bodies:
  if isinstance(b.shape,Sphere):
    O.bodies.erase(b.id)

mn,mx=Vector3(minX,minY,minZ),Vector3(maxX,0.1,maxZ)

## CREATE WALLS AROUND THE PACKING
walls=aabbWalls([mn,mx],thickness=0,material=ContainerId)
wallIds=O.bodies.append(walls)

O.bodies.append(utils.box(center=((maxX-minX)/2,(maxY+HeightBlade/2-InitialPeneterationofBlade),minZ-ThicknessBlade/2+IniDistanceBladefromBoundary),extents=(WidthBlade/2,HeightBlade/2,ThicknessBlade/2),material=ContainerId,fixed=True,color=(0.2,0.5,1),wire=False))
Blade=O.bodies[-1]
Blade.state.blockedDOFs='xyXYZ'

###################################
#####     BIN GENERATION      #####
###################################
Spheres = open('../Adaptive Reduced Order Modeling/GRC-3/Spherical Grains/Importable/CSP-GRC3-B1-A-SphericalGrains-ClumpFormat_Sectioned.txt',"r")
SpheresData = numpy.loadtxt(Spheres, usecols=range(5))
a = 1
NumberofGrains = []
for b in range(1,len(SpheresData)):
  if SpheresData[b,4] == SpheresData[b-1,4]:
    a=a+1
    if (b==(len(SpheresData)-1)):
      counter=numpy.array([a])
      NumberofGrains.append(counter)
  else:
    counter=numpy.array([a])
    NumberofGrains.append(counter)
    a = 1
NumberofGrains=numpy.array(NumberofGrains).reshape(-1)

k = 0
for i in range(0,NBin):
  globals()['Bin' + str(int(i))]=[]
  globals()['PosZ' + str(int(i))]= (i+1)*binSize


for i in range(0,len(NumberofGrains)):
  Data = []
  k = k + NumberofGrains[i]
  l=(k-NumberofGrains[i])
  SpheresVolumn = 0
  ConvexArea = 0
  for j in range (int(l),int(k)):
    X = SpheresData[j,0]
    Y = SpheresData[j,1]
    Z = SpheresData[j,2]
    R = SpheresData[j,3]
    ID = SpheresData[j,4]
    Data.append(numpy.array([float(X),float(Y),float(Z),float(R),'{:d}'.format(int(ID))]))
  Data = numpy.array(Data).reshape(-1,5)
  m=0
  while m<NBin:
    BinLimit = ((globals()['PosZ' + str(int(m))])*numpy.ones(len(Data))).reshape(-1,1)
    BinLimit = BinLimit.astype(float)
    Zpose = Data[:,2].astype(float)
    ZposminC = (BinLimit-Zpose)<=(binSize)
    ZposmaxC = (BinLimit-Zpose)>=0
    if (ZposminC).any() and (ZposmaxC).any():
      for n in range(0,len(Data)):
        globals()['Bin' + str(int(m))].append(numpy.array([Data[n,0],Data[n,1],Data[n,2],Data[n,3],Data[n,4]]))
      break
    m=m+1

for i in range(0,NBin):
  globals()['Bin' + str(int(i))]=numpy.array(globals()['Bin' + str(int(i))]).reshape(-1,5)

  with open(('./'+'Bin'+str(int(i))+'.txt'),'w') as f:   
    f.write('#format x_y_z_r_clumpId')
    for q in range(0,len(globals()['Bin' + str(int(i))])):
      f.write('\n')
      f.write("\n".join(" ".join(map(str, q)) for q in globals()['Bin' + str(int(i))][q].reshape(-1,5)))

for ini in range(0,int(activeDistance/binSize)):
  Bin="./Bin"+str(int(ini))+".txt"
  ymport.textClumps("./Bin"+str(int(ini))+".txt",shift=Vector3(0,0,0),material=SoilId)
  os.remove(Bin)

############################
###   DEFINING ENGINES   ###
############################
gravity=(0,-9.81,0)
## SAVING VTK (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
O.engines=[VTKRecorder(fileName='../Adaptive Reduced Order Modeling/GRC-3/Spherical Grains/VTK/CSP-GRC3-B1-A-SphericalGrains-n0315' ,recorders=['all'],iterPeriod=5000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()], label="collider"),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=normalDamp,betas=shearDamp,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
  ),
  TranslationEngine(ids=[Blade.id],translationAxis=(0.,0.,1.),velocity=HorizentalvelocityofBlade,label='BladeHorizentalMove',dead=True),   
  NewtonIntegrator(damping=numDamp,gravity=gravity),
  PyRunner(command='ImportFirstRigid()',iterPeriod=1,label='checker'),
  PyRunner(iterPeriod=1000,command='history()',label='recorder'),
  ]
O.dt = 1.0e-6

def DataImportRigid(number):

  Bin="./Bin"+str(int(number))+".txt"

## If file exists, delete it ##
  if os.path.isfile(Bin):
    ymport.textClumps("./Bin"+str(int(number))+".txt",shift=Vector3(0,0,0),material=SoilId,mask=3)
    global BinClump, BinSphere, ImportedBin
    BinClump = []
    BinSphere = []
    for b in O.bodies:
      if isinstance(b.shape,Sphere):
        if b.isClumpMember:
          if b.mask==3:
             BinClump.append(b.clumpId)
             BinSphere.append(b.id)
    for n in BinSphere:
      if O.bodies[n].isClumpMember:
        O.bodies[n].groupMask=1
    for o in BinClump:
      if O.bodies[o].isClump:
        O.bodies[o].dynamic=False
        O.bodies[o].state.blockedDOFs='xyzXYZ' 
    os.remove(Bin)

  ImportedBin = int(number)

  return ImportedBin

def ImportFirstRigid():
  global ImportedBin
  ImportedBin = DataImportRigid(int(activeDistance/binSize))
  checker.command='BladeHorizentalMove()'

def BladeHorizentalMove():

  Blade_Center_PosZ=Blade.state.pos[2]
  
  for i in range(int(IniDistanceBladefromBoundary/binSize)+1,NBin):
    temBinClumpActive=[]

    if abs(Blade_Center_PosZ-(i*binSize-minZ))<=1e-5: #(i+1)
      for b in O.bodies:
        if isinstance(b.shape,Sphere):
          if b.isClumpMember:
            if (((i-int(IniDistanceBladefromBoundary/binSize))*binSize-minZ)<b.state.pos[2]<=((i+int(DifferenceWCenterActiveEnd/binSize))*binSize-minZ)):
              temBinClumpActive.append(b.clumpId)

      temBinClumpActive=numpy.unique(numpy.array(temBinClumpActive))
      for r in temBinClumpActive:
        if O.bodies[int(r)].isClump:
          O.bodies[int(r)].dynamic=True
      number= i+int(DifferenceWCenterActiveEnd/binSize)
      ImportedBin= DataImportRigid(int(number))
      
  for i in range(int(IniDistanceBladefromBoundary/binSize)+1,NBin):
    temBinClumpDeactive=[]
    temBinSphereDeactive=[]

    if abs(Blade_Center_PosZ-(i*binSize-minZ))<=1e-5: #0.75*(i)
      for b in O.bodies:
        if isinstance(b.shape,Sphere):
          if b.isClumpMember:
            if (b.state.pos[2]<=((i-int(IniDistanceBladefromBoundary/binSize))*binSize-minZ)):
              if b.clumpId != Blade.id:
                temBinClumpDeactive.append(b.clumpId)
                temBinSphereDeactive.append(b.id)

      temBinSphereDeactive=numpy.unique(numpy.array(temBinSphereDeactive))
      temBinClumpDeactive=numpy.unique(numpy.array(temBinClumpDeactive))
      for p in temBinSphereDeactive:
        if O.bodies[int(p)].isClumpMember:
          O.bodies[int(p)].groupMask=4
      for r in temBinClumpDeactive:
        if O.bodies[int(r)].isClump:
          O.bodies[int(r)].dynamic=False
          O.bodies[int(r)].state.blockedDOFs='xyzXYZ'

      for b in O.bodies:
        if isinstance(b.shape,Sphere):
          if (b.state.pos[2]<=((i-1-int(IniDistanceBladefromBoundary/binSize))*0.5*binSize-minZ)):
            O.bodies.erase(b.id)

  collider.avoidSelfInteractionMask = 4
  #O.engines[4].dead=True
  O.engines[4].dead=False
  print ("Blade_Center_PosZ+ThicknessBlade/2:",Blade_Center_PosZ+ThicknessBlade/2)
  if Blade_Center_PosZ>=5/6*maxZ:
    O.engines[4].dead=True 
    O.pause()

def history():
  global Fx,Fy,Fz,Dx,Dy,Dz
  Fx=0.01*abs(O.forces.f(Blade.id)[0])
  Fy=0.01*abs(O.forces.f(Blade.id)[1])
  Fz=0.01*abs(O.forces.f(Blade.id)[2])
  Dx=Blade.state.pos[0]
  Dy=Blade.state.pos[1]
  Dz=Blade.state.pos[2]+ThicknessBlade/2
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
  ## SAVING NUMERICAL FORCE-DISPLACEMENT (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
  plot.saveDataTxt('../Adaptive Reduced Order Modeling/GRC-3/Spherical Grains/Numerical Output/CSP-GRC3-B1-A-SphericalGrains-n0315.txt')
## declare what is to plot. "None" is for separating y and y2 axis
plot.plots={'Dz':('Fx','Fy','Fz')}
plot.plot()
# Uncomment for non-interactive simulation
#O.run(2000000000000,True)

