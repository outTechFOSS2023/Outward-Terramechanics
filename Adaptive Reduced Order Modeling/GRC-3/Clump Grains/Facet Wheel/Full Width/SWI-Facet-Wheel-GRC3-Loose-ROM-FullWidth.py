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
# This Python script is a typical adaptive reduced order modeling (ROM) of a smooth wheel soil interaction test (SWI) over GRC-3 regolith 
# at a  loose state (n=0.570).
# NOTE: ALL EXTERNAL FILES (DEM GEOMETRY), AND OUTPUT FOLDERES ARE PROVIDED. THE CORRECT ADDRESSES
# MUST BE PROVIDED BY THE USER FROM THEIR OWN LOCAL MACHINES. 
# This simulation has the following characteristics:
# 1- The granular assembly is  made of clump grains per the Outward's grain generation scheme.
# 2- DEM assembly size is 10mm(X-direction)*6.0mm(Y-direction)x20mm(Z-Direction). The simulation is considered as a full-width model (X direction)
# where the wheel penetrates and rolls  along the centerline of the GRC-3 assembly.
# 3- Wall boundaries and the wheel has been generated using YADE built-in geometry generation functions for facets and therefore is not deformable.
# 4- Wheel movement consists of two stages: (a) Penetration of the wheel to 1.0mm from the surface of the granular assembly (-Y direction), and (b) Rolling 
# thorough the regolith in +Z direction.
# 5- The penetration and rolling/movement  proceses have been modeled  using the YADE built-in TranslationEngine and CombinedKinematicEngine respectively. 
# 6- The YADE built-in Hertz-Mindlin constituitive model has been used.
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

############################################
###           INPUT PARAMETERS           ###
############################################
#### UNITS
## LENGTH: m
## VELOCITY: m/s
## ANGULAR VELOCITY: rad/s
## TIME: s
## MASS: kg
## PRESSURE: pa
## DENSITY: kg/m3
## WHEEL GEOMETRY
widthWheel=0.00125 # WIDTH OF THE WHEEL (THIN WHEEL ASSUMPTION)
Wheel_R=0.0025 # RADIUS OF THE WHEEL
Wheel_AdditionalNormalForce=-19.6 # AN ARBITRARY NORMAL FORCE APPLIED TO THE CENTER OF THE WHEEL
# REPRESENTING THE WHEIGHT

##WHELL MOVEMENT
slipRatio=0.3 # THE SLIP RATIO (S) SHOWS THE WHEEL SLIP ACCORDING TO THE LONGITUDINAL VELOCITY OF THE VEHICLE AND ANGULAR VELOCITY OF THE WHEEL. 
# AND 0 ≤ |S| ≤ 1. When THE ANGULAR VELOCITY= 0 THEN |S| = 1 AND THE WHEEL LOCK-UP WILL OCCUR. 
velocity=0.01 #HORIZONTAL VELOCITY OF THE WHEEL MOVEMENT. IN ORDER TO ENSURE THE NUMERICAL STABILITY OF THE MODEL, THE VELOCITY MUST BE 
# LOW ENOUGH TO AVOID DYNAMIC EFFECTS. THE ACCURACY OF THE FORCE-DISPLACEMENT RESPONSE IS HIGHLY DEPENDENT
# ON THE THIS VELOCITY. HERE , THE CUTTING VELOCITY IS SET TO A HIGH VALUE FOR DEMONESTRATION PURPOSES ONLY. 
angularVelocity=velocity/((1-slipRatio)*Wheel_R) # THE ROLLING SEQUENCE HAVE BEEN MODELED USING THE YADE CombinedKinematicEngin. HERE, THE  ROTATION
# TAKES PLACE AND AN ANGULAR VELOCITY MUST BE PROVIDED FOR RotaionEngine OF THIS CombinedKinematicEngin.
gravityAcc=-9.81 #GRAVITY ACCELERATION
deposFricDegree = 40# CONTACT FRICTION
normalDamp=0.7 # NORMAL VISCOUS DAMPING
shearDamp=0.7 # SHEAR VISCOUS DAMPING
youngSoil=0.7e8# CONTACT STIFFNESS FOR SOIL
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.3 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densSoil=2650 # DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.4 # NUMERICAL DAMPING
binSize=0.4*Wheel_R  # BIN SIZE ADDED IN FRONT OF THE ACTIVE ZONE OF THE DEM MODEL
activeDistance=0.01 # ACTIVE DEM ZONE IN THE MODEL.
IniDistanceWheelfromBoundary=0.004 #INITIAL DISTANCE OF THE WHEEL CENTER FROM THE LEFT BOUNDARY DURING PENETRATION
DifferenceWCenterActiveEnd=activeDistance - IniDistanceWheelfromBoundary
InitialPeneterationofWheel=0.001 # PENETRATION DEPTH OF THE WHEEL BEFORE THE START OF THE ROLLING PROCESS
IniWheelDisfromMaxY=0.0001  # INITIAL OFFSET OF THE BOTTOM OF THE WHEEL FROM THE TOP OF THE SAMPLE

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
## CLUMPS
ymport.textClumps("../Terramechanics Delivarable Codes/Adaptive Reduced Order Modeling/GRC-3/Clump Grains/Facet Wheel/Full Width/Importable/CSP-GRC3-B1-A-ClumpGrains-FullWidth-n057.txt",shift=Vector3(0,0,0),material=SoilId)
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
export.textClumps("../Terramechanics Delivarable Codes/Adaptive Reduced Order Modeling/GRC-3/Clump Grains/Facet Wheel/Full Width/Importable/CSP-GRC3-B1-A-ClumpGrains-FullWidth-n057_Sectioned.txt")
for b in O.bodies:
  if isinstance(b.shape,Sphere):
    O.bodies.erase(b.id)

mn,mx=Vector3(minX,minY,minZ),Vector3(maxX,0.01,maxZ)

## CREATE WALLS AROUND THE PACKING
walls=aabbWalls([mn,mx],thickness=0,material=ContainerId)
wallIds=O.bodies.append(walls)
for b in O.bodies:
  if b.state.pos[1]>0.009:
    O.bodies.erase(b.id)

XcenterofWheel = (maxX-minX)/2
YcenterofWheel = Wheel_R+maxY+iniWheelDisfromMaxY
ZcenterofWheel = iniDistanceWheelfromBoundary+minZ

Wheel1 = geom.facetCylinder((XcenterofWheel,YcenterofWheel,ZcenterofWheel),radius=Wheel_R,height=widthWheel,orientation=Quaternion((0,1,0),-pi/2),wallMask=7,segmentsNumber=200,material=ContainerId,wire=True,color=Vector3(255,255,255))
O.bodies.append(Wheel1)
idsr = [w.id for w in Wheel1]
facets = [b for b in O.bodies if isinstance(b.shape,Facet)] # list of facets in simulation
for b in facets:
  O.forces.setPermF(b.id,(0,Wheel_AdditionalNormalForce/800,0))

###################################
#####     BIN GENERATION      #####
###################################
Spheres = open('../Terramechanics Delivarable Codes/Adaptive Reduced Order Modeling/GRC-3/Clump Grains/Facet Wheel/Full Width/Importable/CSP-GRC3-B1-A-ClumpGrains-FullWidth-n057_Sectioned.txt',"r")
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
#print(NumberofGrains)
#print(len(NumberofGrains))

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
## SAVING VTK (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
O.engines=[VTKRecorder(fileName='../Terramechanics Delivarable Codes/Adaptive Reduced Order Modeling/GRC-3/Clump Grains/Facet Wheel/Full Width/VTK/CSP-GRC3-B1-A-ClumpGrains-FullWidth-n057' ,recorders=['all'],iterPeriod=5000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Facet_Aabb()], label="collider"),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=normalDamp,betas=shearDamp,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
  ),
  CombinedKinematicEngine(ids=idsr,label='combEngine') + TranslationEngine(translationAxis=(0,-1,0),velocity=velocity) +\
  RotationEngine(rotationAxis=(1,0,0), angularVelocity=0, rotateAroundZero=True, zeroPoint=(XcenterofWheel,YcenterofWheel,ZcenterofWheel)),
  NewtonIntegrator(damping=numDamp,gravity=(0,gravityAcc,0)),
  PyRunner(command='ImportFirstRigid()',iterPeriod=1,label='checker'),
  PyRunner(iterPeriod=1000,command='history()',label='recorder'),
  PyRunner(iterPeriod=1000,command='eraseOffs()'),
  ]
O.dt = 1e-7

# get TranslationEngine and RotationEngine from CombinedKinematicEngine
transEngine, rotEngine = combEngine.comb[0], combEngine.comb[1]

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
  checker.command='WheelPenetration()'


def WheelPenetration():
  transEngine.velocity = 7*velocity
  rotEngine.zeroPoint += Vector3(0,-1,0)*transEngine.velocity*O.dt
  if (rotEngine.zeroPoint[1]-Wheel_R)<=(maxY-initialPeneterationofWheel):
    checker.command='WheelRolling()'

def WheelRolling():
  transEngine.translationAxis=(0,0,1)  
  transEngine.velocity = velocity
  rotEngine.angularVelocity = angularVelocity
  rotEngine.zeroPoint += Vector3(0,0,1)*velocity*O.dt
  
  for i in range(int(IniDistanceWheelfromBoundary/binSize)+1,NBin):
    temBinClumpActive=[]

    if abs(rotEngine.zeroPoint[2]-(i*binSize-minZ))<=1e-7: #(i+1)
      for b in O.bodies:
        if isinstance(b.shape,Sphere):
          if b.isClumpMember:
            if (((i-int(IniDistanceWheelfromBoundary/binSize))*binSize-minZ)<b.state.pos[2]<=((i+int(DifferenceWCenterActiveEnd/binSize))*binSize-minZ)):
              temBinClumpActive.append(b.clumpId)

      temBinClumpActive=numpy.unique(numpy.array(temBinClumpActive))
      for r in temBinClumpActive:
        if O.bodies[r].isClump:
          O.bodies[r].dynamic=True
      number= i+int(DifferenceWCenterActiveEnd/binSize)
      ImportedBin= DataImportRigid(int(number))
      
  for i in range(int(IniDistanceWheelfromBoundary/binSize)+1,NBin):
    temBinClumpDeactive=[]
    temBinSphereDeactive=[]

    if abs(rotEngine.zeroPoint[2]-(i*binSize-minZ))<=1e-7: #0.75*(i)
      for b in O.bodies:
        if isinstance(b.shape,Sphere):
          if b.isClumpMember:
            if (b.state.pos[2]<=((i-int(IniDistanceWheelfromBoundary/binSize))*binSize-minZ)):
              temBinClumpDeactive.append(b.clumpId)
              temBinSphereDeactive.append(b.id)

      temBinSphereDeactive=numpy.unique(numpy.array(temBinSphereDeactive))
      temBinClumpDeactive=numpy.unique(numpy.array(temBinClumpDeactive))
      for p in temBinSphereDeactive:
        if O.bodies[p].isClumpMember:
          O.bodies[p].groupMask=4
      for r in temBinClumpDeactive:
        if O.bodies[r].isClump:
          O.bodies[r].dynamic=False
          O.bodies[r].state.blockedDOFs='xyzXYZ'

      for b in O.bodies:
        if isinstance(b.shape,Sphere):
          if (b.state.pos[2]<=((i-1-int(IniDistanceWheelfromBoundary/binSize))*binSize-minZ)):
            O.bodies.erase(b.id)

  collider.avoidSelfInteractionMask = 4
  O.engines[5].dead=False

  if rotEngine.zeroPoint[2]>=maxZ-1.05*Wheel_R:
    O.pause()

def eraseOffs():
  for b in O.bodies:
    if isinstance(b.shape,Sphere):
      if b.isClumpMember:
         if ((maxY+Wheel_R*1.5)<b.state.pos[1] or b.state.pos[1]<minY):
            O.bodies.erase(b.id)

def history():
  global Fx,Fy,Fz,Dx,Dy,Dz
  Fx=0
  Fy=0
  Fz=0
  for b in facets:
    Fx+=O.forces.f(b.id,sync=True)[0]
    Fy+=O.forces.f(b.id,sync=True)[1]
    Fz+=O.forces.f(b.id,sync=True)[2]
  Dx=rotEngine.zeroPoint[0]
  Dy=rotEngine.zeroPoint[1]
  Dz=rotEngine.zeroPoint[2] 
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
  ## SAVING NUMERICAL FORCE-DISPLACEMENT (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
  plot.saveDataTxt('../Terramechanics Delivarable Codes/Adaptive Reduced Order Modeling/GRC-3/Clump Grains/Facet Wheel/Full Width/Numerical Output/CSP-GRC3-B1-A-ClumpGrains-FullWidth-n057.txt')
## declare what is to plot. "None" is for separating y and y2 axis
plot.plots={'i':('Fx','Fy','Fz','Dx','Dy','Dz')}
plot.plot()
# Uncomment for non-interactive simulation
#O.run(2000000000000,True)

    
