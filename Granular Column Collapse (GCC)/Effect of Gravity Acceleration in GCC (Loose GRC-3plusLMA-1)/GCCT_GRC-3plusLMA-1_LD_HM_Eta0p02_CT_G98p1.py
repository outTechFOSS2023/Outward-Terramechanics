# -*- coding: utf-8 -*-
#*********************************************************************************************************
#*********************************************************************************************************
#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, #distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY #CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#*********************************************************************************************************
#*********************************************************************************************************                                                                     
#  THIS PROGRAM IS A FREE SOFTWARE; IT IS LICENCED UNDER THE TERMS OF THE 
#  GNU GENERAL PUBLIC LICENCE V2 OR LATER
#*************************************************************************/
## THIS SCRIPT SIMULATES THE COLLAPSE OF A GRANULAR COLUMN OF GRC-3+LMA-1 UNDER THE 
## EFFECT OF GRAVITY. THE COLLAPSE OCCURES IN +X DIRECTION ONLY BY ROTATING
## THE RIGHT WALL FROM A PIVOT POINT ON TOP OF THE CONTAINER.
## THIN SLICES (IN Z DIRECTION) CAN BE COLLAPSED WITH THIS SCRIPT 
## PSD: HALF PSD
## SAMPLE LENGTH AND WIDTH: 2.524mm
## SAMPLE HEIGHT: 10.44mm
## COLLAPSE SYSTEM: ONE SIDE 
## POROSITY: LOOSE (n=0.558)
## NUMERICAL TO EXPERIMENTAL SCALE: 1:15
## PERIODIC VS NON-PERIODIC: NON-PERIODIC
## GRAIN SIZE LIMIT: 90-500mm
## GRAVITY=-98.1m/s
#*************************************************************************/
############################################
#### SEPARATING SPHERES AND CLUMPS      ####
#### FOR ROLLING RESISTANCE OF SPHERES  ####
############################################
AllSpheres = open('/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Samples/GCCT_01NU15EX_GRC-3plusLMA-1_LL90UL500_HPSD_TCO_LD_2p524by10p440by2p524_NP.txt',"r")
SpheresData = numpy.loadtxt(AllSpheres, usecols=range(5))
a = 1
NumberofGrains = []
for b in range(1,len(SpheresData)):
  if SpheresData[b,4] == SpheresData[b-1,4]:
    a=a+1
    if (b==(len(SpheresData)-1)):
      NumberofGrains.append(numpy.array([a]))
  else:
    NumberofGrains.append(numpy.array([a]))
    a = 1
NumberofGrains=numpy.array(NumberofGrains).reshape(-1)

k = 0
Sphere = []
Clumps = []

for i in range(0,len(NumberofGrains)):
  k = k + NumberofGrains[i]
  l=(k-NumberofGrains[i])
  if NumberofGrains[i]==1:
    Sphere.append(numpy.array([SpheresData[l,0],SpheresData[l,1],SpheresData[l,2],SpheresData[l,3],'{:d}'.format(int(SpheresData[l,4]))]))

  if NumberofGrains[i]!=1:
    for p1 in range (int(l),int(k)):
      Clumps.append(numpy.array([SpheresData[p1,0],SpheresData[p1,1],SpheresData[p1,2],SpheresData[p1,3],'{:d}'.format(int(SpheresData[p1,4]))]))

Clumps = numpy.array(Clumps).reshape(-1,5)
Sphere = numpy.array(Sphere).reshape(-1,5)

with open(('/home/ngoudarzi/Desktop/Effect of Gravity Acceleration in GCC /Samples/GCCT_01NU15EX_GRC-3plusLMA-1_LL90UL500_HPSD_TCO_LD_2p524by10p440by2p524_NP_TemporarySphere.txt'),"w") as f:   
  f.write('#format x_y_z_r_clumpId')
  for S in range(0,len(Sphere)):
    f.write('\n')
    f.write("\n".join(" ".join(map(str, S)) for S in Sphere[S]\
          .reshape(-1,5)))
    
with open(('/home/ngoudarzi/Desktop/Effect of Gravity Acceleration in GCC /Samples/GCCT_01NU15EX_GRC-3plusLMA-1_LL90UL500_HPSD_TCO_LD_2p524by10p440by2p524_NP_TemporaryClumps.txt'),"w") as f:   
  f.write('#format x_y_z_r_clumpId')
  for S in range(0,len(Clumps)):
    f.write('\n')
    f.write("\n".join(" ".join(map(str, S)) for S in Clumps[S]\
          .reshape(-1,5)))

from yade.pack import *
from yade import utils
from yade import export
from yade import ymport
from yade import plot
import time
import numpy
import os
import math

############################################
##########   DEFINING VARIABLES    #########
############################################
##** REGOLITH MATERIAL
RegolithFricDegree=30 # INTERPARTICLE FRICTION DURING DEPOSITION
RegolithYoung=1.0e8 # CONTACT STIFFNESS FOR REGOLITH
RegolithPoisson=0.3 # POISSION'S RATIO FOR REGOLITH (CLUMP)
RegolithDensity=2650 # DENSITY FOR REGOLITH
RegolithEta=0.2 # SHAPE PARAMETER IN HM ROLLING RESISTANCE MODEL (0<BETA<0.5)

##** BOUNDARIES MATERIAL
FrictlessBoundFricDegree=0.0 # FRICTIONLESS BOUNDARIES IN THE MODEL
LowFrictFrictionalBoundFricDegree=2.5 # LOW FRICTION FRICTIONAL BOUNDARIES IN THE MODEL
HighFrictFrictionalBoundFricDegree=5.0 # HIGH FRICTION FRICTIONAL BOUNDARIES IN THE MODEL
BoundYoung=210e9 # CONTACT STIFFNESS FOR ALL BOUNDARIES
BoundPoisson=0.2 # POISSION'S RATIO FOR ALL BOUNDARIES
BoundDensity=7850 # DENSITY FOR ALL BOUNDARIES

##** ANALYSIS CONSTANTS
Thick=0.0
NormalViscDamp=0.3 # NORMAL VISCOUS DAMPING (EXP: MAKES UNSTABLE AND IS NOT USED IN COLLAPSE MODELING. IF NEEDED, SHOULD BE ENTERED IN iPhys)
ShearViscDamp=0.3  # SHEAR VISCOUS DAMPING  (EXP: MAKES UNSTABLE AND IS NOT USED IN COLLAPSE MODELING. IF NEEDED, SHOULD BE ENTERED IN iPhys)
NumDamp=0.4 # NUMERICAL DAMPING
GravAcc=(0,-98.1,0) # GRAVITY ACCELERATION 
#********************************************
MN_Boundary,MX_Boundary=Vector3(0,0,0),Vector3(0.025,15.0e-3,2.524e-3) # CORNERS OF THE COLLAPSE RUNWAY (NO WALL ROTATION IS MODELED HERE)
MN_RotWall,Mx_RotWall=Vector3(0,0,0),Vector3(2.524e-3,10.440e-3,2.524e-3) # ROTATING WALL
MN_Benchmark,Mx_MN_Benchmark=Vector3(0,0,0),Vector3(4.5e-3,10.440e-3,2.524e-3) # ROTATING WALL

############################################
###############   MATERIALS   ##############
############################################
O.materials.append(FrictMat(young=RegolithYoung,poisson=RegolithPoisson,frictionAngle=radians(RegolithFricDegree),density=RegolithDensity,label='Regolith'))
O.materials.append(FrictMat(young=BoundYoung,poisson=BoundPoisson,frictionAngle=radians(LowFrictFrictionalBoundFricDegree),density=BoundDensity,label='LowFrictionalFrictionalBoundaries'))
O.materials.append(FrictMat(young=BoundYoung,poisson=BoundPoisson,frictionAngle=radians(HighFrictFrictionalBoundFricDegree),density=BoundDensity,label='HighFrictionalFrictionalBoundaries'))
O.materials.append(FrictMat(young=BoundYoung,poisson=BoundPoisson,frictionAngle=radians(FrictlessBoundFricDegree),density=BoundDensity,label='FrictionlessBoundaries'))

############################################
########## ROTATING WALL GEOMETRY  #########
############################################
RotWall=aabbWalls([MN_RotWall,Mx_RotWall],thickness=Thick,material='FrictionlessBoundaries',oversizeFactor=1.0,color=(1,1,1),wire=False)
id_RotWall=O.bodies.append(RotWall)
O.bodies.erase(0)
O.bodies.erase(3)
O.bodies.erase(2)
O.bodies.erase(4)
O.bodies.erase(5)

############################################
############# HEIGHT BENCHMARK  ############
############################################
HeightBenchmark=aabbWalls([MN_Benchmark,Mx_MN_Benchmark],thickness=Thick,material='FrictionlessBoundaries',oversizeFactor=1.0,color=(1,1,1),wire=False)
id_HeightBenchmark=O.bodies.append(HeightBenchmark)
O.bodies.erase(6)
O.bodies.erase(7)
O.bodies.erase(8)
O.bodies.erase(10)
O.bodies.erase(11)
  
############################################
##########     RUNWAY GEOMETRY     #########
############################################
Floor=aabbWalls([MN_Boundary,MX_Boundary],thickness=Thick,material='HighFrictionalFrictionalBoundaries',oversizeFactor=1.0,color=(1,1,1),wire=False)
id_Floor=O.bodies.append(Floor)
O.bodies.erase(12)
O.bodies.erase(13)
O.bodies.erase(16)
O.bodies.erase(17)

Backwall=aabbWalls([MN_Boundary,MX_Boundary],thickness=Thick,material='LowFrictionalFrictionalBoundaries',oversizeFactor=1.0,color=(1,1,1),wire=False)
id_Backwall=O.bodies.append(Backwall)
O.bodies.erase(19)
O.bodies.erase(20)
O.bodies.erase(21)
O.bodies.erase(22)
O.bodies.erase(23)

Sidewalls=aabbWalls([MN_Boundary,MX_Boundary],thickness=Thick,material='FrictionlessBoundaries',oversizeFactor=1.0,color=(1,1,1),wire=False)
id_Sidewalls=O.bodies.append(Sidewalls)
O.bodies.erase(24)
O.bodies.erase(26)
O.bodies.erase(27)

############################################
##########     REGOLITH COLUMN     #########
############################################
ymport.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Samples/GCCT_01NU15EX_GRC-3plusLMA-1_LL90UL500_HPSD_TCO_LD_2p524by10p440by2p524_NP_TemporarySphere.txt", material='Regolith',color=Vector3(0.7,0.7,0.9))
ymport.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Samples/GCCT_01NU15EX_GRC-3plusLMA-1_LL90UL500_HPSD_TCO_LD_2p524by10p440by2p524_NP_TemporaryClumps.txt", material='Regolith',color=Vector3(0.4,0.1,0.0))
MinX_Import=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxX_Import=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MinY_Import=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxY_Import=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MinZ_Import=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxZ_Import=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
SampleLength=(MaxX_Import-MinX_Import)
SampleHeight=(MaxY_Import-MinY_Import)
SampleWidth=(MaxZ_Import-MinZ_Import)
print('MinX_Import:',MinX_Import,'MaxX_Import:',MaxX_Import,'MinY_Import:',MinY_Import,'MaxY_Import:',MaxY_Import,'MinZ_Import:',MinZ_Import,'MaxZ_Import:',MaxZ_Import)
print('SampleLength:',SampleLength,'SampleHeight:',SampleHeight,'SampleWidth:',SampleWidth)
Porosity=utils.porosity()
V_S=getSpheresVolume()
print('Soild Volume Clculated by YADE:',V_S,'Porosity:',Porosity)
os.remove("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Samples/GCCT_01NU15EX_GRC-3plusLMA-1_LL90UL500_HPSD_TCO_LD_2p524by10p440by2p524_NP_TemporarySphere.txt")
os.remove("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Samples/GCCT_01NU15EX_GRC-3plusLMA-1_LL90UL500_HPSD_TCO_LD_2p524by10p440by2p524_NP_TemporaryClumps.txt")

#########################################################
######################   ENGINES   ######################
#########################################################
Newton=NewtonIntegrator(damping=NumDamp,gravity=GravAcc)
O.engines=[VTKRecorder(fileName='/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/VTK/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_CT_G98p1/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_G98p1_VTK' ,recorders=['all'],iterPeriod=50000), 
  ForceResetter(),
  InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
  InteractionLoop(
    [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
    [Ip2_FrictMat_FrictMat_MindlinPhys(eta=RegolithEta)],
    [Law2_ScGeom_MindlinPhys_Mindlin(includeMoment=True)]),
  Newton,
  #GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.95),
  PyRunner(command='RotateTheRightWall()',iterPeriod=300,label='Checker'),
  RotationEngine(angularVelocity=100,rotationAxis=(0,0,1),rotateAroundZero=True,zeroPoint=(2.524e-3,10.44e-3,0),ids=id_RotWall,label='Rot',dead=True), # ANGULAR VELOCITY IS A CALIBRATION PARAMETER
]
O.dt=1.0e-6
t=0
def RotateTheRightWall():
  global t
  Rot.dead=False
  if O.bodies[1].state.pos[1]>=10.440e-3:#THIS STOPS THE ROTATION OF THE MOVING WALL WELL ABOVE THE COLUMN
    Rot.angularVelocity=0.0
  O.trackEnergy=True
  Ek=utils.kineticEnergy()
  MinX_collapse=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
  MaxX_collapse=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
  MinY_collapse=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
  MaxY_collapse=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
  MaxZ_collapse=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
  UNB=unbalancedForce()
  O.trackEnergy=True
  Ek=utils.kineticEnergy()
  print ("Kinetic Energy: ",Ek, "Unbalanced Force: ",UNB, "MaxZ During the Collapse: ",MaxZ_collapse, "MinX During the Collapse: ",MinX_collapse, "MaxX During the Collapse: ",MaxX_collapse, "MinY During the Collapse: ",MinY_collapse, "MaxY During the Collapse: ",MaxY_collapse)
  if O.iter==30000+(t*30000):
    export.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Numerical Output/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_CT_G98p1/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_G98p1.txt") 
    O.save("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Numerical Output/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_CT_G98p1/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_G98p1.yade.bz2")
    t=t+1
    print ("T: ",t)
  if O.iter>2500000: #AN ADITIONAL CRITERION IS REQUIRED TO ENSURE STABILITY OF THE COLLAPSE
    export.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Numerical Output/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_CT_G98p1/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_G98p1.txt") 
    O.save("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Effect of Gravity Acceleration in GCC (Loose GRC-3plusLMA-1)/Numerical Output/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_CT_G98p1/GCCT_GRC-3plusLMA-1_LD_HM_Eta0p02_G98p1.yade.bz2")
    O.pause()

