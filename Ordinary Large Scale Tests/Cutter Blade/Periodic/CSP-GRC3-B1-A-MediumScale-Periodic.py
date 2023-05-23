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
## THIS SCRIPT MODELS A TYPICAL PERIODIC (LENGTH ONLY) FREE-SURFACE CUTTER BLADE TEST ON GRC-3. THE BED SIZE IS CONSIDERED 
## AS A FULL WIDTH BED. THE DIMENSIONS OF THE BED ARE L=12MM, W=9.0MM, AND H=2.5MM. 
## ALSO. THE SIMULATION IS CONDUCTED UNDER NON-REDUCED CONDITION MEANING THE WHOLE LENGTH OF THE BED
## IS ALWAYS PRESENT DURING THE SIMULATION AND NO ADAPTIVE REDUCED ORDER MODELING TECHNIQUES ARE 
## APPLIED.
## THE BED HAS BEEN PREPARED IN ITS MAXIMUM DRY DENSITY (MINIMUM POROSITY) BY UNIAXIAL COMPRESSION IN 
## IN THE HEIGHT DIRECTION AND ZERO INTERPARTICLE FRICTION. THE FINAL POROSITY 
## AFTER THE BED STABILIZATION IS ALMOST 0.315944326076811.
## GRC-3 PARTICLES ARE SPHERES WITH ROLLING/TWISTING RESISTANCE INCLUDED (Hertz-Mindlin with rolling resistance-- see YADE manual for the activation of the rolling resistance)
## THE BLADE WIDTH, HEIGHT AND THICKNESS ARE 3.5MM, 7MM, AND 0.20MM. THE BLADE PENETRATES AND MOVES AT 
## THE MIDDLE OF THE BED WIDTH. 
## A TYPICAL CB TEST CONSISTS OF TWO STAGES: 1- VERTICAL PENETRATION OF THE BLADE IN Y DIRECTION TO A CERTAIN DEPTH BELOW 
## THE BED TOP SURFACE WITH A PENETRATION VELOCITY, 2- HORIZONTAL MOVEMENT IN THE LENGTH DIRECTION 
## WITH A TRAVEL VELOCITY. HOWEVER IN THIS SIMULATION THE PENETRATION STAGE HAS BBEN SKIPPED AND THE SIMULATION STARTS
## WITH CUTTING FROM THE OUT OF THE BOUNDARY (-X).   
## THE SIMULATION CAN BE RUN IN THE PASSIVE MODE (RECOMMENDED FOR HPC) BY KEEPING 
## OR IN THE INTERACTIVE MODE (VISUALIZATION MODE) BY REMONING THE O.run() COMMAND AT THE END OF THE
## SCRIPT (IN THE LATTER CASE, THE PLAY BOTTUN SHOULD BE PRESSED). IN EITHER CASE THE SIMULATION IS 
## PAUSED WHEN THE MOVEMENT OF THE BLADE CENTER REACHES THE 5/6 OF THE BED LENGTH (20MM).   
## NOTE: 
## 1- THIS SIMULATION IS A NON-MPI SIMULATION WHICH REQUIRES SUPER COMPUTERS TO RUN. DEVELOPMENT OF 
## MPI VERSION OF THE SIMULATION IS UNDER PROGRESS.
## 2- THE STABILITY OF THE SIMULATION IS MAINLY CONTROLLED BY THE CRITIACL TIME STEP (LINE 130).
## IN THE CASE  OF INSTABILITY OF GRIANS (SUCH AS FLYING AROUND), A SMALLER CRITICAL TIME STEP SHOULD 
## BE PROVIDED TO STABILIZE THE SIMULATION.
## 3- THE MATERIAL PROPERTIES FOR THE GRAINS, CONTAINER, AND BLADE HAVE BEEN CHOSEN ARBITRARILY AND CAN BE
## CHANGED FOR CALIBRATION PURPOSES. 
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
deposFricDegree = 21.42 # INITIAL CONTACT FRICTION FOR SAMPLE PREPARATION
normalDamp=0.7 # NUMERICAL DAMPING
shearDamp=0.7
youngSoil=0.75e8 # CONTACT STIFFNESS FOR SOIL
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.44 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densSoil=1801.11858943976 # BULK DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.4
IniDistanceBladefromBoundary = 0.002
HeightBlade=0.007
WidthBlade=0.0035
ThicknessBlade=0.00020
InitialPeneterationofBlade = 0.0010
InitialDistanceofBladefromTopSoil = 0
HorizentalvelocityofBlade = 0.5
PenetrationvelocityofBlade = 0.04

O.periodic=True
length=0.013
height=0.015
width=0.00768
thickness=0.01

O.cell.hSize=Matrix3(length, 0, 0,
      0 ,height, 0,
      0, 0, width)


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
## IMPORTING THE ALREADY COMPACTED BIN
O.bodies.append(ymport.text("/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Periodic/RelaxedBin_CSP-GRC3-B1-A.txt",shift=Vector3(0.001,0.003,0),material=SoilId,color=Vector3(0,0.5,0.7)))
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"maxX:",maxX,"minY:",minY,"maxY:",maxY,"minZ:",minZ,"maxZ:",maxZ)

mn,mx=Vector3(0.001,0.003,minZ),Vector3(maxX,0.01,maxZ)

## CREATE WALLS AROUND THE PACKING
walls=aabbWalls([mn,mx],thickness=0,material=ContainerId,oversizeFactor=1.0)
wallIds=O.bodies.append(walls)
# for b in O.bodies:
#   if b.state.pos[1]>0.009:
#     O.bodies.erase(b.id)

O.bodies.append(utils.box(center=((minX-ThicknessBlade/2),(maxY+HeightBlade/2-InitialPeneterationofBlade),(maxZ-minZ)/2),extents=(ThicknessBlade/2,HeightBlade/2,WidthBlade/2),material=ContainerId,fixed=True,color=(0.8,0.1,0.2),wire=False))
Blade=O.bodies[-1]
Blade.state.blockedDOFs='yzXYZ'

############################
###   DEFINING ENGINES   ###
############################
gravity=(0,-9.81,0)
O.engines=[VTKRecorder(fileName='home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Periodic/VTK/VTKPeriodic' ,recorders=['all'],iterPeriod=10000),

ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()],allowBiggerThanPeriod=True,label="collider"),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=normalDamp,betas=shearDamp,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
  ),
  ## WE WILL USE THE GLOBAL STIFFNESS OF EACH BODY TO DETERMINE AN OPTIMAL TIMESTEP (SEE HTTPS://YADE-DEM.ORG/W/IMAGES/1/1B/CHAREYRE&VILLARD2005_LICENSED.PDF)
  TranslationEngine(ids=[Blade.id],translationAxis=(1.,0.,0.),velocity=HorizentalvelocityofBlade,label='BladeHorizentalMove',dead=True),  
  NewtonIntegrator(damping=numDamp,gravity=gravity),
  PyRunner(command='BladeHorizentalMove()',iterPeriod=1,label='checker'),
  #PyRunner(iterPeriod=1000,command='eraseOffs()'),
  PyRunner(iterPeriod=1000,command='history()',label='recorder'),
  ]
O.dt = 1.2e-7

def BladeHorizentalMove():
  #O.engines[4].dead=True
  O.engines[4].dead=False
  Blade_Center_PosX=Blade.state.pos[0]
  print ("Blade_Center_PosX:",Blade_Center_PosX)
  if Blade_Center_PosX>=5/6*maxX:
    O.engines[4].dead=True 
    O.pause()

# def eraseOffs():
#   for b in O.bodies:
#     if isinstance(b.shape,Sphere):
#       if b.isClumpMember:
#          if ((Blade.state.pos[1]+0.005)<b.state.pos[1] or b.state.pos[1]<minY):
#             O.bodies.erase(b.id)

def history():
  global Fx,Fy,Fz,Dx,Dy,Dz
  Fx=O.forces.f(Blade.id,sync=True)[0]
  Fy=O.forces.f(Blade.id,sync=True)[1]
  Fz=O.forces.f(Blade.id,sync=True)[2]
  Dx=Blade.state.pos[0]
  Dy=Blade.state.pos[1]
  Dz=Blade.state.pos[2]
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
  ## In that case we can still save the data to a text file at the the end of the simulation, with:
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Periodic/Output/CSP-GRC3-B1-A-Periodic.txt')
  
## declare what is to plot. "None" is for separating y and y2 axis
plot.plots={'Dx':('Fx','Fy','Fz')}
plot.plot()
  


#O.run(2000000000000,True)

