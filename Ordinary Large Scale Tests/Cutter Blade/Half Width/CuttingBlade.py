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
## THIS SCRIPT MODELS A TYPICAL FREE-SURFACE CUTTER BLADE TEST ON GRC-3. THE BED SIZE IS CONSIDERED 
## AS A HALF WIDTH BED. THE DIMENSIONS OF THE BED ARE L=20MM, W=5MM, AND H=5.5MM. 
## ALSO. THE SIMULATION IS CONDUCTED UNDER NON-REDUCED CONDITION MEANING THE WHOLE LENGTH OF THE BED
## IS ALWAYS PRESENT DURING THE SIMULATION AND NO ADAPTIVE REDUCED ORDER MODELING TECHNIQUES ARE 
## APPLIED.
## THE BED HAS BEEN PREPARED IN ITS MINIMUM DRY DENSITY (MAXIMUM POROSITY) BY GRAVITY DEPOSITION 
## OF GRAINS WITH THEIR ORIGINAL INTERPARTICLE FRICTION (HERE 40 DEGREES). THE FINAL POROSITY 
## AFTER THE BED TRIMMING IS ALMOST 0.57.
## GRC-3 PARTICLES ARE CLUMPS WITH NO ROLLING/TWISTING RESISTANCE INCLUDED (IT IS BELIEVED THAT THE
## CLUMPS COMPLEXITY ACCOUNT FOR SUCH RESISTANCES).
## THE BLADE WIDTH, HEIGHT AND THICKNESS ARE 2.5MM, 5MM, AND 0.15MM. THE BLADE PENETRATES AND MOVES ALONG
## THE RIGHTMOST BOUNDARY OF THE BED.
## THE TEST CONSISTS OF TWO STAGES: 1- VERTICAL PENETRATION OF THE BLADE IN Y DIRECTION TO 1.5MM BELOW 
## THE BED TOP SURFACE WITH A PENETRATION VELOCITY OF 0.07M/S, 2- HORIZONTAL MOVEMENT IN Z DIRECTION 
## WITH A TRAVEL VELOCITY OF 0.025M/S.  
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
deposFricDegree = 40 # INITIAL CONTACT FRICTION FOR SAMPLE PREPARATION
normalDamp=0.7 # NUMERICAL DAMPING
shearDamp=0.7
youngSoil=0.7e8# CONTACT STIFFNESS FOR SOIL
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.3 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densSoil=2650 # DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.4
IniDistanceBladefromBoundary = 0.001
HeightBlade=0.005
WidthBlade=0.0025
ThicknessBlade=0.00015
InitialPeneterationofBlade = 0.0015
InitialDistanceofBladefromTopSoil = 0
HorizentalvelocityofBlade = 0.07
PenetrationvelocityofBlade = 0.07

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
#'minX:', -2.927345865710862e-18, 'maxX:', 0.00506834, 'minY:', 0.0, 'maxY:', 0.0061737499999999996, 'minZ:', -1.3552527156068805e-20, 'maxZ:', 0.020000200000000003
## IMPORTING THE ALREADY COMPACTED BIN
ymport.textClumps("/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Half Width/5by6by20_n057_Relaxed_FreeSurface.txt",shift=Vector3(0,0,0),material=SoilId)
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"maxX:",maxX,"minY:",minY,"maxY:",maxY,"minZ:",minZ,"maxZ:",maxZ)

mn,mx=Vector3(minX,minY,minZ),Vector3(maxX,0.01,maxZ)

## CREATE WALLS AROUND THE PACKING
walls=aabbWalls([mn,mx],thickness=0,material=ContainerId)
wallIds=O.bodies.append(walls)
for b in O.bodies:
  if b.state.pos[1]>0.009:
    O.bodies.erase(b.id)

O.bodies.append(utils.box(center=(WidthBlade/2+minX,(maxY+HeightBlade/2+InitialDistanceofBladefromTopSoil),minZ+IniDistanceBladefromBoundary),extents=(WidthBlade/2,HeightBlade/2,ThicknessBlade/2),material=ContainerId,fixed=True,color=(0.2,0.5,1),wire=False))
Blade=O.bodies[-1]
Blade.state.blockedDOFs='xXYZ'

############################
###   DEFINING ENGINES   ###
############################
gravity=(0,-9.81,0)
O.engines=[VTKRecorder(fileName='home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Half Width/VTK/Output' ,recorders=['all'],iterPeriod=5000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()], label="collider"),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=normalDamp,betas=shearDamp,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
  ),
  ## WE WILL USE THE GLOBAL STIFFNESS OF EACH BODY TO DETERMINE AN OPTIMAL TIMESTEP (SEE HTTPS://YADE-DEM.ORG/W/IMAGES/1/1B/CHAREYRE&VILLARD2005_LICENSED.PDF)
  #GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.45),
  TranslationEngine(ids=[Blade.id],translationAxis=(0.,1.,0.),velocity=-1*PenetrationvelocityofBlade,label='BladePenetration',dead=True),
  TranslationEngine(ids=[Blade.id],translationAxis=(0.,0.,1.),velocity=HorizentalvelocityofBlade,label='BladeHorizentalMove',dead=True),  
  NewtonIntegrator(damping=numDamp,gravity=gravity),
  PyRunner(command='BladePenetration()',iterPeriod=1,label='checker'),
  PyRunner(iterPeriod=1000,command='eraseOffs()'),
  PyRunner(iterPeriod=1000,command='history()',label='recorder'),
  ]
O.dt = 1e-7

def BladePenetration():
  O.engines[4].dead=False
  Blade_Bottom_PosY = Blade.state.pos[1]-(HeightBlade/2)
  print ("Blade_Bottom_PosY:",Blade_Bottom_PosY)
  if Blade_Bottom_PosY<=(maxY-InitialPeneterationofBlade):
    O.engines[4].dead=True
    checker.command='BladeHorizentalMove()'

def BladeHorizentalMove():

  O.engines[4].dead=True
  O.engines[5].dead=False
  Blade_Center_PosZ=Blade.state.pos[2]
  print ("Blade_Center_PosZ:",Blade_Center_PosZ)
  if Blade_Center_PosZ>=5/6*maxZ:
    O.engines[5].dead=True 
    O.pause()

def eraseOffs():
  for b in O.bodies:
    if isinstance(b.shape,Sphere):
      if b.isClumpMember:
         if ((Blade.state.pos[1]+0.005)<b.state.pos[1] or b.state.pos[1]<minY):
            O.bodies.erase(b.id)

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
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Half Width/Output/CuttingBlade_OutputData_5by6by20_n057_Relaxed_FreeSurface.txt')
  
## declare what is to plot. "None" is for separating y and y2 axis
plot.plots={'i':('Fx','Fy','Fz','Dx','Dy','Dz')}
plot.plot()
  


#O.run(2000000000000,True)

