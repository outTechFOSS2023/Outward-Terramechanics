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
## THIS SCRIPT MODELS A TYPICAL FREE-SURFACE SMOOTH WHEEL TEST ON GRC-3. THE BED SIZE IS CONSIDERED 
## AS A HALF WIDTH BED. THE DIMENSIONS OF THE BED ARE L=20MM, W=5MM, AND H=5.5MM. 
## ALSO. THE SIMULATION IS CONDUCTED UNDER NON-REDUCED CONDITION MEANING THE WHOLE LENGTH OF THE BED
## IS ALWAYS PRESENT DURING THE SIMULATION AND NO ADAPTIVE REDUCED ORDER MODELING TECHNIQUES ARE 
## APPLIED.
## THE BED HAS BEEN PREPARED IN ITS MINIMUM DRY DENSITY (MAXIMUM POROSITY) BY GRAVITY DEPOSITION 
## OF GRAINS WITH THEIR ORIGINAL INTERPARTICLE FRICTION (HERE 40 DEGREES). THE FINAL POROSITY 
## AFTER THE BED TRIMMING IS ALMOST 0.57.
## GRC-3 PARTICLES ARE CLUMPS WITH NO ROLLING/TWISTING RESISTANCE INCLUDED (IT IS BELIEVED THAT THE
## CLUMPS COMPLEXITY ACCOUNT FOR SUCH RESISTANCES).
## THE WHEEL WIDTH, AND RADIUS ARE 1.261MM AND 2.5MM. IT HAS BEEN MODELED AS A FACET CYLINDER OF 
## 200 FACETS. THE WHEEL PENETRATES AND MOVES ALONG THE RIGHTMOST BOUNDARY OF THE BED. 
## THE TEST CONSISTS OF TWO STAGES: 1- VERTICAL PENETRATION OF THE WHEEL IN Y DIRECTION TO 1.0MM BELOW 
## THE BED TOP SURFACE WITH A PENETRATION VELOCITY OF 0.07M/S, 2- HORIZONTAL MOVEMENT IN Z DIRECTION 
## WITH A TRAVEL HORIZONTAL VELOCITY OF 0.01M/S AND A SLIDING RATIO OF 0.3 FROM WHICH THE WHEEL ANGULAR
## VELOCITY IS CALCULATED. THE WHHEL MASS IS CONSIDERED TO BE 0.5 KG AND ADDITIONAL NORMAL FORCE OF 19.6
## N IS APPLIED TO IT DURING THE SIMULATION.  
## THE SIMULATION CAN BE RUN IN THE PASSIVE MODE (RECOMMENDED FOR HPC) BY KEEPING 
## OR IN THE INTERACTIVE MODE (VISUALIZATION MODE) BY REMONING THE O.run() COMMAND AT THE END OF THE
## SCRIPT (IN THE LATTER CASE, THE PLAY BOTTUN SHOULD BE PRESSED). IN EITHER CASE THE SIMULATION IS 
## PAUSED WHEN THE MOVEMENT OF THE WHEEL CENTER REACHES THE LENGTH OF THE BED-1.05 WHEEL RADIUS
## NOTE: 
## 1- THIS SIMULATION IS A NON-MPI SIMULATION WHICH REQUIRES SUPER COMPUTERS TO RUN. DEVELOPMENT OF 
## MPI VERSION OF THE SIMULATION IS UNDER PROGRESS.
## 2- THE STABILITY OF THE SIMULATION IS MAINLY CONTROLLED BY THE CRITIACL TIME STEP (LINE 148).
## IN THE CASE  OF INSTABILITY OF GRIANS (SUCH AS FLYING AROUND), A SMALLER CRITICAL TIME STEP SHOULD 
## BE PROVIDED TO STABILIZE THE SIMULATION.
## 3- THE MATERIAL PROPERTIES FOR THE GRAINS, CONTAINER, AND WHEEL HAVE BEEN CHOSEN ARBITRARILY AND CAN BE
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
############################################
###           INPUT PARAMETERS           ###
############################################
## Wheel PROPERTIES
widthWheel=0.00125
Wheel_R=0.0025
Wheel_AdditionalNormalForce=-19.6

##WHELL MOVEMENT
slidingratio=0.3
velocity=0.1 #m0.01 m/s IS THE MOST STABLE VELOCITY
angularVelocity=velocity/((1-slidingratio)*Wheel_R) #RAD/S
gravityAcc=-9.81

deposFricDegree = 28.5 # INITIAL CONTACT FRICTION FOR SAMPLE PREPARATION
normalDamp=0.7 # NUMERICAL DAMPING
shearDamp=0.7
youngSoil=0.7e8# CONTACT STIFFNESS FOR SOIL
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.3 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densSoil=2650 # DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.4
binSize=0.4*Wheel_R
activeDistance=0.01
iniDistanceWheelfromBoundary=0.004
differenceWCenterActiveEnd=activeDistance - iniDistanceWheelfromBoundary
initialPeneterationofWheel=0.001
iniWheelDisfromMaxY=0.0001

############################################
###   DEFiniNG VARIABLES AND MATERIALS   ###
############################################
SoilId=FrictMat(young=youngSoil,poisson=poissonSoil,frictionAngle=radians(deposFricDegree),density=densSoil)
O.materials.append(SoilId)


ContainerId=FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=radians(0),density=densContainer)
O.materials.append(ContainerId)

###################################
#####   CREATING GEOMETRIES   #####
###################################

## IMPORTING THE ALREADY COMPACTED BIN
ymport.textClumps("/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/SWI/5by6by20_n057_Relaxed_FreeSurface.txt",shift=Vector3(0,0,0),material=SoilId)
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

XcenterofWheel = widthWheel/2+minX
YcenterofWheel = Wheel_R+maxY+iniWheelDisfromMaxY
ZcenterofWheel = iniDistanceWheelfromBoundary+minZ

Wheel1 = geom.facetCylinder((XcenterofWheel,YcenterofWheel,ZcenterofWheel),radius=Wheel_R,height=widthWheel,orientation=Quaternion((0,1,0),-pi/2),wallMask=7,segmentsNumber=200,material=ContainerId,wire=True,color=Vector3(255,255,255))
O.bodies.append(Wheel1)
idsr = [w.id for w in Wheel1]
facets = [b for b in O.bodies if isinstance(b.shape,Facet)] # list of facets in simulation
for b in facets:
  O.forces.setPermF(b.id,(0,Wheel_AdditionalNormalForce/800,0))

############################
###   DEFiniNG ENGINES   ###
############################
O.engines=[VTKRecorder(fileName='/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/SWI/VTK/Output' ,recorders=['all'],iterPeriod=5000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Facet_Aabb(),Bo1_Subdomain_Aabb()], label="collider"),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=normalDamp,betas=shearDamp,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
  ),
  ## WE WILL USE THE GLOBAL STIFFNESS OF EACH BODY TO DETERMINE AN OPTIMAL TIMESTEP (SEE HTTPS://YADE-DEM.ORG/W/IMAGES/1/1B/CHAREYRE&VILLARD2005_LICENSED.PDF)
  #GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.25),
  CombinedKinematicEngine(ids=idsr,label='combEngine') + TranslationEngine(translationAxis=(0,-1,0),velocity=velocity) +\
  RotationEngine(rotationAxis=(1,0,0), angularVelocity=0, rotateAroundZero=True, zeroPoint=(XcenterofWheel,YcenterofWheel,ZcenterofWheel)),
  NewtonIntegrator(damping=numDamp,gravity=(0,gravityAcc,0)),
  PyRunner(iterPeriod=1, command="WheelPenetration()" ,label='checker'),
  PyRunner(iterPeriod=100,command='history()',label='recorder'),
  PyRunner(iterPeriod=1000,command='eraseOffs()'),
  ]
O.dt = 1e-7

# get TranslationEngine and RotationEngine from CombinedKinematicEngine
transEngine, rotEngine = combEngine.comb[0], combEngine.comb[1]

def WheelPenetration():
  transEngine.velocity = 7*velocity
  rotEngine.zeroPoint += Vector3(0,-1,0)*transEngine.velocity*O.dt
  print ("rotEngine.zeroPoint[1]-Wheel_R):",(rotEngine.zeroPoint[1]-Wheel_R),"Iteration:",O.iter)
  if (rotEngine.zeroPoint[1]-Wheel_R)<=(maxY-initialPeneterationofWheel):
    checker.command='WheelRolling()'

def WheelRolling():
  transEngine.translationAxis=(0,0,1)  
  transEngine.velocity = velocity
  rotEngine.angularVelocity = angularVelocity
  rotEngine.zeroPoint += Vector3(0,0,1)*velocity*O.dt
  print ("rotEngine.zeroPoint[2]):",(rotEngine.zeroPoint[2]),"Iteration:",O.iter)
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
  ## In that case we can still save the data to a text file at the the end of the simulation, with:
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/SWI/Output/Wheel_OutputData_5by6by20_n057_Relaxed_FreeSurface_Vel0.1.txt')
  
## declare what is to plot. "None" is for separating y and y2 axis
plot.plots={'i':('Fx','Fy','Fz','Dx','Dy','Dz')}
plot.plot()
#O.run(2000000000000,True)

    
