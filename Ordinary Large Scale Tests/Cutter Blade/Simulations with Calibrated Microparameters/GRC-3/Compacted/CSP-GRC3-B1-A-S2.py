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
## This script models a typical free-surface cutter blade (CB) test on GRC-3. The bed size is considered 
## as a full width bed. The dimensions of the bed are l=7.68cm, w=7.68cm, and h=4.0cm.
## The simulation is conducted under non-reduced condition meaning the whole length of the bed
## is always present during the simulation and no adaptive reduced order modeling techniques are 
## applied.
## the bed has been prepared in its maximum dry density (minimum porosity) by uniaxial compression in 
## in the height direction and zero interparticle friction. The final porosity 
## after the bed stabilization is almost 0.315944326076811.
## GRC-3 particles are spheres with rolling/twisting resistance included (Hertz-Mindlin resistance model
## --see YADE manual for the activation of the rolling resistance).
## This CB test consists of two stages: 1- vertical penetration of the blade in the height direction 
## (-Y here) to the desired depth and with a penetration velocity, and 2- horizontal movement in the length direction 
## (+Z here) with a predefined travel velocity.Both penetration and cutting processes are conducted using the YADE 
## built-in TranslationEngine. 
## the simulation can be run in the passive mode (recommended for hpc) by keeping 
## or in the interactive mode (visualization mode) by remoning the O.run() command at the end of the
## script (in the latter case, the play bottun should be pressed). In either case the simulation is 
## paused when the movement of the blade center plus the half its thickness  reaches the 5/6 of the of the sample size
## in the +Z direction.   
## NOTE: 
## 1- this simulation is a non-mpi simulation which requires super computers to run. Development of 
## mpi version is the subject of future investigations.
## 2- the stability of the simulation is mainly controlled by the critiacl time step and the blade movement.
## velocity. In the case  of instability of grians (such as flying around), a smaller critical time step and/or 
## smaller travel velocity must be provided to stabilize the simulation. Additionally, the accuracy of force-
## displacement response (F-D) is highly dependent upon the velocity of the movement. Generraly, values close to
## experimental values result in more matching F-D responses. 
## 3- the material properties for the grains, are those determined from GA calibration of relevant triaxial tests.
## Note: the microparameter calibration of the used constituitive model (HM) does not guarantee 
## the extensibilty for large scale tests. The consistancy of such results is also highly dependent on the absence 
## of boundary effects necessitating the use of large DEM models (reduction in computation efficiency) and infinitesimal 
## blade movement velocity. 
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
deposFricDegree = 21.38 # CONTACT FRICTION
normalDamp=0.7 # NORMAL VISCOUS DAMPING
shearDamp=0.7 # SHEAR VISCOUS DAMPING
youngSoil=0.689e8 # CONTACT STIFFNESS FOR SOIL
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.330 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densSoil=1801.11858943976 # BULK DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
etaSoil=0.1 # HM ROLLING RESISTANCE PARAMETER
numDamp=0.4 # NUMERICAL DAMPING
IniDistanceBladefromBoundary = 0.0 # INITIAL DISTANCE OF THE RIGHT BOUNDARY OF THE CUTTER BLADE FROM THE LEFT BOUNDARY OF THE COMPACTED BED (-X). 
HeightBlade=0.07 # CUTTER BLADE HEIGHT (Y DIRECTION)
WidthBlade=0.02 # CUTTER BLADE WIDTH (Z DIRECTION) 
ThicknessBlade=0.0010 # CUTTER BLADE THICKNESS (X DIRECTION)
InitialPeneterationofBlade = 0.010 # PENETRATION DEPTH OF THE BLADE. NOTE: HERE, THE BLADE STARTS CUTTING FROM THE OUTSIDE OF THE BOUNDARY (-X)
# AND THIS PARAMETER IS ONLY USED FOR ITS INITIAL PLACEMENT.
InitialDistanceofBladefromTopSoil = 0 # INITIAL DISTANCE OF THE BOTTOM OF THE CUTTER BLADE FROM THE SURFACE OF THE COMPACTED BED. SINCE THE ASSEMBLY
# IS COMPLETELY SMOOTH, THIS PARAMETER CAN BE SET TO ZERO.
# NOTE: THE NUMERICAL STABILITY AND RESPNSE ACCURACY IS HIGHLY DEPENDENT UPON THE SIZE OF THE ACTIVE ZONE (THE LARGER, THE BETTER) AND THE BIN SIZE
# (THE SMALLER, THE BETTER).
HorizentalvelocityofBlade = 0.25e-2 #10 TIMES THE EXPERIMENTAL VELOCITY
WidthSample=0.05 # DEM SAMPLE WIDTH (Z)
HeightSample=0.03 # DEM SAMPLE HEIGHT (Y)
LengthSample=0.05 # DEM SAMPLE LENGTH (X)
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
## USE YOUR OWN LOCAL MACHINE ADDRESS
O.bodies.append(ymport.text("/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Simulations with Calibrated Microparameters/GRC-3/Compacted/Samples/CSP-GRC3-B1-A-S2.txt",shift=Vector3(0,0,0),material=SoilId,color=Vector3(0,0.5,0.7)))
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"maxX:",maxX,"minY:",minY,"maxY:",maxY,"minZ:",minZ,"maxZ:",maxZ)

mn,mx=Vector3(minX,minY,minZ),Vector3(maxX,0.03,maxZ)

## CREATE WALLS AROUND THE PACKING
walls=aabbWalls([mn,mx],thickness=0,material=ContainerId)
wallIds=O.bodies.append(walls)

O.bodies.append(utils.box(center=((maxX-minX)/2,(maxY+HeightBlade/2+InitialDistanceofBladefromTopSoil),minZ+IniDistanceBladefromBoundary),extents=(WidthBlade/2,HeightBlade/2,ThicknessBlade/2),material=ContainerId,fixed=True,color=(0.2,0.5,1),wire=False))
Blade=O.bodies[-1]
Blade.state.blockedDOFs='xXYZ'

############################
###   DEFINING ENGINES   ###
############################
gravity=(0,-9.81,0)
## RECORDING SIMULATION VTK 
## USE YOUR OWN LOCAL MACHINE ADDRESS
O.engines=[VTKRecorder(fileName='/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Simulations with Calibrated Microparameters/GRC-3/Compacted/VTK/CSP-GRC3-B1-A-S2/CSP-GRC3-B1-A-S2' ,recorders=['all'],iterPeriod=5000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()], label="collider"),
InteractionLoop(
      [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
      [Ip2_FrictMat_FrictMat_MindlinPhys(eta=etaSoil,betan=normalDamp,betas=shearDamp,label='ContactModel')],
      [Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
  ),
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
  Blade_PosZ=Blade.state.pos[2]+ThicknessBlade/2
  print ("Blade_PosZ:",Blade_PosZ)
  if Blade_PosZ>=5/6*maxZ:
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
  Dz=Blade.state.pos[2]+ThicknessBlade/2
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
  ## RECORDING NUMERICAL F-D HISTORY
  ## ## USE YOUR OWN LOCAL MACHINE ADDRESS  
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Simulations with Calibrated Microparameters/GRC-3/Compacted/Numerical Output/CSP-GRC3-B1-A-S2/CSP-GRC3-B1-A-S2.txt')
  
## DECLARE WHAT IS TO PLOT ON SCREEN. "None" IS FOR SEPARATING Y AND Y2 AXIS.
plot.plots={'i':('Fx','Fy','Fz','Dx','Dy','Dz')}
plot.plot()
# UNCOMMENT FOR NON-INTERACTIVE SIMULATION
#O.run(2000000000000,True)

