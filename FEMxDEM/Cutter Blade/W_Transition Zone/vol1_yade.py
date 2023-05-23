#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

######################################################################
# This Python script is the YADE (DEM) part of a FEMxDEM simulation of a cutter blade (CB) test.
# The FEM part of the simulation is defined in a diferent Python script entitled vol1_oofem.py.
# There are also some other dependencies which link the DEM and FEM for such adjustments of  as output request and number of iterations.
# NOTE: ALL EXTERNAL FILES (DEM GEOMETRY), AND OUTPUT FOLDERES ARE PROVIDED. THE CORRECT ADDRESSES
# MUST BE PROVIDED BY THE USER FROM THEIR OWN LOCAL MACHINES. 
# The DEM model has the following chracteristics:
# 1- No microparameter calibrations were conducted and the model has been setup to show the runability of  FEMxDEM.
# 2- The DEM sample is an assembly of a mixture of clumps and spherical perticles used to respresent  a very loose state of  GRC-3+LMA-1 (n=0.57)
# 3- DEM assembly size is 5mm(X-direction)*2.5mm(Y-direction)x20mm(Z-Direction) supported by a 3mm thick underlying FEM mesh. The first bottom 
# 1mm of the DEM geometry is composed of purely spherical grains as a transition zone between the FEM and clumps (upper 1.5mm). However, the use of
# this transition zone is not necessary and clumps can come into direct contact with FEM. 
# 4- Surface coupling (SC) FEM*DEM has been implemented. SC has been shown to run faster than volumetic coupling (VC). 
# 5- Wall boundaries and the cutter blade has been generated using YADE built-in geometry generation functions for facets and therefore are not deformable. 
# However, it is possible to extend the FEM to all external objects.
# 6-  A simple FrictMat constituitive law has been used for both regolith, boundaries, and the cutter blade. 
# 7-  The whole process starts with the penetration of the cutter blade to 1mm into the regolith and then the
# cutting of the regolith in the +Z direction.
# 8- The penetration and cutting processes are simulated using YADE CombinedKinematicEngine. 
######################################################################
import math
from math import sqrt
from math import pi
from libyade import yade
from yade import *
from yade import pack, plot 
import numpy
from pyquaternion import Quaternion
from yade import pack, export, Vector3
from yade import ymport
import sys
import os
import os.path
import IPython
import __builtin__
############################################
###           INPUT PARAMETERS           ###
############################################
#### UNITS
## LENGTH: m
## VELOCITY: m/s
## TIME: s
## PRESSURE: pa
## DENSITY: kg/m3
deposFricDegree = 28.5 #  CONTACT FRICTION
IniDistanceBladefromBoundary = 0.001 # INITIAL DISTANCE OF THE BLADE FROM THE LEFT BOUNDARY DURING PENETRATION
HeightBlade=0.005 # HEIGHT OF THE CUTTER BLADE (m) (Y-DIRECTION)
WidthBlade=0.0025 # WIDTH OF THE CUTTER BLADE (m) (X-DIRECTION)
ThicknessBlade=0.00015 # THICKNESS OF THE CUTTER BLADE (m) (Z-DIRECTION)
InitialPeneterationofBlade = 0.001 # PENETRATION DEPTH OF THE BLADE BEFORE THE START OF THE CUTTING PROCESS
InitialDistanceofBladefromTopSoil = 0 # INITIAL OFFSET OF THE BOTTOM OF THE BLADE FROM THE TOP OF THE SAMPLE
velocity = 0.050 # CUTTING VELOCITY. IN ORDER TO ENSURE THE NUMERICAL STABILITY OF THE MODEL, THE VELOCITY MUST BE 
# LOW ENOUGH TO AVOID DYNAMIC EFFECTS. THE ACCURACY OF THE FORCE-DISPLACEMENT RESPONSE IS HIGHLY DEPENDENT
# ON THE CUTTER BLADE VELOCITY. HERE , THE CUTTING VELOCITY IS SET TO A HIGH VALUE FOR DEMONESTRATION PURPOSES ONLY. 
angularVelocity=0 # THE CUTTING SEQUENCE HAVE BEEN MODELED USING THE YADE CombinedKinematicEngin. SINCE THERE IS NO ROTATION
# FOR THE BLADE,  THE ANGULAR VELOCITY OF THE RotaionEngine OF THIS CombinedKinematicEngin IS SET TO ZERO.
PenetrationvelocityofBlade = 0.07 # PENETRATION VELOCITY. IN ORDER TO ENSURE THE NUMERICAL STABILITY OF THE MODEL, THE VELOCITY MUST BE 
# LOW ENOUGH TO AVOID DYNAMIC EFFECTS. THE ACCURACY OF THE FORCE-DISPLACEMENT RESPONSE IS HIGHLY DEPENDENT
# ON THE CUTTER BLADE VELOCITY. . HERE , THE PENETRATION VELOCITY IS SET TO A HIGH VALUE FOR DEMONESTRATION PURPOSES ONLY. 

############################################
###   DEFINING VARIABLES AND MATERIALS   ###
############################################
SoilId=FrictMat(young=youngSoil,poisson=poissonSoil,frictionAngle=math.radians(deposFricDegree),density=densSoil)
O.materials.append(SoilId)

ContainerId=FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=math.radians(0),density=densContainer)
O.materials.append(ContainerId)

###################################
#####   CREATING GEOMETRIES   #####
###################################
## DEM BOUNDARY
XcenterofContainer = 0.0025
YcenterofContainer = 0.0054383384
BottomOfDEMGeometry=0.00241406
ZcenterofContainer = 10.0e-3
Container=yade.geom.facetBox((XcenterofContainer,YcenterofContainer,ZcenterofContainer), (XcenterofContainer,(YcenterofContainer-BottomOfDEMGeometry),ZcenterofContainer),wallMask=51,material=ContainerId)
O.bodies.append(Container)
## IMPORTING THE ALREADY COMPACTED BIN (USE YOUR OWN LOCAL MACHINE ADDRESS)
## CLUMPS
## 
ymport.textClumps("../Cutter Blade/W_Transition Zone/Importable/DEMPClumps.txt",shift=Vector3(0,0,0),material=SoilId)
## TRANSITION
O.bodies.append(ymport.text("../Cutter Blade/W_Transition Zone/Importable/DEMSpheres.txt",shift=Vector3(0,0,0),material=SoilId,color=Vector3(0.6,0.6,0.6)))
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"maxX:",maxX,"minY:",minY,"maxY:",maxY,"minZ:",minZ,"maxZ:",maxZ)

XcenterofBlade = WidthBlade/2
YcenterofBlade = maxY+(HeightBlade/2)
ZcenterofBlade = minZ+IniDistanceBladefromBoundary
Blade=yade.geom.facetBox((XcenterofBlade,YcenterofBlade,ZcenterofBlade), (WidthBlade/2,HeightBlade/2,ThicknessBlade/2),material=ContainerId)
O.bodies.append(Blade)

idsr = [w.id for w in Blade]
facets = [b for b in O.bodies if isinstance(b.shape,Facet)] # list of facets in simulation

############################
###   DEFINING ENGINES   ###
############################
gravityAcc=-9.81
## SAVING VTK (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
O.engines=[VTKRecorder(fileName='../Cutter Blade/W_Transition Zone/VTK/Output_CB_SC_W_Transition_GRC3A_n057' ,recorders=['all'],iterPeriod=5000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Facet_Aabb()], label="collider"),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=normalDamp,betas=shearDamp,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
  ),
  CombinedKinematicEngine(ids=idsr,label='combEngine') + TranslationEngine(translationAxis=(0,-1,0),velocity=velocity) +\
  RotationEngine(rotationAxis=(1,0,0), angularVelocity=0, rotateAroundZero=True, zeroPoint=(XcenterofBlade,YcenterofBlade,ZcenterofBlade)),  
  NewtonIntegrator(damping=numDamp,gravity=(0,gravityAcc,0)),
  PyRunner(command='BladePenetration()',iterPeriod=1,label='checker'),
  PyRunner(iterPeriod=100,command='history()',label='recorder'),
  PyRunner(iterPeriod=1000,command='eraseOffs()'),
  ]
O.dt = 1e-7
O.step()
# Get TranslationEngine and RotationEngine from CombinedKinematicEngine
transEngine, rotEngine = combEngine.comb[0], combEngine.comb[1]

def BladePenetration():
  transEngine.velocity = 2*velocity
  rotEngine.zeroPoint += Vector3(0,-1,0)*transEngine.velocity*O.dt
  print ("rotEngine.zeroPoint[1]-HeightBlade/2):",(rotEngine.zeroPoint[1]-HeightBlade/2))
  print ("maxY-InitialPeneterationofBlade:",maxY-InitialPeneterationofBlade)
  # Stop the penetration process
  if rotEngine.zeroPoint[1]-HeightBlade/2<=maxY-InitialPeneterationofBlade:
    checker.command='BladeHorizentalMove()'
__builtin__.BladePenetration=BladePenetration

def BladeHorizentalMove():
  transEngine.translationAxis=(0,0,1)
  transEngine.velocity = velocity
  rotEngine.angularVelocity = angularVelocity
  rotEngine.zeroPoint += Vector3(0,0,1)*velocity*O.dt
  print ("rotEngine.zeroPoint[2]+ThicknessBlade/2):",(rotEngine.zeroPoint[2]+ThicknessBlade/2))
  # Stop the cutting process
  if rotEngine.zeroPoint[2]+ThicknessBlade/2>=5/6*maxZ:
    O.pause()
__builtin__.BladeHorizentalMove=BladeHorizentalMove

def eraseOffs():
  for b in O.bodies:
    if isinstance(b.shape,Sphere):
      if b.isClumpMember:
         if ((maxY+HeightBlade*1.5)<b.state.pos[1] or b.state.pos[1]<minY):
            O.bodies.erase(b.id)
__builtin__.eraseOffs=eraseOffs

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
  Dz=rotEngine.zeroPoint[2] +ThicknessBlade/2
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
## SAVING NUMERICAL FORCE-DISPLACEMENT (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
  plot.saveDataTxt('../Cutter Blade/W_Transition Zone/Numerical Output/Output_CB_SC_W_Transition_GRC3A_n057.txt')
__builtin__.history = history

def vtkExport(i):
  from yade import export
  ## PLEASE DO NOT CHANGE THIS
  name = '/tmp/vol1_yade'
  export.VTKExporter(name,i).exportSpheres()
  export.VTKExporter(name,i).exportFacets()




