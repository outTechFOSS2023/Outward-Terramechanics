#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

######################################################################
# This Python script is the YADE (DEM) part of a FEMxDEM simulation of a soil wheel interaction (SWI) test.
# The FEM part of the simulation is defined in a diferent Python script entitled vol1_oofem.py.
# There are also some other dependencies which link the DEM and FEM for such adjustments of  as output request and number of iterations.
# NOTE: ALL EXTERNAL FILES (DEM GEOMETRY), AND OUTPUT FOLDERES ARE PROVIDED. THE CORRECT ADDRESSES
# MUST BE PROVIDED BY THE USER FROM THEIR OWN LOCAL MACHINES. 
# The DEM model has the following chracteristics:
# 1- No microparameter calibrations were conducted and the model has been setup to show the runability of  FEMxDEM.
# 2- The DEM sample is an assembly of  clumps used to respresent  a medium dense state of  GRC-3+LMA-1 (n=0.418).
# 3- DEM assembly size is 5cm(X-direction)*3.0cm(Y-direction)x10cm(Z-Direction) supported by a 3mm thick underlying FEM mesh. No 
# transition  zone between the FEM and clumps was used. NOTE: To accelerate the simulation, a scale factor of 10 has been applied during the 
# import of the assembly. 
# 4- Surface coupling (SC) FEM*DEM has been implemented. SC has been shown to run faster than volumetic coupling (VC). 
# 5- Wall boundaries and the smooth wheel has been generated using YADE built-in geometry generation functions for facets and therefore are not deformable. 
# However, it is possible to extend the FEM to all external objects.
# 6-  An advanced rolling /twisting resistance model (Hertz-Mindlin with rolling resistance--see YADE manual for the activation of the rolling resistance) constituitive law has been used for both regolith, boundaries, and the cutter blade where
# the rolling resistance parameters have been disabled for the boundaries and the wheel. 
# 7-  The whole process starts with the penetration of the wheel to 1mm into the regolith and then the
# rolling over the regolith in the +Z direction.
# 8- The penetration and rolling processes are simulated using YADE CombinedKinematicEngine. 
######################################################################
import math
from math import sqrt
from math import pi
from libyade import yade
from yade import *
from yade import pack, plot 
import numpy
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
## ANGULAR VELOCITY: rad/s
## TIME: s
## PRESSURE: pa
## DENSITY: kg/m3
## Wheel PROPERTIES
widthWheel=0.0125 #WIDTH OF THE SMOOTH WHEEL
Wheel_R=0.015 # RADIUS OF THE SMOOTH WHEEL
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
deposFricDegree = 20.6  # CONTACT FRICTION
normalDamp=0.7 # NORMAL VISCOUS DAMPING
shearDamp=0.7 # SHEAR VISCOUS DAMPING
youngSoil=0.489e8 # CONTACT STIFFNESS FOR SOIL
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.393 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densSoil=1801.11858943976 # DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.4 # NUMERICAL DAMPING
initialPeneterationofWheel=0.01 # PENETRATION DEPTH OF THE WHEEL BEFORE THE START OF THE ROLLING PROCESS
InitialDistanceofWheelfromTopSoil=0.0001 # INITIAL OFFSET OF THE BOTTOM OF THE WHEEL FROM THE TOP OF THE SAMPLE
iniDistanceWheelfromBoundary=0.015 #INITIAL DISTANCE OF THE WHEEL CENTER FROM THE LEFT BOUNDARY DURING PENETRATION
WidthSample=0.05 # X DIRECTION
HeightSample=0.03 # Y DIRECTION
LengthSample=0.1 # Z Direction
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
XcenterofContainer = WidthSample/2
YcenterofContainer = HeightSample
#BottomOfDEMGeometry=0.00241406
ZcenterofContainer = LengthSample/2
Container=yade.geom.facetBox((XcenterofContainer,YcenterofContainer,ZcenterofContainer), (XcenterofContainer,YcenterofContainer,ZcenterofContainer),wallMask=51,material=ContainerId)
O.bodies.append(Container)
## IMPORTING THE ALREADY COMPACTED BIN TWICE
## CLUMPS
## IMPORTING THE ALREADY COMPACTED BIN (USE YOUR OWN LOCAL MACHINE ADDRESS)
ymport.textClumps("../SWI/WO_Transition Zone/Importable/CSP-GRC3A-B1-A.txt",scale=10,shift=Vector3(0,0,0),material=SoilId)
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
ymport.textClumps("../SWI/WO_Transition Zone/Importable/CSP-GRC3A-B1-A.txt",scale=10,shift=Vector3(0,0,maxZ),material=SoilId)
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"maxX:",maxX,"minY:",minY,"maxY:",maxY,"minZ:",minZ,"maxZ:",maxZ)
## WHEEL
XcenterofWheel = WidthSample/2
YcenterofWheel = Wheel_R+maxY+InitialDistanceofWheelfromTopSoil
ZcenterofWheel = iniDistanceWheelfromBoundary+minZ
Wheel1 = yade.geom.facetCylinder((XcenterofWheel,YcenterofWheel,ZcenterofWheel),radius=Wheel_R,height=widthWheel,wallMask=7,segmentsNumber=200,material=ContainerId,orientation=yade.Quaternion((0,1,0),-pi/2),wire=True,color=Vector3(255,255,255))
O.bodies.append(Wheel1)
idsr = [w.id for w in Wheel1]
facets = [b for b in O.bodies if isinstance(b.shape,Facet)] # list of facets in simulation
for b in facets:
  O.forces.setPermF(b.id,(0,Wheel_AdditionalNormalForce/800,0))
############################
###   DEFINING ENGINES   ###
############################
O.engines=[VTKRecorder(fileName='../SWI/WO_Transition Zone/VTK/Output_SWI_SC_W_Transition_GRC3A_n0418' ,recorders=['all'],iterPeriod=5000),
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
  PyRunner(iterPeriod=1, command="WheelPenetration()" ,label='checker'),
  PyRunner(iterPeriod=100,command='history()',label='recorder'),
  ]
O.dt = 1e-7
O.step()

# Get TranslationEngine and RotationEngine from CombinedKinematicEngine
transEngine, rotEngine = combEngine.comb[0], combEngine.comb[1]

def WheelPenetration():
  transEngine.velocity = 2*velocity
  rotEngine.zeroPoint += Vector3(0,-1,0)*transEngine.velocity*O.dt
  print ("rotEngine.zeroPoint[1]-Wheel_R):",(rotEngine.zeroPoint[1]-Wheel_R))
  print ("maxY-initialPeneterationofWheel):",(maxY-initialPeneterationofWheel))
  if (rotEngine.zeroPoint[1]-Wheel_R)<=(maxY-initialPeneterationofWheel):
    checker.command='WheelRolling()'
__builtin__.WheelPenetration = WheelPenetration
def WheelRolling():
  transEngine.translationAxis=(0,0,1)  
  transEngine.velocity = velocity
  rotEngine.angularVelocity = angularVelocity
  rotEngine.zeroPoint += Vector3(0,0,1)*velocity*O.dt
  print ("rotEngine.zeroPoint[2]+Wheel_R):",(rotEngine.zeroPoint[2]+Wheel_R),"Iteration:",O.iter)
  if rotEngine.zeroPoint[2]>=maxZ-1.05*Wheel_R:
    O.pause()
__builtin__.WheelRolling = WheelRolling


def history():
  global Fx,Fy,Fz,Dx,Dy,Dz
  Fx=0
  Fy=0
  Fz=0
  for b in facets:
    Fx+=abs(O.forces.f(b.id,sync=True)[0])
    Fy+=abs(O.forces.f(b.id,sync=True)[1])
    Fz+=abs(O.forces.f(b.id,sync=True)[2])
  Dx=rotEngine.zeroPoint[0]
  Dy=rotEngine.zeroPoint[1]
  Dz=rotEngine.zeroPoint[2]+ Wheel_R
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
  ## In that case we can still save the data to a text file at the the end of the simulation, with:
  plot.saveDataTxt('../SWI/WO_Transition Zone/Numerical Output/Output_SWI_SC_W_Transition_GRC3A_n0418.txt')
__builtin__.history = history

def vtkExport(i):
  from yade import export
  name = '/tmp/vol1_yade'
  export.VTKExporter(name,i).exportSpheres()
  export.VTKExporter(name,i).exportFacets()




