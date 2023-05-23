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
# 2- The DEM sample is an assembly of  clumps used to respresent  a medium dense state of  GRC-3+LMA-1 (n=0.418).
# 3- DEM assembly size is 7.68cm(X-direction)*4.0cm(Y-direction)x7.68cm(Z-Direction) supported by a 3mm thick underlying FEM mesh. No 
# transition  zone between the FEM and clumps was used. NOTE: To accelerate the simulation, a scale factor of 10 has been applied during the 
# import of the assembly. 
# 4- Surface coupling (SC) FEM*DEM has been implemented. SC has been shown to run faster than volumetic coupling (VC). 
# 5- Wall boundaries and the cutter blade has been generated using YADE built-in geometry generation functions for facets and therefore are not deformable. 
# However, it is possible to extend the FEM to all external objects.
# 6-  An advanced rolling /twisting resistance model (Hertz-Mindlin with rolling resistance--see YADE manual for the activation of the rolling resistance) constituitive law has been used for both regolith, boundaries, and the cutter blade where
# the rolling resistance parameters have been disabled for the boundaries and the cutter blade. 
# 7-  The curring process starts immidiately in the +Z direction meaning that the penetration step has been skipped to accelerate the simulation. 
# 8- The cutting process is simulated using YADE CombinedKinematicEngine where only the TranslationEngine component in +Z direction
# is active.  
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
deposFricDegree = 21.38 # CONTACT FRICTION DURING THE CUTTING PROCESS
normalDamp=0.7  # NORMAL VISCOUS DAMPING
shearDamp=0.7 # SHEAR VISCOUS DAMPING
youngSoil=0.489e8# CONTACT STIFFNESS FOR SOIL
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonSoil=0.393 # POISSION'S RATIO FOR SOIL
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
crushingSoil=2.1 # CRUSHING PARAMETER FOR SOIL
betaSoil=0.0 # SHAPE PARAMETER FOR SOIL
densSoil=1801.11858943976 # DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.4  # NUMERICAL DAMPING
WidthSample=0.0768 # X DIRECTION
HeightSample=0.0400 # Y DIRECTION
LengthSample=0.0768 # Z Direction
HeihtFEM=0.03 # CONSIDERED THICKNESS FOR THE UNDERLYING FEM
HeightBlade=0.08 # HEIGHT OF THE CUTTER BLADE (m) (Y-DIRECTION)
WidthBlade=0.02  # WIDTH OF THE CUTTER BLADE (m) (X-DIRECTION)
ThicknessBlade=0.000020 # THICKNESS OF THE CUTTER BLADE (m) (Z-DIRECTION)
InitialPeneterationofBlade = 0.01 # PENETRATION DEPTH OF THE BLADE BEFORE THE START OF THE CUTTING PROCESS
InitialDistanceofBladefromTopSoil = 0 # INITIAL OFFSET OF THE BOTTOM OF THE BLADE FROM THE TOP OF THE SAMPLE
velocity = 0.025 # CUTTING VELOCITY. IN ORDER TO ENSURE THE NUMERICAL STABILITY OF THE MODEL, THE VELOCITY MUST BE 
# LOW ENOUGH TO AVOID DYNAMIC EFFECTS. THE ACCURACY OF THE FORCE-DISPLACEMENT RESPONSE IS HIGHLY DEPENDENT
# ON THE CUTTER BLADE VELOCITY. HERE , THE CUTTING VELOCITY IS SET TO A HIGH VALUE FOR DEMONESTRATION PURPOSES ONLY. 
angularVelocity=0  # THE CUTTING SEQUENCE HAVE BEEN MODELED USING THE YADE CombinedKinematicEngin. SINCE THERE IS NO ROTATION
# FOR THE BLADE,  THE ANGULAR VELOCITY OF THE RotaionEngine OF THIS CombinedKinematicEngin IS SET TO ZERO.

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
ZcenterofContainer = LengthSample/2
## BOUNDARIES
Container=yade.geom.facetBox((XcenterofContainer,YcenterofContainer,ZcenterofContainer), (XcenterofContainer,YcenterofContainer,ZcenterofContainer),wallMask=51,material=ContainerId)
O.bodies.append(Container)
## IMPORTING THE ALREADY COMPACTED BIN
## CLUMPS
## IMPORTING THE ALREADY COMPACTED BIN (USE YOUR OWN LOCAL MACHINE ADDRESS)
ymport.textClumps("../Cutter Blade/WO_Transition Zone/7p68cmx4p0cmx7p68cm/Importable/CSP-GRC3A-B1-A.txt",scale=10,shift=Vector3(0,0,0),material=SoilId)
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
## BLADE
XcenterofBlade = WidthSample/2
YcenterofBlade = maxY+(HeightBlade/2)-InitialPeneterationofBlade
ZcenterofBlade = -ThicknessBlade/2
Blade=yade.geom.facetBox((XcenterofBlade,YcenterofBlade,ZcenterofBlade), (WidthBlade/2,HeightBlade/2,ThicknessBlade/2),material=ContainerId)
O.bodies.append(Blade)
idsr = [w.id for w in Blade]
facets = [b for b in O.bodies if isinstance(b.shape,Facet)] # list of facets in simulation

############################
###   DEFINING ENGINES   ###
############################
gravityAcc=-9.81
## SAVING VTK (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
O.engines=[VTKRecorder(fileName='../Cutter Blade/WO_Transition Zone/7p68cmx4p0cmx7p68cm/VTK/Output_CB_SC_WO_Transition_GRC3A_n0418' ,recorders=['all'],iterPeriod=5000),
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
  PyRunner(command='BladeHorizentalMove()',iterPeriod=100,label='checker'),
  PyRunner(iterPeriod=1000,command='history()',label='recorder'),
  ]
O.dt = 1e-7
O.step()
# Get TranslationEngine and RotationEngine from CombinedKinematicEngine
transEngine, rotEngine = combEngine.comb[0], combEngine.comb[1]

def BladeHorizentalMove():
  transEngine.translationAxis=(0,0,1)
  transEngine.velocity = velocity
  rotEngine.angularVelocity = angularVelocity
  rotEngine.zeroPoint += Vector3(0,0,1)*velocity*O.dt
  print ("rotEngine.zeroPoint[2]+ThicknessBlade/2):",(rotEngine.zeroPoint[2]+ThicknessBlade/2))
  if rotEngine.zeroPoint[2]+ThicknessBlade/2>=5/6*maxZ:
    O.pause()
__builtin__.BladeHorizentalMove=BladeHorizentalMove

def history():
  global Fx,Fy,Fz,Dx,Dy,Dz
  Fx=0
  Fy=0
  Fz=0
  for b in facets:
    Fx+=O.forces.f(b.id,sync=True)[0]
    Fy+=O.forces.f(b.id,sync=True)[1]
    Fz+=abs(O.forces.f(b.id,sync=True)[2])
  Dx=rotEngine.zeroPoint[0]
  Dy=rotEngine.zeroPoint[1]
  Dz=rotEngine.zeroPoint[2]+ThicknessBlade/2 
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
## SAVING NUMERICAL FORCE-DISPLACEMENT (USE YOUR OWN LOCAL MACHINE ADDRESS)--NAMING IS ARBITRARY
  plot.saveDataTxt('../Cutter Blade/W_Transition Zone/Numerical Output/7p68cmx4p0cmx7p68cm/Output_CB_SC_WO_Transition_GRC3A_n0418.txt')
__builtin__.history = history

def vtkExport(i):
  from yade import export
  name = '/tmp/vol1_yade'
  export.VTKExporter(name,i).exportSpheres()
  export.VTKExporter(name,i).exportFacets()




