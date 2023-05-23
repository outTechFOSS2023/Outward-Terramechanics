# -*- coding: utf-8 -*-
#*********************************************************************************************************
#*********************************************************************************************************
# SBIR Rights Notice (DEC 2007)
# These SBIR data are furnished with SBIR rights under Contract No. 80NSSC20C0213. 
#For a period of 20 years, unless extended in accordance with FAR 27.409(h), after acceptance 
#of all items to be delivered under this contract, the Government will use these data for Government 
#purposes only, and they shall not be disclosed outside the Government (including disclosure for 
#procurement purposes) during such period without permission of the Contractor, except that, 
#subject to the foregoing use and disclosure prohibitions, these data may be disclosed for use by 
#support Contractors. After the protection period, the Government has a paid-up license to use, 
#and to authorize others to use on its behalf, these data for Government purposes, but is relieved 
#of all disclosure prohibitions and assumes no liability for unauthorized use of these data by third parties. 
#This notice shall be affixed to any reproductions of these data, in whole or in part.
# (End of Notice)
#*********************************************************************************************************
#*********************************************************************************************************
## This script models a typical free-surface cutter blade (CB) test on GRC-3. The bed size is considered 
## as a full width bed. The dimensions of the bed are l=5.0cm, w=5.0cm, and h=3.0cm.
## The simulation is conducted under non-reduced condition meaning the whole length of the bed
## is always present during the simulation and no adaptive reduced order modeling techniques are 
## applied.
## the bed has been prepared in its minimum dry density (maximum porosity) by uniaxial compression in 
## in the height direction and zero interparticle friction. The final porosity 
## after the bed stabilization is almost 0.428993526925169.
## The GRC-3+LMA-1 particles are clumps with no rolling/twisting resistance enabled (to include natural shape irregularities 
## by clumps
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
densSoil=1503.46004360603 # BULK DENSITY FOR SOIL
densContainer=7850 # DENSITY FOR CONTAINER
etaSoil=0.1 # HM ROLLING RESISTANCE PARAMETER
numDamp=0.4 # NUMERICAL DAMPING
IniDistanceBladefromBoundary = 0.01 # INITIAL DISTANCE OF THE RIGHT BOUNDARY OF THE CUTTER BLADE FROM THE LEFT BOUNDARY OF THE COMPACTED BED (-Z). 
HeightBlade=0.07 # CUTTER BLADE HEIGHT (Y DIRECTION)
WidthBlade=0.02 # CUTTER BLADE WIDTH (Z DIRECTION) 
ThicknessBlade=0.0015 # CUTTER BLADE THICKNESS (X DIRECTION)
InitialPeneterationofBlade = 0.010 # PENETRATION DEPTH OF THE BLADE.
InitialDistanceofBladefromTopSoil = 0 # INITIAL DISTANCE OF THE BOTTOM OF THE CUTTER BLADE FROM THE SURFACE OF THE COMPACTED BED. SINCE THE ASSEMBLY
# IS COMPLETELY SMOOTH, THIS PARAMETER CAN BE SET TO ZERO.
# NOTE: THE NUMERICAL STABILITY AND RESPNSE ACCURACY IS HIGHLY DEPENDENT UPON THE SIZE OF THE ACTIVE ZONE (THE LARGER, THE BETTER) AND THE BIN SIZE
# (THE SMALLER, THE BETTER).
HorizentalvelocityofBlade = 0.25e-2 #10 TIMES THE EXPERIMENTAL VELOCITY
PenetrationvelocityofBlade = 0.25e-2 #10 TIMES THE EXPERIMENTAL VELOCITY
WidthSample=0.050 # DEM SAMPLE WIDTH (Z)
HeightSample=0.030 # DEM SAMPLE HEIGHT (Y)
LengthSample=0.050 # DEM SAMPLE LENGTH (X)

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
ymport.textClumps("/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Simulations with Calibrated Microparameters/GRC-3/Loose/Samples/LSP-GRC3-B1-A-S2.txt",scale=10,shift=Vector3(0,0,0),material=SoilId,color=Vector3(0,0.5,0.7))
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
ymport.textClumps("/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Simulations with Calibrated Microparameters/GRC-3/Loose/Samples/LSP-GRC3-B1-A-S2.txt",scale=10,shift=Vector3(0,0,maxZ),material=SoilId,color=Vector3(0,0.5,0.7))
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"maxX:",maxX,"minY:",minY,"maxY:",maxY,"minZ:",minZ,"maxZ:",maxZ)

mn,mx=Vector3(minX,minY,minZ),Vector3(maxX,0.3,maxZ)

## CREATE WALLS AROUND THE PACKING
walls=aabbWalls([mn,mx],thickness=0,material=ContainerId,oversizeFactor=1.0)
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
O.engines=[VTKRecorder(fileName='/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Simulations with Calibrated Microparameters/GRC-3/Loose/VTK/LSP-GRC3-B1-A-S2/CSP-GRC3-B1-A-S2' ,recorders=['all'],iterPeriod=5000),
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
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Ordinary Large Scale Tests/Cutter Blade/Simulations with Calibrated Microparameters/GRC-3/Loose/Numerical Output/LSP-GRC3-B1-A-S2/LSP-GRC3-B1-A-S2.txt')
  
## DECLARE WHAT IS TO PLOT ON SCREEN. "None" IS FOR SEPARATING Y AND Y2 AXIS.
plot.plots={'i':('Fx','Fy','Fz','Dx','Dy','Dz')}
plot.plot()
# UNCOMMENT FOR NON-INTERACTIVE SIMULATION
#O.run(2000000000000,True)

