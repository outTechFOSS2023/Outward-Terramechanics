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
## THIS SCRIPT COMPRESS THE ALREADY GENERATED CLOUDS USING OUTWARD'S GRAIN GENERATION SCHEME (LOCATED IN INITIAL CLOUD GENERATION FOLDER)
## TO THE DESIRED SAMPLE SIZE FOR THE TRIAXIAL TEST.
## NOTE THAT THE INITIAL CLOUD IS GENERATED USING WEIGHT-VOLUME RELATIONSHIP FOR THE DESIRED SAMPLE POROSITY WHERE THE SOLID VOLUME
## IS GIVEN TO THE InputData (IN INITIAL CLOUD GENERATION FOLDER)

############################################
###########   IMPORT LIBRARIES   ###########
############################################
from yade import pack
from yade import plot
from yade import ymport
from yade import export
# import os
# from yade import mpy as mp

############################################
##########   DEFINING VARIABLES    #########
############################################
nRead=readParamsFromTable(
  sigma=100000,
  damp=0.4,
  stabilityThreshold=0.01,
  targetPorosity=0.42954934,
  RequiredV_Solid=0.00000005112800344827,
  compFricDegree = 0,
  beta=0.27,
  crushing=2.1,
  young=0.80e8,
  poisson=0.42,
  density=1672,
  unknownOk=True,
)

from yade.params import table
sigma=table.sigma # TO BE CHANGED ON QUEST BETWEEN FOUR VALUSE 0.25e5, 0.5e5, 1.0e5
targetPorosity = table.targetPorosity# THE POROSITY WE NEED FOR THE PACKING
RequiredV_Solid=table.RequiredV_Solid # REQUIRED SOLID VOLUME FOR THE DESIRED POROSITY
compFricDegree = table.compFricDegree #  CONSTANT CONTACT FRICTION DURING THE SAMPLE PREPARATION 
beta = table.beta #  SHAPE PARAMETER DURING SAMPLE PREPARATION (WILL CHANGE FOR DEVIATORIC LOADING)
crushing= table.crushing #  CRUSHING PARAMETER DURING SAMPLE PREPARATION (WILL CHANGE FOR DEVIATORIC LOADING)
damp=table.damp # NUMERICAL DAMPING COEFFICIENT
stabilityThreshold=table.stabilityThreshold # WE TEST UNBALANCEDFORCE AGAINST THIS VALUE IN DIFFERENT LOOPS (SEE BELOW)
young=table.young # CONTACT STIFFNESS FOR SOIL-SOIL INTERACTION
poisson=table.poisson # POISSON'S RATIO
density=table.density # DENSITY OF REGOLITH

############################################
#################  MATERIALS  ##############
############################################
## A SIMPLE FrictMat MATERIAL IS USED FOR BOTH THE REGOLITH AND CONTAINER
O.materials.append(JiangMat(young=young,poisson=poisson,frictionAngle=radians(compFricDegree),density=density,xIc=crushing,beta=beta,label='regolith'))

############################################
##################  WALLS  #################
############################################
## NOTE: THESE TWO VECTORS ARE DEFINED FROM MAXS AND MINS (LINES 68 T0 73) AS THE SIZES OF THE LOOSE CLOUD
mn,mx=Vector3(0.000409207,0.00051028,0.000494709),Vector3(0.122897053,0.12297747799999999,0.12290813499999999) 
wallIds=aabbWalls([mn,mx],thickness=0.005)
O.bodies.append(wallIds)

############################################
#################  GRAINS  #################
############################################
ymport.textClumps("/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Typical Samples/Sample_GRC-3PlusLMA-1_Triaxial_OutwardPSD_CLumps_GCSC2p5_TestA_Initial.txt",color=(0.8,0.6,0.3),shift=Vector3(0,0,0),material='regolith')
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
sampleLength=(maxX-minX)
sampleheight=(maxY-minY)
sampleWidth=(maxZ-minZ)
V_total_MinMax=(sampleLength*sampleheight*sampleWidth)
NPor=1-(getSpheresVolume()/V_total_MinMax)
print ("minX:",minX,"minY:",minY,"minZ:",minZ,"maxX:",maxX,"maxY:",maxY,"maxZ:",maxZ)
print ("sampleLength:",sampleLength,"sampleheight:",sampleheight,"sampleWidth:",sampleWidth,"V_total_MinMax:",V_total_MinMax,"Provided_Vs:",getSpheresVolume(),"NPor:",NPor)


############################
###   DEFINING ENGINES   ###
############################
triax=TriaxialStressController(
    goal1=-sigma,
    goal2=-sigma,
    goal3=-sigma,
    max_vel=30,
    dead=True,
    stressMask = 7,
    internalCompaction=False,   
)


from yade import qt
newton=NewtonIntegrator(damping=damp,gravity=(0,0,0))
O.engines=[
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
    InteractionLoop(
      [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
      [Ip2_JiangMat_JiangMat_JiangPhys()],
      [Law2_ScGeom_JiangPhys_Jiang(includeRollResistMoment=False,includeTwistResistMoment=False,label='rollingResistance')]
    ),
    ## We will use the global stiffness of each body to determine an optimal timestep (see https://yade-dem.org/w/images/1/1b/Chareyre&Villard2005_licensed.pdf)
    #GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.35),
    triax,
    #TriaxialStateRecorder(iterPeriod=100,file='WallStresses'),
    newton,
    PyRunner(command='sampleAlignment()',iterPeriod=100,label='checker'),
]
O.dt=0.7e-7
Gl1_Sphere.stripes=0

#######################################
##########   SAMPLE ALIGNMENT   #######
#######################################
def sampleAlignment():
    newton.gravity=(0,0,0)
    triax.dead=False
    ## DA DENOTES FOR DURING ALIGNMENT COMPACTION  
    meanStress_DA=triax.meanStress
    s11_DA=-triax.stress(triax.wall_right_id)[0]
    s22_DA=-triax.stress(triax.wall_top_id)[1]
    s33_DA=-triax.stress(triax.wall_front_id)[2]
    unb_DA=unbalancedForce()
    minX_DA=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
    minY_DA=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
    minZ_DA=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
    maxX_DA=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
    maxY_DA=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
    maxZ_DA=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
    sampleLength_DA=(maxX_DA-minX_DA)
    sampleheight_DA=(maxY_DA-minY_DA)
    sampleWidth_DA=(maxZ_DA-minZ_DA)
    V_total_DA=(sampleLength_DA*sampleheight_DA*sampleWidth_DA)
    NPor_DA=utils.porosity()
    print ("sampleLength_DA:",sampleLength_DA,"sampleheight_DA:",sampleheight_DA,"sampleWidth_DA:",sampleWidth_DA,"NPor_DA:",NPor_DA)
    if sampleLength_DA<=0.003494255:
        triax.wall_left_activated=False
        triax.wall_right_activated=False
    if sampleWidth_DA<0.003494255:
        triax.wall_back_activated=False
        triax.wall_front_activated=False
    if sampleheight_DA<=0.0073406:
        triax.wall_top_activated=False
        triax.wall_bottom_activated=False
    if sampleWidth_DA<=0.003494255 and sampleLength_DA<=0.003494255 and sampleheight_DA<=0.0073406:
        triax.dead=True
        ## NOW A TRUE TO SIZE SAMPLE IS EXPORTEDs
        export.textClumps("/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Typical Samples/Sample_GRC-3PlusLMA-1_Triaxial_OutwardPSD_CLumps_GCSC2p5_TestA_Aligned.txt")
        O.pause()













