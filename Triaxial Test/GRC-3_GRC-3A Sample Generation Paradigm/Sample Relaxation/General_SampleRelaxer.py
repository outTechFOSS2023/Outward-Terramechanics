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
## THIS SCRIPT RELAXES THE EXCESS STRESSES GENERATED WITHIN THE TRIAXIAL SAMPLE DUE TO MULTIAXIAL COMPACTION SO THAT
## THE FINAL PRODUCT WILL BE AT LOWER STRESS STATES MAKING IT SUITABLE FOR LOWER CONFINING PRESSURES. 
#*************************************************************************/

from yade import pack,export
from yade import pack, plot 
from yade import ymport

############################################
##########   DEFINING VARIABLES    #########
############################################
nRead=readParamsFromTable(
  sigma=0,
  NSTEPS=100000,
  damp=0.4,
  stabilityThreshold=0.01,
  targetPorosity=0.368686869,
  RequiredV_Solid=0.00000006437685095369,
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
NSTEPS=table.NSTEPS # NUMBER OF STEPS TO STABILIZE
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
mn,mx=Vector3(0,0,0),Vector3(0.003604511,0.0078486,0.003604511) # CORNERS OF THE INITIAL PACKING.
############################################
#################  MATERIALS  ##############
############################################
## A SIMPLE FrictMat MATERIAL IS USED FOR BOTH THE REGOLITH AND CONTAINER
O.materials.append(JiangMat(young=young,poisson=poisson,frictionAngle=radians(compFricDegree),density=density,xIc=crushing,beta=beta,label='regolith'))
#O.materials.append(FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=radians(0),density=densContainer,label='walls'))

############################################
##################  WALLS  #################
############################################
wallIds=aabbWalls([mn,mx],thickness=0.0)
O.bodies.append(wallIds)
ymport.textClumps("/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Typical Samples/Sample_GRC-3PlusLMA-1_Triaxial_OutwardPSD_CLumps_GCSC2p5_TestA_Aligned.txt",color=(0.3,0.8,0.1),shift=Vector3(0,0,0),material='regolith')
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
V_T=(mx[0]*mx[1]*mx[2])
V_S=getSpheresVolume()
Porosity=1-(V_S/V_T)
# sampleLength=(maxX-minX)
# sampleheight=(maxY-minY)
# sampleWidth=(maxZ-minZ)
# V_total_MinMax=(sampleLength*sampleheight*sampleWidth)
print ("minX:",minX,"minY:",minY,"minZ:",minZ,"maxX:",maxX,"maxY:",maxY,"maxZ:",maxZ)
print ("V_T:",V_T,"V_S:",V_S,"Porosity:",Porosity)
xPos=((maxX+minX)/2)
yPos=((maxY+minY)/2)
zPos=((maxZ+minZ)/2)
for b in O.bodies:
  if isinstance(b.shape,Sphere):
    O.bodies.erase(b.id)
ymport.textClumps("/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Typical Samples/Sample_GRC-3PlusLMA-1_Triaxial_OutwardPSD_CLumps_GCSC2p5_TestA_Aligned.txt",color=(0.3,0.8,0.1),shift=Vector3((mx[0]/2-xPos),(mx[1]/2-yPos),(mx[2]/2-zPos)),material='regolith')

############################
###   DEFINING ENGINES   ###
############################
Triax=TriaxialStressController(
  ## TRIAXIALSTRESSCONTROLLER WILL BE USED TO CONTROL STRESS AND STRAIN. IT CONTROLS PARTICLES SIZE AND PLATES POSITIONS.
  thickness = 0.00,
  ## SWITCH STRESS/STRAIN CONTROL USING A BITMASK. WHAT IS A BITMASK, HUH?!
  ## SAY X=1 IF STESS IS CONTROLLED ON X, ELSE X=0. SAME FOR FOR Y AND Z, WHICH ARE 1 OR 0.
  ## THEN AN INTEGER UNIQUELY DEFINING THE COMBINATION OF ALL THESE TESTS IS: MASK = X*1 + Y*2 + Z*4
  ## TO PUT IT DIFFERENTLY, THE MASK IS THE INTEGER WHOSE BINARY REPRESENTATION IS XYZ, I.E.
  ## "100" (1) MEANS "X", "110" (3) MEANS "X AND Y", "111" (7) MEANS "X AND Y AND Z", ETC.
  stressMask = 7,
  internalCompaction=False, # IF TRUE THE CONFINING PRESSURE IS GENERATED BY GROWING PARTICLES (INTERNAL COMPACTION). WE AVOID THIS
  wall_left_activated=False,
  wall_right_activated=False,
  wall_back_activated=False,
  wall_front_activated=False,
  wall_top_activated=False,
  wall_bottom_activated=False,
  )

Newton=NewtonIntegrator(damping=damp)
O.engines=[
  ForceResetter(),
  InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
  InteractionLoop(
    [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
    [Ip2_JiangMat_JiangMat_JiangPhys()],
    [Law2_ScGeom_JiangPhys_Jiang(includeRollResistMoment=False,includeTwistResistMoment=False,label='rollingResistance')]
  ),
  ## WE WILL USE THE GLOBAL STIFFNESS OF EACH BODY TO DETERMINE AN OPTIMAL TIMESTEP (SEE https://yade-dem.org/w/images/1/1b/Chareyre&Villard2005_licensed.pdf)
  GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.5),
  Triax,
  #TriaxialStateRecorder(iterPeriod=100,file='WallStresses'+table.key),
  Newton,
  PyRunner(command='confiningPressure()',iterPeriod=100,label='checker'),
]
Triax.goal1=Triax.goal2=Triax.goal3=-sigma
#######################################
###   APPLYING CONFINING PRESSURE   ###
#######################################
def confiningPressure():
  global V_T
  global V_S
  Porosity=1-(V_S/V_T)
  unb=unbalancedForce()
  print ('unbalanced force:',unb,' mean stress: ',Triax.meanStress)
  print ("ProvidedV_Solid:",utils.getSpheresVolume(),"Porosity:",Porosity)
  if O.iter>NSTEPS and unb<stabilityThreshold:
    #O.save('/home/ngoudarzi/Desktop/GRC-3 Calibration/Outward PSD/Triaxial/MTL with_Internal Confining Pressure/4p5x9p0x4p5_IncludingSG/confinedState_'+key+'.yade.gz')
    print ("###      ISOTROPIC STATE SAVED       ###")
    export.textClumps("/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Typical Samples/Sample_GRC-3PlusLMA-1_Triaxial_OutwardPSD_CLumps_GCSC2p5_TestA_Relaxed.txt.txt")
    O.pause()


 

