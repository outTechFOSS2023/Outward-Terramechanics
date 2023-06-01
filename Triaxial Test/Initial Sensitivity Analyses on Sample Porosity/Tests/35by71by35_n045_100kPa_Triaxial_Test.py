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
## THIS SCRIPT SIMULATES THE TRIAXIAL COMPRESSION TESTS OF A GRC-3 SAMPLE
## AS A PART OF A SENSITVITY ANALYSIS ON THE SAMPLE POROSITY ON THE FORCE-DISPLACEMENT
## HISTORY. 
## MODEL GEOMETRY: 3D
## PARTICLE GEOMETRY: CLUMPS
## SAMPLE LENGTH: 3.50mm
## SAMPLE HEIGHT: 7.10mm
## SAMPLE WIDTH:  3.50mm
## POROSITY: n=0.45
## CONFINING PRESSURE: 100KPA
## MATERIAL MODEL: HERTZ MINDLIN MODEL (NO ROLLING RESISTANCE)
#*************************************************************************/

############################################
###########   IMPORT LIBRARIES   ###########
############################################
from yade import pack
from yade import plot
from yade import ymport
from yade import export

############################################
############  DEFINING VARIABLES  ##########
############################################
nRead=readParamsFromTable(
	FricDegree = 30, # CONTACT FRICTION DURING COMPACTION
  FinalFricDegree=30, # CONTACT FRICTION DURING SHEAR
  key='_triax_base_', # PUT THE SIMULATION'S TITLE HERE
	unknownOk=True)
from yade.params import table
key=table.key
targetPorosity = 0.45 #THE TARGET POROSITY FOR THE MEDIUM DENSE PACKING
FricDegree = table.FricDegree # CONTACT FRICTION FOR SAMPLE PREPARATION
FinalFricDegree=table.FinalFricDegree
rate=-0.1 # STRAIN RATE DURING DEVIATORIC LOADING (CHANGE THIS DEPENDING ON COMPUTATIONAL NEEDS)
stabilityThreshold=0.01 # STABILITY THERESHOLD FOR APPLICATION OF CONFINEMENT
youngRegolith=0.7e8 # CONTACT STIFFNESS FOR REGOLITH
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonRegolith=0.3 # POISSION'S RATIO FOR REGOLITH
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densRegolith=2650 # DENSITY FOR REGOLITH
densContainer=7850 # DENSITY FOR CONTAINER
normalDamp=0.7 # HERTZ-MINDLIN NORMAL DAMPING
shearDamp=0.7 # HERTZ-MINDLIN SHEAR DAMPING
numDamp=0.4 # NUMERICAL DAMPING
sigma=100000 # CONFINING PRESSURE FOR SAMPLE PREPARATION, CONFINEMENT, AND DEVIATORIC LOADING
gravAcc=(0,0,0) # NO GRAVITY IS APPLIED TO PREVENT SEGRAGATION
mn,mx=Vector3( 0,0,0),Vector3(0.0035,0.0071,0.0035) # CORNERS OF THE TRIAXIAL CONTAINER

############################################
#################  MATERIALS  ##############
############################################

## A SIMPLE FrictMat MATERIAL IS USED FOR BOTH THE REGOLITH AND CONTAINER
O.materials.append(FrictMat(young=youngRegolith,poisson=poissonRegolith,frictionAngle=radians(FricDegree),density=densRegolith,label='regolith'))
O.materials.append(FrictMat(young=youngRegolith,poisson=poissionContainer,frictionAngle=radians(0),density=densContainer,label='container'))

############################################
############   MODEL GEOMETRY   ############
############################################
## CREATING WALLS AROUND THE PACKING
walls=aabbWalls([mn,mx],material='container')
wallIds=O.bodies.append(walls)

## IMPORTING ALIGNED LOOSE PACKING OF THE REGOLITH
ymport.textClumps("/home/ngoudarzi/Desktop/Triaxial Test/Initial Sensitivity Analyses on Sample Porosity/Samples/35by71by35_n045_Relaxed_Triaxial_Sample",shift=Vector3(0,0,0),material='regolith')
## BSC DENOTES FOR BEFORE SAMPLE COMPACTION 
minX_BSC=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY_BSC=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ_BSC=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX_BSC=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY_BSC=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ_BSC=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
sampleLength_BSC=(maxX_BSC-minX_BSC)
sampleheight_BSC=(maxY_BSC-minY_BSC)
sampleWidth_BSC=(maxZ_BSC-minZ_BSC)
V_Solid=utils.getSpheresVolume ()
NPor_BSC=utils.porosity ()
print ("minX_BSC:",minX_BSC,"minY_BSC:",minY_BSC,"minZ_BSC:",minZ_BSC,"maxX_BSC:",maxX_BSC,"maxY_BSC:",maxY_BSC,"maxZ_BSC:",maxZ_BSC)
print ("sampleLength_BSC:",sampleLength_BSC,"sampleheight_BSC:",sampleheight_BSC,"sampleWidth_BSC:",sampleWidth_BSC)
print ("V_Solid:",V_Solid,"NPor_BSC:",NPor_BSC)
#DISPLAY SPHERES WITH 2 COLORS FOR BETTER VISULAZIATION 
Gl1_Sphere.stripes=0
if nRead==0: yade.qt.Controller(), yade.qt.View()

############################
###   DEFINING ENGINES   ###
############################

triax=TriaxialStressController(
  ## TRIAXIALSTRESSCONTROLLER WILL BE USED TO CONTROL STRESS AND STRAIN. IT CONTROLS PARTICLES SIZE AND PLATES POSITIONS.
  maxMultiplier=0, # SPHERES GROWING FACTOR (FAST GROWTH). WE SET TO ZERO BECAUSE THERE IS NO INTERNAL COMPACTION HERE
  finalMaxMultiplier=0, # SPHERES GROWING FACTOR (SLOW GROWTH). WE SET TO ZERO BECAUSE THERE IS NO INTERNAL COMPACTION HERE
  thickness = 0,
  goal1=-sigma,
  goal2=-sigma,
  goal3=-sigma,
  ## SWITCH STRESS/STRAIN CONTROL USING A BITMASK. WHAT IS A BITMASK, HUH?!
  ## SAY X=1 IF STESS IS CONTROLLED ON X, ELSE X=0. SAME FOR FOR Y AND Z, WHICH ARE 1 OR 0.
  ## THEN AN INTEGER UNIQUELY DEFINING THE COMBINATION OF ALL THESE TESTS IS: MASK = X*1 + Y*2 + Z*4
  ## TO PUT IT DIFFERENTLY, THE MASK IS THE INTEGER WHOSE BINARY REPRESENTATION IS XYZ, I.E.
  ## "100" (1) MEANS "X", "110" (3) MEANS "X AND Y", "111" (7) MEANS "X AND Y AND Z", ETC.
  stressMask = 7,
  wall_top_activated=False,
  wall_back_activated=False,
  wall_front_activated=False,
  wall_left_activated=False,
  wall_right_activated=False,
  wall_bottom_activated=False, 
  internalCompaction=False, # IF TRUE THE CONFINING PRESSURE IS GENERATED BY GROWING PARTICLES (INTERNAL COMPACTION). WE AVOID THIS
  dead=True,
  )
newton=NewtonIntegrator(damping=numDamp,gravity=gravAcc)
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(en=normalDamp,es=shearDamp)],
		[Law2_ScGeom_MindlinPhys_Mindlin()]),
	GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.45),
	triax,
	TriaxialStateRecorder(iterPeriod=100,file='WallStresses'+table.key),
	newton,
	PyRunner(command='stressRelaxation()',iterPeriod=100,label='checker')]

#######################################
########   STRESS RELAXATION   ########
#######################################
def stressRelaxation():
  triax.dead=False
  sampleLength=abs(O.bodies[0].state.pos[0]-O.bodies[1].state.pos[0])
  sampleHeight=abs(O.bodies[2].state.pos[1]-O.bodies[3].state.pos[1])
  sampleWidth=abs(O.bodies[4].state.pos[2]-O.bodies[5].state.pos[2])
  totalVolume=sampleLength*sampleHeight*sampleWidth
  #solidVolume=9.730844402040118e-08
  numericalPorosity=triax.porosity
  meanStress=triax.meanStress
  unb=unbalancedForce()
  sxx=-triax.stress(triax.wall_right_id)[0],
  syy=-triax.stress(triax.wall_top_id)[1],
  szz=-triax.stress(triax.wall_front_id)[2],
  print ("Numerical Porosity:",numericalPorosity,"Target Porosity:",targetPorosity)
  print ("Unbalanced Force: ",unb,"Mean Stress:",meanStress,"SXX:",sxx,"SYY:",syy,"SZZ:",szz)
  print ("sampleLength:",sampleLength,"sampleWidth:",sampleWidth,"sampleHeight:",sampleHeight)
  if unb<=0.001:
    checker.command='confiningPressure()'
#######################################
###   APPLYING CONFINING PRESSURE   ###
#######################################
def confiningPressure():
  ## DC DENOTES FOR DURING CONFINEMENT
  setContactFriction(radians(FinalFricDegree))
  triax.wall_top_activated=True
  triax.wall_bottom_activated=True
  triax.wall_right_activated=True
  triax.wall_left_activated=True
  triax.wall_front_activated=True
  triax.wall_back_activated=True
  meanStress_DC=triax.meanStress
  unb_DC=unbalancedForce()
  NPor_DC=triax.porosity
  print ("unb_DC: ",unb_DC,"meanStress_DC:",meanStress_DC,"NPor_DC:",NPor_DC)
  if unb_DC<stabilityThreshold and abs(-sigma-meanStress_DC)/sigma<0.001:
    ## AC DENOTES FOR AFTER CONFINEMENT
    NPor_AC=triax.porosity
    global e11_AC,e22_AC,e33_AC  
    e11_AC=-triax.strain[0]
    e22_AC=-triax.strain[1]
    e33_AC=-triax.strain[2]
    print ("NPor_AC:",NPor_AC,"initPor:",NPor_BSC)
    print ("e11_AC:",e11_AC,"e22_AC:",e22_AC,"e33_AC:",e33_AC)
    print ("###      SAMPLE CONFINED      ###")
    O.engines=O.engines+[PyRunner(command='addPlotData()',iterPeriod=100)]
    checker.command='deviatoricLoading()'

############################################
###   DEVIATORIC LOADING-LOADING PHASE   ###
############################################
def deviatoricLoading():
  triax.stressMask = 5
  triax.goal2=rate
  triax.goal1=-sigma
  triax.goal3=-sigma
  ## DL DENOTES FOR DEVIATORIC LOADING
  NPor_DL=triax.porosity
  e11_DL=-triax.strain[0]-e11_AC
  e22_DL=-triax.strain[1]-e22_AC
  e33_DL=-triax.strain[2]-e22_AC
  s11_DL=-triax.stress(triax.wall_right_id)[0]
  s22_DL=-triax.stress(triax.wall_top_id)[1]
  s33_DL=-triax.stress(triax.wall_front_id)[2]
  DevStress_DL=(s22_DL-s33_DL)
  print ("e22_DL:",e22_DL,"NPor_DL:",NPor_DL,"DevStress_DL:",DevStress_DL)
  if e22_DL>=0.4:
    triax.wall_top_activated=False
    triax.wall_bottom_activated=False
    triax.dead=True
    O.pause()

############################################
##########   RECORDING VARIABLES  ##########
############################################
def addPlotData():
  ## DL DENOTES FOR DEVIATORIC LOADING
  global e11_DL,e22_DL,e33_DL,s11_DL,s22_DL,s33_DL,DevStress_DL,NPor_DL,NVoid_DL
  i=O.iter
  s11_DL=-triax.stress(triax.wall_right_id)[0]
  s22_DL=-triax.stress(triax.wall_top_id)[1]
  s33_DL=-triax.stress(triax.wall_front_id)[2]
  DevStress_DL=abs(s22_DL-s33_DL)
  e11_DL=-triax.strain[0]-e11_AC
  e22_DL=-triax.strain[1]-e22_AC
  e33_DL=-triax.strain[2]-e22_AC
  NPor_DL=triax.porosity  
  NVoid_DL=NPor_DL/(1-NPor_DL)
  yade.plot.addData({'e22':e22_DL,'DevStress':DevStress_DL,'NPor_DL':NPor_DL,})
  plot.saveDataTxt("/home/ngoudarzi/Desktop/Triaxial Test/Initial Sensitivity Analyses on Sample Porosity/Output/35by71by35_n045_100kPa_Triaxial_Output.txt") 
## DECLARE WHAT IS TO PLOT. "NONE" IS FOR SEPARATING Y AND Y2 AXIS
plot.plots={'e22':('DevStress',None,'NPor_DL')}
# DISPLAY ON THE SCREEN (DOESN'T WORK ON VMWARE IMAGE IT SEEMS)
plot.plot() 
#####  PLAY THE SIMULATION HERE WITH "PLAY" BUTTON OR WITH THE COMMAND O.run(N)  #####
#O.run(5000000000,True)        
    
