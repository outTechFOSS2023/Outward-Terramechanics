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
## THIS SCRIPT SIMULATES THE TRIAXIAL COMPRESSION TESTS OF A FILLITE SAMPLE
## AS A PART OF A SENSITVITY ANALYSIS ON THE SAMPLE SIZE ON THE FORCE-DISPLACEMENT
## HISTORY. 
## MODEL GEOMETRY: 3D
## PARTICLE GEOMETRY: SPHERICAL
## SAMPLE LENGTH: 3.200mm
## SAMPLE HEIGHT: 6.400mm
## SAMPLE WIDTH:  3.200mm
## POROSITY: n=0.45
## MATERIAL MODEL: JIANG ROLLING RESISTANCE MODEL
#*************************************************************************/
from yade import pack
from yade import ymport
import numpy
import math
from yade import plot

############################################
###   DEFINING VARIABLES AND MATERIALS   ###
############################################
deposFricDegree = 30 # CONTACT FRICTION DEGREE
stabilityThreshold=0.01 # STABILITY THERESHOLD FOR APPLICATION OF CONFINEMENT
youngRegolith=0.7e8# CONTACT STIFFNESS FOR REGOLITH
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonRegolith=0.3 # POISSION'S RATIO FOR REGOLITH
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densRegolith=2650 # DENSITY FOR REGOLITH
densContainer=7850 # DENSITY FOR CONTAINER
gravAcc=(0,-9.81,0) # NO GRAVITY IS APPLIED TO PREVENT SEGRAGATION
numDamp=0.4 # NUMERICAL DAMPING
targetPorosity=0.45 # TARGET POROSITY FOR THE SAMPLE PREPARATION
sigma=100000 # CONFINING PRESSURE
rateDeviatoric=-0.05 # STRAIN RATE FOR THE DEVIATORIC LOADING (SUBJECT TO SENSITIVITY ANALYSIS TO ENSURE QUASI-STATIC SIMULATION)
damp=0.4 # NUMERICAL DAMPING
stabilityThreshold=0.01 # STABILITY THERESHOLD FOR APPLICATION OF CONFINING PRESSURE
mn,mx=Vector3(0.003394630176313126,0.001779094823135192,0.0033714370025745254),Vector3(0.006594025517912851, 0.008174811752523948,0.0065711956304041535) # CORNERS OF THE INITIAL PACKING

## MATERIALS FOR GRAINS AND BOUNDARIES
O.materials.append(JiangMat(young=youngRegolith,poisson=poissonRegolith,frictionAngle=radians(deposFricDegree),density=densRegolith,xIc=2.1,beta=0.0,label='spheres'))

## BOUNDARIES
walls=aabbWalls([mn,mx],thickness=0)
wallIds=O.bodies.append(walls)

O.bodies.append(ymport.text("/home/ngoudarzi/Desktop/Triaxial Test/Size Effect for Fillite Sample/Samples/3p2mmx6p4mmx3p2mm/3p2mmx6p4mmx3p2mm_n045_Fillite.txt",shift=Vector3(0,0,0),material='spheres'))
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"minY:",minY,"minZ:",minZ,"maxX:",maxX,"maxY:",maxY,"maxZ:",maxZ)
n_YADE=utils.porosity ()
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("\r n_YADE: ",n_YADE,"maxY: ",maxY)
############################
###   DEFINING ENGINES   ###
############################

triax=TriaxialStressController(
  ## TriaxialStressController will be used to control stress and strain. It controls particles size and plates positions.
  ## this control of boundary conditions was used for instance in http://dx.doi.org/10.1016/j.ijengsci.2008.07.002
  thickness = 0,
  ## switch stress/strain control using a bitmask. What is a bitmask, huh?!
  ## Say x=1 if stess is controlled on x, else x=0. Same for for y and z, which are 1 or 0.
  ## Then an integer uniquely defining the combination of all these tests is: mask = x*1 + y*2 + z*4
  ## to put it differently, the mask is the integer whose binary representation is xyz, i.e.
  ## "100" (1) means "x", "110" (3) means "x and y", "111" (7) means "x and y and z", etc.
  stressMask = 7,
  internalCompaction=False, # If true the confining pressure is generated by growing particles
)

newton=NewtonIntegrator(damping=numDamp)

O.engines=[
  ForceResetter(),
  InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
  InteractionLoop(
    [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
    [Ip2_JiangMat_JiangMat_JiangPhys()],
    [Law2_ScGeom_JiangPhys_Jiang(includeRollResistMoment=False,includeTwistResistMoment=False)]
  ),
  ## We will use the global stiffness of each body to determine an optimal timestep (see https://yade-dem.org/w/images/1/1b/Chareyre&Villard2005_licensed.pdf)
  GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8),
  triax,
  TriaxialStateRecorder(iterPeriod=100,file='WallStresses'),
  newton,
  PyRunner(command='confiningPressure()',iterPeriod=100,label='checker'),
]
#Display spheres with 2 colors for seeing rotations better
#Gl1_Sphere.stripes=0
#if nRead==0: yade.qt.Controller(), yade.qt.View()
#######################################
###   APPLYING CONFINING PRESSURE   ###
#######################################
def confiningPressure():
  global e11_AC, e22_AC, e33_AC
  meanStress=triax.meanStress
  triax.goal1=triax.goal2=triax.goal3=-sigma
  unb=unbalancedForce()
  print ("Unbalanced Force: ",unb,"Mean Stress:",meanStress)
  if unb<stabilityThreshold and abs(-sigma-triax.meanStress)/sigma<0.001:
    NPor_AC=triax.porosity
    e11_AC=-triax.strain[0]
    e22_AC=-triax.strain[1]
    e33_AC=-triax.strain[2]
    print ("Numerical porosity after confinement is:",NPor_AC,"Target porosity after sample preparation is:",targetPorosity)
    print ("Sample axial strain after confinement:",e22_AC)
    print ("###      SAMPLE CONFINED      ###")
    O.engines=O.engines+[PyRunner(command='addPlotData()',iterPeriod=100)]
    #O.engines=O.engines+[PyRunner(command='addContourData()',iterPeriod=100)]
    checker.command='deviatoricLoading()'
############################################
###   DEVIATORIC LOADING-LOADING PHASE   ###
############################################
def deviatoricLoading():
  triax.stressMask = 5
  triax.goal2=rateDeviatoric
  triax.goal1=-sigma
  triax.goal3=-sigma
  NPor_L=triax.porosity
  e11_L=-triax.strain[0]-e11_AC
  e22_L=-triax.strain[1]-e22_AC
  e33_L=-triax.strain[2]-e33_AC
  s11_L=-triax.stress(triax.wall_right_id)[0]
  s22_L=-triax.stress(triax.wall_top_id)[1]
  s33_L=-triax.stress(triax.wall_front_id)[2]
  DevStress_L=(s22_L-s33_L)
  print ('e22_L:',e22_L)
  print ("Numerical porosity during deviatoric loading is:",NPor_L,"Deviator stress during loading is:",DevStress_L)
  if e22_L>0.35:
    print ("###      SAMPLE SHEARED     ###")
    wall_top_activated=False
    wall_bottom_activated=False
    triax.dead=True
    O.pause()
############################################
##########   RECORDING VARIABLES  ##########
############################################
def addPlotData():
  global e11,e22,e33,s11,s22,s33,DevStress,MobFrict,NPor,Nvoid
  i=O.iter
  s11=-triax.stress(triax.wall_right_id)[0]
  s22=-triax.stress(triax.wall_top_id)[1]
  s33=-triax.stress(triax.wall_front_id)[2]
  DevStress=(s22-s33)
  MobFrict=(s22-s33)/(s22+s33)
  e11=-triax.strain[0]
  e22=-triax.strain[1]
  e33=-triax.strain[2]
  NPor=triax.porosity
  Nvoid=NPor/(1-NPor)
  CN=avgNumInteractions()
  yade.plot.addData({'e11':e11,'e22_dev':e22,'e22_MobFrict':e22,'e22_NPor':e22,'e22_CN':e22,'e33':e33,'s11':s11,'s22':s22,'s33':s33,'DevStress':DevStress,'MobFrict':MobFrict,'NPor':NPor,'CN':CN,'i':i,})
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Triaxial Test/Size Effect for Fillite Sample/Numerical Output/3p2mmx6p4mmx3p2mm/3p2mmx6p4mmx3p2mm_n045_Fillite_DigitalCurves_Triaxial.txt')
  
## DECLARE WHAT IS TO PLOT. "NONE" IS FOR SEPARATING Y AND Y2 AXIS
plot.plots={'i':('e11','e22','e33',None,'s11','s22','s33')}
## THE TRADITIONAL TRIAXIAL CURVES WOULD BE MORE LIKE THIS:
plot.plots={'e22_dev':('DevStress',),'e22_MobFrict':('MobFrict',),'e22_NPor':('NPor',),'e22_CN':('CN',),}
# DISPLAY ON THE SCREEN (DOESN'T WORK ON VMWARE IMAGE IT SEEMS)
plot.plot() 
#####  PLAY THE SIMULATION HERE WITH "PLAY" BUTTON OR WITH THE COMMAND O.run(N)  #####
#O.run(5000000000,True)





