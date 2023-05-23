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
#  THIS PROGRAM IS A FREE SOFTWARE; IT IS LICENCED UNDER THE TERMS OF THE 
#  GNU GENERAL PUBLIC LICENCE V2 OR LATER
#*************************************************************************/
#*****************************************************************************************************************
## THIS SCRIPT IS An INPUT FILE DETAILING THE SIMULATION OF NON-PERIODIC  PENETRATION AND DRAG IN GRC-3 USING A BUCKET 
## EXCAVATOR IN MD STATE (n~0.44) USING YADE. THE ROTATIONAL EXCAVATION STARTS WITHOUT PEE-PENETRATION OF THE 
## BUCKET INTO THE BED. THE DESIRED PENETRATION AND DRAG VELOCITY CAN BE ADJUSTED IN LINE 28 BUT FOR STABILITY 
## PURPOSES, LOWR VALUES<1 M/S ARE RECOMMENDED. THE USED CONSTITUIIVE LAW IS HERTZ-MINDLIN WITH SOME ROLLING 
## RESISTANCE IMPLEMENTED FOR BETTER STABILITY (LINE 68), SINCE SOME SPHERICAL PARTICLES ARE PRESENT IN THE MODEL. 
## FORCE HISTORY IN ALL THREE X, Y, AND Z ARE RECORDED AND VTK OF THE DEMOFRMATION IS STORED AS WELL. A COMBINED 
## KINEMATIC ENGINE CONTROLS THE INITIAL VERTICAL AND HORIZONTAL MOVEMENT 
#*****************************************************************************************************************

from yade import pack
from yade import plot
from yade import ymport
from yade import export
############################################
###           INPUT PARAMETERS           ###
############################################
velocity=1 #m/s THAT SHOULD BE 0.01 m/s
gravityAcc=-9.81
atRestFricDegree = 40 # INITIAL CONTACT FRICTION FOR SAMPLE PREPARATION
normalDamp=0.7 # NUMERICAL DAMPING
shearDamp=0.7
youngRegolith=0.7e8# CONTACT STIFFNESS FOR Regolith
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonRegolith=0.3 # POISSION'S RATIO FOR Regolith
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densRegolith=2650 # DENSITY FOR Regolith
densContainer=7850 # DENSITY FOR CONTAINER
numDamp=0.6
InitialPeneterationofExcavator=2.5e-3

#O.materials.append(FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=radians(0),density=densContainer,label='excavator'))
O.materials.append(FrictMat(young=youngRegolith,poisson=poissonRegolith,frictionAngle=radians(atRestFricDegree),density=densRegolith,label='regolith'))

ymport.textClumps("/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Penetration and Drag GRC-3 (Medium Dense)/Non-Periodic/Samples/128by64by64_n045_GRC-3_GCSC4p0_FinalBed_Relaxed_NonPeriodic_(4).txt",shift=Vector3(-0.001,0,0),material='regolith')
for b in O.bodies:
  if isinstance(b.shape,Sphere):
    if (b.state.pos[1]+b.shape.radius)>0.0054:
      O.bodies.erase(b.id)
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("maxY:",maxY)
mn,mx=aabbExtrema()
walls=aabbWalls([mn,mx],thickness=0.0001, oversizeFactor=1.5)
#for w in walls: w.shape.wire=False
O.bodies.append(walls[:3]+walls[4:]) #don't insert top wall



Excavator = O.bodies.append(ymport.stl('/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Penetration and Drag GRC-3 (Medium Dense)/Non-Periodic/External Objects/Bucket Excavator_HalfWidth_New.stl',wire=True,color=Vector3(1.0,1.0,0.0)))
ExcavatorID = [b for b in O.bodies if isinstance(b.shape,Facet)] # list of facets in simulation


O.engines = [VTKRecorder(fileName='/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Penetration and Drag GRC-3 (Medium Dense)/Non-Periodic/VTK/Output_N_Slow' ,recorders=['all'],iterPeriod=1000),
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Facet_Aabb()], label="collider"),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(eta=0.4,betan=normalDamp,betas=shearDamp,label='ContactModel')],
		[Law2_ScGeom_MindlinPhys_Mindlin(includeMoment=True,label='Mindlin')]
	),

	PyRunner(iterPeriod=1, command="penetration()" ,label='checker'),
	# construct CombinedKinematicEngine with label, add (with +) TranslationEngine and RotationEngine
	CombinedKinematicEngine(ids=Excavator,label='combEngine') + TranslationEngine(translationAxis=(0,-1,0),velocity=velocity) + RotationEngine(rotationAxis=(0,0,1), angularVelocity=0, rotateAroundZero=True, zeroPoint=(0.005171,0.008629,0.003200)),
	NewtonIntegrator(damping=numDamp,gravity=(0,gravityAcc,0)),
	PyRunner(iterPeriod=100,command='history()',label='recorder'),
]

O.dt=0.05*PWaveTimeStep()

# get TranslationEngine and RotationEngine from CombinedKinematicEngine
transEngine, rotEngine = combEngine.comb[0], combEngine.comb[1]
initialZeroPointVertical=rotEngine.zeroPoint[1]
print ("Zero Point Initial Vertical Position:",initialZeroPointVertical)


def penetration():
 transEngine.velocity = velocity
 print ("Zero Point Vertical Position:",rotEngine.zeroPoint[1])
 rotEngine.zeroPoint += Vector3(0,-1,0)*velocity*O.dt
 if rotEngine.zeroPoint[1]<=(initialZeroPointVertical-InitialPeneterationofExcavator):
 	checker.command='updateKinematicEngines()'

def updateKinematicEngines():
 transEngine.translationAxis=(1,0,0)	
 transEngine.velocity = velocity
 print ("Zero Point Horizontal Position:",rotEngine.zeroPoint[0])
 rotEngine.angularVelocity = 0
 rotEngine.zeroPoint += Vector3(1,0,0)*velocity*O.dt

def history():
  global Fx,Fy,Fz
  Fx=0
  Fy=0
  Fz=0
  for b in ExcavatorID:
    Fx+=O.forces.f(b.id,sync=True)[0]
    Fy+=O.forces.f(b.id,sync=True)[1]
    Fz+=O.forces.f(b.id,sync=True)[2]
  Dx=rotEngine.zeroPoint[0]
  Dy=rotEngine.zeroPoint[1]
  Dz=rotEngine.zeroPoint[2]
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,'Dx':Dx,'Dy':Dy,'Dz':Dz,})
  ## In that case we can still save the data to a text file at the the end of the simulation, with:
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Penetration and Drag GRC-3 (Medium Dense)/Non-Periodic/Numerical Output/penetrationandDrag_OutputData_128by50by64_n044_GRC-3_N_Slow.txt')
  
## declare what is to plot. "None" is for separating y and y2 axis
plot.plots={'i':('Fx','Fy','Fz','Dx','Dy','Dz')}
plot.plot()

