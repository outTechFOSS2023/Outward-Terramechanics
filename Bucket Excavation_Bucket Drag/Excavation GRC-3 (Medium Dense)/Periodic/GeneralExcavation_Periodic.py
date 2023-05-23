# -*- coding: utf-8 -*-
#*********************************************************************************************************
#*********************************************************************************************************
#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.#*********************************************************************************************************
#*********************************************************************************************************                                                                     
#  THIS PROGRAM IS A FREE SOFTWARE; IT IS LICENCED UNDER THE TERMS OF THE 
#  GNU GENERAL PUBLIC LICENCE V2 OR LATER
#*************************************************************************/
#***************************************************************************************************************
## THIS SCRIPT IS An INPUT FILE DETAILING THE SIMULATION OF PERIODIC  EXCAVATION OF GRC-3 USING A BUCKET 
## EXCAVATOR IN MD STATE (n~0.44) USING YADE. THE ROTATIONAL EXCAVATION STARTS WITHOUT PEE-PENETRATION OF THE 
## BUCKET INTO THE BED. THE DESIRED ANGULAR VELOCITY CAN BE ADJUSTED IN LINE 29 BUT FOR STABILITY PURPOSES, 
## LOWR VALUES<1 RAD/S ARE RECOMMENDED. THE USED CONSTITUIIVE LAW IS HERTZ-MINDLIN WITH SOME ROLLING RESISTANCE
## IMPLEMENTED FOR BETTER STABILITY (LINE 80), SINCE SOME SPHERICAL PARTICLES ARE PRESENT IN THE MODEL. 
## FORCE HISTORY IN ALL THREE X, Y, AND Z ARE RECORDED AND VTK OF THE DEMOFRMATION IS STORED AS WELL. A COMBINED 
## KINEMATIC ENGINE CONTROLS THE INITIAL VERTICAL MOVEMENT (WHICH GET THE BUCKET CLOSER TO THE REGOLITH SURFACE)
## AND ROTATIONAL EXCAVATION. 
#***************************************************************************************************************
from yade import pack
from yade import plot
from yade import ymport
from yade import export
############################################
###           INPUT PARAMETERS           ###
############################################
velocity=50 #m/s THAT SHOULD BE 0.01 m/s
angularVelocity=1 #RAD/S
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
InitialPeneterationofExcavator=0.40e-3


O.periodic=True
O.cell.setBox(0.012814050300000001,0.015,0.0064090263)


#O.materials.append(FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=radians(0),density=densContainer,label='excavator'))
O.materials.append(FrictMat(young=youngRegolith,poisson=poissonRegolith,frictionAngle=radians(atRestFricDegree),density=densRegolith,label='regolith'))
mn,mx=Vector3(0,0.001,0),Vector3(0.012814050300000001,(0.0063999958+0.001),0.0064090263) # CORNERS OF THE TRIAXIAL CONTAINER
walls=aabbWalls([mn,mx],thickness=0.0001, oversizeFactor=1.5)
wallIds=O.bodies.append(walls)
O.bodies.erase(0)
O.bodies.erase(1)
O.bodies.erase(3)
O.bodies.erase(4)
O.bodies.erase(5)



ymport.textClumps("/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Excavation GRC-3 (Medium Dense)/Periodic/Samples/128by64by64_n045_GRC-3_GCSC4p0_FinalBed_Relaxed_NonPeriodic_resetToZero_(4).txt",shift=Vector3(0,0.001,0),material='regolith')
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("minX:",minX,"minY:",minY,"minZ:",minZ,"maxX:",maxX,"maxY:",maxY,"maxZ:",maxZ)


Excavator = O.bodies.append(ymport.stl('/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Excavation GRC-3 (Medium Dense)/Periodic/External Objects/Bucket Excavator_HalfWidth_Displaced.stl',wire=True,color=Vector3(0.9,0.0,0.4)))
ExcavatorID = [b for b in O.bodies if isinstance(b.shape,Facet)] # list of facets in simulation


O.engines = [VTKRecorder(fileName='/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Excavation GRC-3 (Medium Dense)/Periodic/VTK/Output_P' ,recorders=['all'],iterPeriod=1000),
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Facet_Aabb()],allowBiggerThanPeriod=True,label="collider"),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(eta=0.4,betan=normalDamp,betas=shearDamp,label='ContactModel')],
		[Law2_ScGeom_MindlinPhys_Mindlin(includeMoment=True,label='Mindlin')]
	),

	PyRunner(iterPeriod=1, command="penetration()" ,label='checker'),
	# construct CombinedKinematicEngine with label, add (with +) TranslationEngine and RotationEngine
	CombinedKinematicEngine(ids=Excavator,label='combEngine') + TranslationEngine(translationAxis=(0,-1,0),velocity=velocity) + RotationEngine(rotationAxis=(0,0,1), angularVelocity=0, rotateAroundZero=True, zeroPoint=( 0.006171,0.009629,0.003200)),
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
 #rotEngine.angularVelocity = angularVelocity
 rotEngine.zeroPoint += Vector3(0,-1,0)*velocity*O.dt
 if rotEngine.zeroPoint[1]<=(initialZeroPointVertical-InitialPeneterationofExcavator):
 	checker.command='updateKinematicEngines()'

def updateKinematicEngines():
 transEngine.translationAxis=(0,0,1)	
 transEngine.velocity = 0
 print ("Zero Point Vertical Position:",rotEngine.zeroPoint[1])
 rotEngine.angularVelocity = 1
 #rotEngine.zeroPoint += Vector3(0,0,1)*velocity*O.dt

def history():
  global Fx,Fy,Fz
  Fx=0
  Fy=0
  Fz=0
  for b in ExcavatorID:
    Fx+=O.forces.f(b.id,sync=True)[0]
    Fy+=O.forces.f(b.id,sync=True)[1]
    Fz+=O.forces.f(b.id,sync=True)[2]
  yade.plot.addData({'i':O.iter,'Fx':Fx,'Fy':Fy,'Fz':Fz,})
  ## In that case we can still save the data to a text file at the the end of the simulation, with:
  plot.saveDataTxt('/home/ngoudarzi/Desktop/Bucket Excavation_Bucket Drag/Excavation GRC-3 (Medium Dense)/Periodic/Numerical Output.txt')
  
## declare what is to plot. "None" is for separating y and y2 axis
plot.plots={'i':('Fx','Fy','Fz')}
plot.plot()

