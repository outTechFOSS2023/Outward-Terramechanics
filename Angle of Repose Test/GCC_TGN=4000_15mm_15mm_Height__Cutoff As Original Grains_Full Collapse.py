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
## THIS SCRIPT SIMULATES THE COLLAPSE OF A GRANULAR COLUMN OF GRC-3+LMA-1 UNDER THE 
## EFFECT OF GRAVITY. THE COLLAPSE OCCURES IN XY PLANE ON A CYLINDRICAL PLATFORM 
## DIRECTION BY THE MOVEMENT OF WHOLE ASSEMBLY CONTAINER IN THE HEIGHT DIRECTION.
## PSD: ORIGINAL WITH SINGE SPHERES FOR GRAINS SMALLER THAN THE CUTOFF (90 micrometer)
## SAMPLE LENGTH: 15mm
## SAMPLE WIDTH: 15mm
## SAMPLE HEIGHT: 149.4mm
## COLLAPSE SYSTEM: FULL TRANSLATION OF THE CONTAINER IN THE HEIGHT DIRECTION
## POROSITY: LOOSE (n=0.713963257)
## LENGTH TO HEIGHT SCALE: 1:10
## PERIODIC VS NON-PERIODIC VS THIN SLICE: NON-PERIODIC
## GRAIN SIZE LIMIT: 90-800 micrometer
## ANGLE OF REPOSE CALCULATION: YES
## TGN: 4000
## CUTOFF: LESS THAN 90 micrometers AS SINGLE SPHERES
#*********************************************************************************************************
from yade import utils
from yade import export
from yade import ymport
from yade import plot
from yade import qt
import time
import numpy
import math
#########################################################
#####################   MATERIALS   #####################
#########################################################
regolith=O.materials.append(FrictMat(young=1e8,poisson=0.25,frictionAngle=radians(35),density=2600))
plateform=O.materials.append(FrictMat(young=1e8,poisson=0.30,frictionAngle=radians(15),density=7850))
walls=O.materials.append(FrictMat(young=1e8,poisson=0.30,frictionAngle=radians(5),density=7850))

#########################################################
###################  TEST BOX   ####################
#########################################################
id_walls=O.bodies.append(geom.facetBox(center=(0,0,0.08), extents=(7.5e-3,7.5e-3,0.08), orientation=Quaternion((1, 0, 0), 0), wallMask=15,material=walls))

#########################################################
###################  Bed Platform   #####################
#########################################################
RBedPlatform = 45.0e-3
id_BedPlatform=O.bodies.append(geom.facetCylinder(center=(0,0,-0.0075),radius=RBedPlatform,height=0.015,\
  wallMask=7,orientation=Quaternion((1, 0, 0), 0),segmentsNumber=50,wire=False,color=(1,1,1),material=plateform))
#########################################################
#####################  PARTICLES   ######################
#########################################################
yade.ymport.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Samples/GCC_TGN=4000_11mm*11mm*1000mm_<Cutoff As Original Grains_Splitted Assembly_deposited.txt",material=regolith)
for b in O.bodies:
	if isinstance(b.shape,Sphere):
		if b.state.pos[0]>0.0075 or b.state.pos[0]<-0.0075 or b.state.pos[1]>0.0075 or b.state.pos[1]<-0.0075 or b.state.pos[2]>0.16 or b.state.pos[2]<0:
			O.bodies.erase(b.id)
MinX_import=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxX_import=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MinY_import=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxY_import=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MinZ_import=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxZ_import=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print('MinX_import:',MinX_import,'MaxX_import:',MaxX_import,'MinY_import:',MinY_import,'MaxY_import:',MaxY_import,'MinZ_import:',MinZ_import,'MaxZ_import:',MaxZ_import)
print('Soild Volume Clculated by YADE:',getSpheresVolume())
#########################################################
######################   ENGINES   ######################
#########################################################
newton=NewtonIntegrator(damping=0.5,gravity=(0,0,-9.81),dead=True)

O.engines=[VTKRecorder(fileName='/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/VTK Collection/GCC_TGN=4000_15mm*15mm*Height_<Cutoff As Original Grains_Full Collapse' ,recorders=['spheres','facets','boxes'],iterPeriod=10000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()]),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.30,betas=0.30,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
),
newton,
GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.4),
PyRunner(command='calmPlate()',iterPeriod=2000,label='checker'),
TranslationEngine(translationAxis=[0,0,1],velocity=200,ids=id_walls,label='Trans',dead=True),
]

def calmPlate():
	maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]) 
	O.bodies.append(utils.box(center=(0,0,(maxZ+5.0e-3)),extents=(7.5e-3,7.5e-3,5.0e-3),material=walls,fixed=True,color=(0.2,0.5,1),wire=False))
	global topPlate_calm
	topPlate_calm=O.bodies[-1]
	O.engines=O.engines+[PyRunner(iterPeriod=10000,command='calm()',label='calmRunner')]#No.8
	checker.command='calm()'

def calm():
	O.trackEnergy=True
	Ek=utils.kineticEnergy()
	print ("Kinetic Energy: ",Ek)
	O.engines[6].dead=False
	if Ek<1e-10:
		O.engines[6].dead=True
		#O.bodies.erase(101889)
		newton.dead=False
		newton.gravity=(0,0,-9.81)
		Trans.dead = False
		checker.command='liftTheContainer()'

def liftTheContainer():
	global RBedPlatform
	O.trackEnergy=True
	Ek=utils.kineticEnergy()
	MinX_collapse=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxX_collapse=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MinY_collapse=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxY_collapse=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxZ_collapse=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	# Eliminate far grains
	for b in O.bodies:
		if isinstance(b.shape,Sphere):
			X = abs(b.state.pos[0])
			Y = abs(b.state.pos[1])
			Z = b.state.pos[2]
			if (((X+b.shape.radius)*(X+b.shape.radius)+(Y+b.shape.radius)*(Y+b.shape.radius))>=(RBedPlatform*RBedPlatform)):
				O.bodies.erase(b.id)
	UNB=unbalancedForce()
	O.trackEnergy=True
	Ek=utils.kineticEnergy()
	print ("Kinetic Energy: ",Ek, "Unbalanced Force: ",UNB, "MaxZ During the Collapse: ",MaxZ_collapse, "MinX During the Collapse: ",MinX_collapse, "MaxX During the Collapse: ",MaxX_collapse, "MinY During the Collapse: ",MinY_collapse, "MaxY During the Collapse: ",MaxY_collapse)
	if O.iter>1100000 and unbalancedForce()<1e-3:
	# Eliminate far grains
		for b in O.bodies:
			if isinstance(b.shape,Sphere):
				X = abs(b.state.pos[0])
				Y = abs(b.state.pos[1])
				Z = b.state.pos[2]
				if (((X+b.shape.radius)*(X+b.shape.radius)+(Y+b.shape.radius)*(Y+b.shape.radius))>=(RBedPlatform*RBedPlatform)) or \
				(Z-b.shape.radius <= 0.0) or (Z+b.shape.radius >= 0.1):
					O.bodies.erase(b.id)	
			# Slope Calculation
		MaxZ_endofcollapse=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		for b in O.bodies: 
			if isinstance(b.shape,Sphere):
				if b.state.pos[2]+b.shape.radius == MaxZ_endofcollapse:
					CenterX = b.state.pos[0]
					CenterY = b.state.pos[1]
					CenterZ = b.state.pos[2]+b.shape.radius
		n=2
		allAngles =[]
		for i in range(n):
			Bottom =[]
			Top =[]
			Height = MaxZ_endofcollapse/n
			for b in O.bodies:
				if isinstance(b.shape,Sphere):
					if i*Height-2*b.shape.radius <= b.state.pos[2] <= i*Height+2*b.shape.radius:
						BottomPos = numpy.array([b.state.pos[0], b.state.pos[1]])
						Bottom.append(BottomPos)

					if (i+1)*Height-2*b.shape.radius <= b.state.pos[2] <= (i+1)*Height+2*b.shape.radius:
						TopPos = numpy.array([b.state.pos[0], b.state.pos[1]])
						Top.append(TopPos)

			Bottom = numpy.array(Bottom).reshape(-1,2)
			Top = numpy.array(Top).reshape(-1,2)

			MinXBottom = abs(min(Bottom[:,0]))
			MaxXBottom = abs(max(Bottom[:,0]))
			MinYBottom = abs(min(Bottom[:,1]))
			MaxYBottom = abs(max(Bottom[:,1]))

			MinXTop = abs(min(Top[:,0]))
			MaxXTop = abs(max(Top[:,0]))
			MinYTop = abs(min(Top[:,1]))
			MaxYTop = abs(max(Top[:,1]))

			angle_MaxX = math.atan(Height/(MaxXBottom-MaxXTop))
			angle_MinX = math.atan(Height/(MinXBottom-MinXTop))
			angle_MaxY = math.atan(Height/(MaxYBottom-MaxYTop))
			angle_MinY = math.atan(Height/(MinYBottom-MinYTop))

			angles = numpy.array([angle_MaxX,angle_MinX,angle_MaxY,angle_MinY])
			allAngles.append(angles)

		allAngles = numpy.array(allAngles).reshape(-1,4)
		meanAngle = numpy.mean(allAngles, axis=0).reshape(-1,4)
		angleRepose = numpy.concatenate((allAngles, meanAngle), axis=0).reshape(-1,4)

		with open(('/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Automatic Angle of Repose (AOR) Calculations/'+'GCC_TGN=4000_15mm*15mm*Height_<Cutoff As Original Grains_Full Collapse'+'.dat'),"w") as f:   
			f.write('# angle_MaxX  angle_MinX  angle_MaxY  angle_MinY\n')
			for k in range(0,len(angleRepose)):
				f.write('\n')
				f.write("\n".join(" ".join(map(str, k)) for k in angleRepose[k]\
					.reshape(-1,4)))
		O.save('/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Saved Simulations and Final Configurations/GCC_TGN=4000_15mm*15mm*Height_<Cutoff As Original Grains_Full Collapse.yade.bz2')
		export.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Saved Simulations and Final Configurations/GCC_TGN=4000_15mm*15mm*Height_<Cutoff As Original Grains_Full Collapse.txt")			
		O.pause()