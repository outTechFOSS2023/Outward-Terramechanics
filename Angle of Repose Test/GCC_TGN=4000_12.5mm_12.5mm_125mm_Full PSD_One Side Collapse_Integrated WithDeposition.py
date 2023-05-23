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
## EFFECT OF GRAVITY. SAMPLE PREPARATION IS INTEGRATED. THE COLLAPSE OCCURES IN +X 
## DIRECTION BY THE MOVEMENT OF THE RIGHT BOUNDARY (+X) IN THE HEIGHT DIRECTION.
## PSD: FULL PSD
## SAMPLE LENGTH: 12.5mm
## SAMPLE WIDTH: 12.5mm
## COLLAPSE SYSTEM: ONE SIDE TRANSLATION IN THE HEIGHT DIRECTION
## POROSITY: THIS GCC IS INTEGRATED WITH SAMPLE PREPARATION (BY GRAVITY DEPOSITION). THE POROSITY OF THE COLUMN
## IS TO BE CALCULATED AFTER THE DEPOSITION USING SOLID VOLUME OF THE GRAINS (=9.708973078531903e-06) AND TOTAL VOLUME 
## OF THE COLUMN (=SAMPLE LENGTH (12.5MM)*SAMPLE WIDTH (12.5MM)*SAMPLE HEIGHT AFTER DEPOSITION (SEE BELOW--gravityDeposition)
## LENGTH TO HEIGHT SCALE: 1:10 (DESIRED)
## PERIODIC VS NON-PERIODIC VS THIN SLICE: NON-PERIODIC
## GRAIN SIZE LIMIT: 90-800 micrometer
## ANGLE OF REPOSE CALCULATION: YES
## TGN: 4000
## NOTE: THIS ANLYSIS IS TIME CONSUMING BECAUSE OF THE INCLUSION OF SAMPLE PREPARATION
#*************************************************************************/
from yade import ymport
from yade import export
from yade import pack
import shutil
import time
import collections
import numpy
import os
O.timingEnabled=True
#####MATERIAL#######
regolith=O.materials.append(FrictMat(young=1e8,poisson=0.25,frictionAngle=radians(35),density=3100))
container=O.materials.append(FrictMat(young=1e8,poisson=0.30,frictionAngle=radians(5),density=7850))
runway=O.materials.append(FrictMat(young=1e8,poisson=0.30,frictionAngle=radians(15),density=7850))
#####GEOMETRY########
#####CONTAINER#######
mn_container,mx_container=Vector3(0,0,0),Vector3(12.5e-3,12.5e-3,1.7) # CORNERS OF THE INITIAL PACKING. NOTE THAT FOR DIFFERENT CONFIGURATIONS, DIFFERENT HEIGHTS MUSY BE USED
container=aabbWalls([mn_container,mx_container],thickness=1.0e-3,oversizeFactor=1,material=container)
containerIds=O.bodies.append(container)
O.bodies.erase(containerIds[0])
O.bodies.erase(containerIds[2])
O.bodies.erase(containerIds[3])
O.bodies.erase(containerIds[4])
O.bodies.erase(containerIds[5])
######RUNWAY########
mn_runway,mx_runway=Vector3(0,0,0),Vector3((0.15+12.5e-3),12.5e-3,1.7) # CORNERS OF THE INITIAL PACKING. NOTE THAT FOR DIFFERENT CONFIGURATIONS, DIFFERENT HEIGHTS MUSY BE USED
runway=aabbWalls([mn_runway,mx_runway],thickness=1.0e-3,oversizeFactor=1,material=runway)
runwayIds=O.bodies.append(runway)
O.bodies.erase(runwayIds[5])
#####GRAINS########
yade.ymport.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Samples/GCC_TGN=4000_OriginalGrains_8.0mm*8.0mm_Combined Clumps and Spheres.txt",material=regolith,shift=Vector3(2.25e-3,2.25e-3,0))
maxX=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxY=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minX=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minY=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
minZ=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
print ("maxX: ",maxX,"maxY: ",maxY,"maxZ: ",maxZ,"minX: ",minX,"minY: ",minY,"minZ: ",minZ)
print('Soild Volume Clculated by YADE:',getSpheresVolume())
###GENERAL ENGINES###
newton=NewtonIntegrator(damping=0.5,gravity=(0,0,-9.81),dead=True)
O.engines=[
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb(),Bo1_Wall_Aabb(),Bo1_Box_Aabb()]),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom(),Ig2_Wall_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.3,betas=0.3,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
),
newton,
GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=10,timestepSafetyCoefficient=0.4),
PyRunner(command='calmPlate()',iterPeriod=2000,label='checker'),
TranslationEngine(translationAxis=[0,0,1],velocity=200,ids=containerIds,label='Trans',dead=True),
]
#####FUNCTIONS#######
def calmPlate():
	maxZ=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]) 
	O.bodies.append(utils.box(center=(6.25e-3,6.25e-3,(maxZ+1.0e-3)),extents=(6.25e-3,6.25e-3,1.0e-3),fixed=True,color=(0.2,0.5,1),wire=False))
	global topPlate
	topPlate=O.bodies[-1]
	O.engines=O.engines+[PyRunner(iterPeriod=10000,command='calm()',label='calmRunner1')]#No.8
	O.engines=O.engines+[TranslationEngine(ids=[topPlate.id],translationAxis=(0.,0.,-1.),velocity=0.2,label='accelerate',dead=True)]#No.8
	checker.command='calm()'

def calm():
	O.trackEnergy=True
	Ek=utils.kineticEnergy()
	print ("Kinetic Energy: ",Ek)
	O.engines[7].dead=False
	if Ek<1e-10:
		O.engines[7].dead=True
		O.engines[8].dead=False
		newton.dead=False
		newton.gravity=(0,0,-9.81)
		checker.command='gravityDeposition()'

def gravityDeposition():
	UNB=unbalancedForce()
	maxZ_depos=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	print ("maxZ_depos: ",maxZ_depos,"unbalanced force:",UNB)
	if maxZ_depos<0.2:
		O.engines[8].dead=True
	if O.iter>410000 and unbalancedForce()<8e-5:
		newton.dead=True
		maxZ_depos=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		minZ_depos=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		sampleHeightAfterDeposition=abs((maxZ_depos)-minZ_depos)
		sampleLengthAfterDeposition=12.5e-3
		sampleWidthAfterDeposition=12.5e-3
		V_T=sampleHeightAfterDeposition*sampleLengthAfterDeposition*sampleWidthAfterDeposition
		V_S=getSpheresVolume()
		GCPorosity=1-(V_S/V_T)
		export.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Samples/GCC_TGN=4000_12.5mm*12.5mm*Height_Full PSD_Integrated With Collapse_Deposited.txt")
		f= open("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Timing_Height_Unbalanced Force Info/GCC_TGN=4000_12.5mm*12.5mm*Height_Full PSD_Integrated With Collapse_Deposition Info.txt", 'w+')
		f.write('%s %s %s %s %s' % ("RealTime","SimTime","TimeStep","DeposHeight","UnbalancedForce"))
		f.write('\n')
		f.write('%f %f %f %f %f' % (O.realtime,O.time,O.iter,maxZ_depos,UNB))
		f.write('\n')
		f.close()
		for b in O.bodies:
			if b.isClumpMember:
				if b.state.pos[2]>0.125:
					O.bodies.erase(b.id)
		newton.dead=False
		newton.gravity=(0,0,-9.81)
		Trans.dead = False
		O.engines=O.engines+[VTKRecorder(fileName='/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/VTK Collection/GCC_TGN=4000_12.5mm*12.5mm*125mm_Full PSD_Integrated With Deposition_One Side Collapse' ,label='VTKRecorder',recorders=['spheres','facets','boxes'],iterPeriod=10000,dead=True)]
		checker.command='liftTheContainer()'
def liftTheContainer():
	O.engines[9].dead=False
	MinX_collapse=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxX_collapse=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MinY_collapse=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxY_collapse=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxZ_collapse=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	UNB=unbalancedForce()
	O.trackEnergy=True
	Ek=utils.kineticEnergy()
	print ("Unbalanced Force: ",UNB, "MaxZ During the Collapse: ",MaxZ_collapse, "MinX During the Collapse: ",MinX_collapse, "MaxX During the Collapse: ",MaxX_collapse, "MinY During the Collapse: ",MinY_collapse, "MaxY During the Collapse: ",MaxY_collapse)
	if O.iter>1610000 and unbalancedForce()<0.8e-5:
	# Eliminate far grains
		for b in O.bodies:
			if isinstance(b.shape,Sphere):
				X = abs(b.state.pos[0])
				if (X+b.shape.radius)>=0.08:
					O.bodies.erase(b.id)				
			# Slope Calculation
		MaxZ_endofcollapse=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		for b in O.bodies: 
			if isinstance(b.shape,Sphere):
				if b.state.pos[2]+b.shape.radius == MaxZ_endofcollapse:
					CenterX = b.state.pos[0]
					CenterY = b.state.pos[1]
					CenterZ = b.state.pos[2]+b.shape.radius
		n=10
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

			MaxXBottom = abs(max(Bottom[:,0]))

			MaxXTop = abs(max(Top[:,0]))

			angle_MaxX = math.atan(Height/(MaxXBottom-MaxXTop))

			angles = numpy.array([angle_MaxX])
			allAngles.append(angles)

		allAngles = numpy.array(allAngles).reshape(-1,1)
		meanAngle = numpy.mean(allAngles, axis=0).reshape(-1,1)
		angleRepose = numpy.concatenate((allAngles, meanAngle), axis=0).reshape(-1,1)

		with open(('/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Automatic Angle of Repose (AOR) Calculations/'+'GCC_TGN=4000_12.5mm*12.5mm*125mm_Full PSD_With Deposition_One Side Collapse'+'.dat'),"w") as f:   
			f.write('# angle_MaxX \n')
			for k in range(0,len(angleRepose)):
				f.write('\n')
				f.write("\n".join(" ".join(map(str, k)) for k in angleRepose[k]\
					.reshape(-1,1)))
		maxZ_collapse=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		maxX_collapse=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		O.save('/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Saved Simulations and Final Configurations/GCC_TGN=4000_12.5mm*12.5mm*125mm_Full PSD_Integrated With Deposition_One Side Collapse.yade.bz2')
		export.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Saved Simulations and Final Configurations/GCC_TGN=4000_12.5mm*12.5mm*125mm_Full PSD_Integrated With Deposition_One Side Collapse.txt")			
		f= open("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Angle of Repose Test/Timing_Height_Unbalanced Force Info/GCC_TGN=4000_12.5mm*12.5mm*125mm_Full PSD_One Side Collapse_Integrated With Deposition_Collapse Info.txt", 'w+')
		f.write('%s %s %s %s %s %s' % ("RealTime","SimTime","TimeStep","CollapseHeight","CollapseRunout","UnbalancedForce"))
		f.write('\n')
		f.write('%f %f %f %f %f %f' % (O.realtime,O.time,O.iter,maxZ_collapse,maxX_collapse,UNB))
		f.write('\n')
		f.close()
		O.pause()


