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
## EFFECT OF GRAVITY. THE COLLAPSE OCCURES IN +X DIRECTION ONLY BY ROTATING
## THE RIGHT WALL FROM A PIVOT POINT ON TOP OF THE CONTAINER.
## THIN SLICES (IN Z DIRECTION) CAN BE COLLAPSED WITH THIS SCRIPT 
## PSD: HALF PSD
## SAMPLE LENGTH : 2.524mm
## SAMPLE WIDTH : 0.6mm (THIN SLICE)
## SAMPLE HEIGHT: 10.091mm
## COLLAPSE SYSTEM: ONE SIDE 
## POROSITY: LOOSE (n=0.51898238)
## NUMERICAL TO EXPERIMENTAL SCALE: 1:15
## PERIODIC VS NON-PERIODIC VS THIN SLICE: NON PERIODIC/THIN SLICE
## GRAIN SIZE LIMIT: 90-800mm
## ANGLE OF REPOSE CALCULATED: YES
#*************************************************************************/
from yade.pack import *
from yade import utils
from yade import export
from yade import ymport
from yade import plot
import time
import numpy
import math
############################################
##########   DEFINING VARIABLES    #########
############################################
regolithFricDegree= 35 # INITIAL CONTACT FRICTION FOR SAMPLE PREPARATION
containFricDegree=0
runwayFrictDegree=5
normalDamp=0.3 # NUMERICAL DAMPING
shearDamp=0.3
youngRegolith=1.0e8# CONTACT STIFFNESS FOR REGOLITH
youngContainer=210e9 # CONTACT STIFFNESS FOR CONTAINER
poissonRegolith=0.3 # POISSION'S RATIO FOR REGOLITH
poissionContainer=0.25 # POISSION'S RATIO FOR CONTAINER
densRegolith=3100 # DENSITY FOR REGOLITH
densContainer=7850 # DENSITY FOR CONTAINER
gravAcc=(0,-9.810,0) # NO GRAVITY IS APPLIED TO PREVENT SEGRAGATION
numDamp=0.5
thick=0.0001
mn,mx=Vector3(0,0,0),Vector3(2.524e-3,0.015,2.524e-3) # CORNER OF DEPOSITION BOX

#########################################################
#####################   MATERIALS   #####################
#########################################################
O.materials.append(FrictMat(young=youngRegolith,poisson=poissonRegolith,frictionAngle=radians(regolithFricDegree),density=densRegolith,label='regolith'))
O.materials.append(FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=radians(containFricDegree),density=densContainer,label='walls'))
O.materials.append(FrictMat(young=youngContainer,poisson=poissionContainer,frictionAngle=radians(runwayFrictDegree),density=densContainer,label='runWay'))
mnContain,mxContain=Vector3(0,0,0),Vector3(2.524e-3,0.015,0.6e-3) # CORNER OF THE TEST CONTAINER
mnRunway,mxRunway=Vector3(0,0,0),Vector3(0.025,0.015,0.6e-3) # CORNERS OF THE RUNWAY

#########################################################
###################      TEST BOX    ####################
#########################################################

mainBox=aabbWalls([mnContain,mxContain],thickness=thick,material='walls',oversizeFactor=1.0,color=(1,1,1),wire=False)
id_walls=O.bodies.append(mainBox)
O.bodies.erase(0)
O.bodies.erase(3)
O.bodies.erase(2)
O.bodies.erase(4)
O.bodies.erase(5)


#########################################################
###################  Bed Platform   #####################
#########################################################
runWay=aabbWalls([mnRunway,mxRunway],thickness=thick,material='runWay',oversizeFactor=1.0,color=(1,1,1),wire=False)
id_runWay=O.bodies.append(runWay)


#########################################################
#####################  PARTICLES   ######################
#########################################################
yade.ymport.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Thin Slice GCC (GRC-3plusLMA-1)/Samples/GCC_2p524byFREEby2p524_nMDD_GRC-3PlusLMA-1_Trimmed_Column_HalfPSD.txt",shift=Vector3(0,0,0),material='regolith',color=(0.6,0.3,0.4))
for b in O.bodies:
	if isinstance(b.shape,Sphere):
		if (b.state.pos[2]+b.shape.radius)>=0.6e-3:
			O.bodies.erase (b.id)
MinX_import=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxX_import=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MinY_import=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxY_import=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MinZ_import=min([b.state.pos[2]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
MaxZ_import=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
SampleLength=(MaxX_import-MinX_import)
SampleHeight=(MaxY_import-MinY_import)
SampleWidth=(MaxZ_import-MinZ_import)
print('MinX_import:',MinX_import,'MaxX_import:',MaxX_import,'MinY_import:',MinY_import,'MaxY_import:',MaxY_import,'MinZ_import:',MinZ_import,'MaxZ_import:',MaxZ_import)
print('SampleLength:',SampleLength,'SampleHeight:',SampleHeight,'SampleWidth:',SampleWidth)
V_solid=getSpheresVolume()
print('Soild Volume Clculated by YADE:',V_solid)
#########################################################
######################   ENGINES   ######################
#########################################################
newton=NewtonIntegrator(damping=numDamp,gravity=gravAcc,dead=False)

O.engines=[VTKRecorder(fileName='/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Thin Slice GCC (GRC-3plusLMA-1)/VTK/Thin Slice GCC' ,recorders=['all'],iterPeriod=5000),
ForceResetter(),
InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()],allowBiggerThanPeriod=True),
InteractionLoop(
[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
[Ip2_FrictMat_FrictMat_MindlinPhys(betan=normalDamp,betas=shearDamp,label='ContactModel')],
[Law2_ScGeom_MindlinPhys_Mindlin(label='Mindlin')]
),
newton,
#GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.55),
PyRunner(command='rotateTheRightWall()',iterPeriod=300,label='checker'),
#PyRunner(command='liftTheContainer()',iterPeriod=2000,label='checker'),
RotationEngine(angularVelocity=5000,rotationAxis=(0,0,1),rotateAroundZero=True,zeroPoint=(2.524e-3,0.015,0),ids=id_walls,label='Rot',dead=True),
]
#O.dt=1.0e-7
O.dt=0.1*PWaveTimeStep()


def rotateTheRightWall():
	Rot.dead=False
	global RBedPlatform
	O.trackEnergy=True
	Ek=utils.kineticEnergy()
	MinX_collapse=min([b.state.pos[0]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxX_collapse=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MinY_collapse=min([b.state.pos[1]-b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxY_collapse=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	MaxZ_collapse=max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
	UNB=unbalancedForce()
	O.trackEnergy=True
	Ek=utils.kineticEnergy()
	if O.iter>5000:
		Rot.dead=True
		Rot.angularVelocity=0.0
		O.bodies.erase(1)
	print ("Kinetic Energy: ",Ek, "Unbalanced Force: ",UNB, "MaxZ During the Collapse: ",MaxZ_collapse, "MinX During the Collapse: ",MinX_collapse, "MaxX During the Collapse: ",MaxX_collapse, "MinY During the Collapse: ",MinY_collapse, "MaxY During the Collapse: ",MaxY_collapse)
	if O.iter>2000000 and unbalancedForce()<0.001:
	# Eliminate far grains
		for b in O.bodies:
			if isinstance(b.shape,Sphere):
				X = abs(b.state.pos[0])
				if (X+b.shape.radius)>=0.02:
					O.bodies.erase(b.id)				
			# Slope Calculation
		MaxY_endofcollapse=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		for b in O.bodies: 
			if isinstance(b.shape,Sphere):
				if b.state.pos[1]+b.shape.radius == MaxY_endofcollapse:
					CenterX = b.state.pos[0]
					CenterY = b.state.pos[1]+b.shape.radius
					CenterZ = b.state.pos[2]
		n=200
		allAngles =[]
		for i in range(n):
			Bottom =[]
			Top =[]
			Height = MaxY_endofcollapse/n
			for b in O.bodies:
				if isinstance(b.shape,Sphere):
					if i*Height-2*b.shape.radius <= b.state.pos[1] <= i*Height+2*b.shape.radius:
						BottomPos = numpy.array([b.state.pos[0], b.state.pos[2]])
						Bottom.append(BottomPos)

					if (i+1)*Height-2*b.shape.radius <= b.state.pos[1] <= (i+1)*Height+2*b.shape.radius:
						TopPos = numpy.array([b.state.pos[0], b.state.pos[2]])
						Top.append(TopPos)

			Bottom = numpy.array(Bottom).reshape(-1,2)
			Top = numpy.array(Top).reshape(-1,2)

			MinXBottom = abs(min(Bottom[:,0]))
			MaxXBottom = abs(max(Bottom[:,0]))
			MinZBottom = abs(min(Bottom[:,2]))
			MaxZBottom = abs(max(Bottom[:,2]))

			MinXTop = abs(min(Top[:,0]))
			MaxXTop = abs(max(Top[:,0]))
			MinZTop = abs(min(Top[:,2]))
			MaxZTop = abs(max(Top[:,2]))

			angle_MaxX = math.atan(Height/(MaxXBottom-MaxXTop))
			angle_MinX = math.atan(Height/(MinXBottom-MinXTop))
			angle_MaxZ = math.atan(Height/(MaxZBottom-MaxZTop))
			angle_MinZ = math.atan(Height/(MinZBottom-MinZTop))

			angles = numpy.array([angle_MaxX,angle_MinX,angle_MaxZ,angle_MinZ])
			allAngles.append(angles)

		allAngles = numpy.array(allAngles).reshape(-1,4)
		meanAngle = numpy.mean(allAngles, axis=0).reshape(-1,4)
		angleRepose = numpy.concatenate((allAngles, meanAngle), axis=0).reshape(-1,4)

		with open(('/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Thin Slice GCC (GRC-3plusLMA-1)/Numerical Output/'+'Output'+'.dat'),"w") as f:   
			f.write('# angle_MaxX  angle_MinX  angle_MaxZ  angle_MinZ\n')
			for k in range(0,len(angleRepose)):
				f.write('\n')
				f.write("\n".join(" ".join(map(str, k)) for k in angleRepose[k]\
					.reshape(-1,4)))
		maxY_collapse=max([b.state.pos[1]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		maxX_collapse=max([b.state.pos[0]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
		O.save('/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Thin Slice GCC (GRC-3plusLMA-1)/Numerical Output/Output.yade.bz2')
		export.textClumps("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Thin Slice GCC (GRC-3plusLMA-1)/Numerical Output/Output_Particles.txt")			
		f= open("/home/ngoudarzi/Desktop/Granular Column Collapse (GCC)/Thin Slice GCC (GRC-3plusLMA-1)/Numerical Output/Output_Info.txt", 'w+')
		f.write('%s %s %s %s %s %s' % ("RealTime","SimTime","TimeStep","CollapseHeight","CollapseRunout","UnbalancedForce"))
		f.write('\n')
		f.write('%f %f %f %f %f %f' % (O.realtime,O.time,O.iter,maxY_collapse,maxX_collapse,UNB))
		f.write('\n')
		f.close()
		O.pause()

