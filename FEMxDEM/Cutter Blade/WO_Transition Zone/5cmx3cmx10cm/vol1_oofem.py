#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

######################################################################
# This Python script is the OOFEM (FEM) part of a FEMxDEM simulation of a cutter blade (CB) test.
# The DEM part of the simulation is defined in a diferent Python script entitled vol1_yade.py.
# There are also some other dependencies which link the DEM and FEM for such adjustments of  as output request and number of iterations.
# NOTE: ALL EXTERNAL FILES (DEM GEOMETRY), AND OUTPUT FOLDERES ARE PROVIDED. THE CORRECT ADDRESSES
# MUST BE PROVIDED BY THE USER FROM THEIR OWN LOCAL MACHINES. 
# The FEM model has the following chracteristics:
# 1- No microparameter calibrations were conducted and the model has been setup to show the runability of  FEMxDEM.  A simple  isotropic linear  elastic 
# model (isoLE) has been used. isoLE needs  material density d (kg/m3), material Young modulus E (pa) , material Poisson ratio n (unitless), and 
# material thermal dilatation coefficient tAlpha  (unitless). However, OOFEM is a highly extendable  and new constituitive theory can be used to model 
# a more realistic response of the FEM part. 
# 2- FEM  size is 5cm(X-direction)*3.0cm(Y-direction)x10cm(Z-Direction) covered by a 3mm thick clump GRC-3+LMA-1
# DEM particles. 
# 3- number of elements in each direction  is adjustable by the user depending on tthe desired precision. (Here: nx, ny, nz = 5, 4, 15)
# 4- Surface coupling (SC) FEM*DEM has been implemented. SC has been shown to run faster than volumetic coupling (VC). 
# 5- Wall boundaries and the cutter blade has been generated using YADE built-in geometry generation functions in vol1_yade.py for facets and therefore are not deformable. 
# However, it is possible to extend the FEM to all external objects.
# 6-  Please do not change other sections to avoid malfunction.
######################################################################
import liboofem
import numpy
import itertools

class Node:
	def __init__(self,id,coords):
		self.id = id
		self.coords = list(coords)
		self.bcs = [0,0,0]
	def toOofem(self):
		kw = dict(coords=self.coords)
		if any(self.bcs):
			kw["bc"] = self.bcs
		return liboofem.node(self.id,domain,**kw)

class Brick:
	def __init__(self,id,nodes):
		self.id = id
		self.nodes = list(nodes)
	def toOofem(self):
		return liboofem.element("LSpace",self.id, domain, nodes=self.nodes, mat=1, crossSect=1, nlgeo=1)

# THIS IS THE FEM MESH SIZE. THE DIMENSIONS IN THE X AND Z DIRECTIONS ARE EXACTLY
# THOSE CONSIDERED FOR THE DEM BUT THE HEIGHT CAN BE CHANGED. HOWEVER, WITH
# ANY CHANGES IN THE DIMENSION, THE FACET BOUNDARIES IN THE DEM MUST BE ADJUSTED.
dx, dy, dz = 0.05, 0.03, 0.1 # dimensions
# CHANGE THE NUMBER OF ELEMENTS IN EACH DIRECTION IF DESIRED
nx, ny, nz = 5, 4, 15 # number of elements per dimension
ymin=0.03

def createOofemNodes():
	xs, ys, zs = [[i*d/n for i in xrange(n+1)] for d,n in zip((dx,dy,dz),(nx,ny,nz))]
	ys=[y-ymin for y in ys]

	nodes = [Node(-1,(x,y,z)) for z in zs for y in ys for x in xs]
	for i,n in enumerate(nodes):
		n.id = i+1
		if n.coords[0] == 0. or n.coords[0] == dx or n.coords[1] == -1*ymin or n.coords[2] == 0. or n.coords[2] == dz :
			n.bcs = [1,1,1]
	return [n.toOofem() for n in nodes]

def createOofemElems():
	elems = []
	nx1 = nx + 1
	nx1ny1 = (nx+1)*(ny+1)
	for ix in xrange(nx):
		for iy in xrange(ny):
			for iz in xrange(nz):				
				n5 = 1 + ix + nx1*iy + nx1ny1*iz
				n8 = n5 + 1
				n6 = n5 + nx1
				n7 = n6 + 1
				n1 = n5 + nx1ny1
				n2 = n6 + nx1ny1
				n3 = n7 + nx1ny1
				n4 = n8 + nx1ny1
				elems.append(Brick(-1,(n1,n2,n3,n4,n5,n6,n7,n8)))
	for i,e in enumerate(elems):
		e.id = i+1
	return [e.toOofem() for e in elems]


############################################################


############################################################

# engngModel
problem = liboofem.engngModel("nldeidynamic",1,nSteps=5,dumpCoef=.1,deltaT=10000,outFile="/tmp/vol1_oofem.out")

# domain
domain = liboofem.domain(1, 1, problem, liboofem.domainType._3dMode, tstep_step=10000, dofman_all=True, element_all=True)
problem.setDomain(1, domain, True)

# vtkxml
vtkxmlModule = liboofem.vtkxml(1,problem,tstep_step=10000,vars=[1,4],primvars=[1])
exportModuleManager = problem.giveExportModuleManager()
exportModuleManager.resizeModules(1)
exportModuleManager.setModule(1,vtkxmlModule)

# boundary condition and time function
ltf = liboofem.loadTimeFunction("constantFunction",1,f_t=0)
bc = liboofem.boundaryCondition(1, domain, loadTimeFunction=1, prescribedValue=0.0)

# nodes
nodes = createOofemNodes()

# material and cross section
mat = liboofem.isoLE(1, domain, d=1550, E=0.7e8, n=0.3, tAlpha=0)
cs  = liboofem.simpleCS(1, domain)

# elements
elems = createOofemElems()

# add eveything to domain (resize container first)
domain.resizeDofManagers(len(nodes))
for n in nodes:
	domain.setDofManager(n.number, n)
domain.resizeElements(len(elems))
for e in elems:
	domain.setElement(e.number, e)
domain.resizeMaterials(1)
domain.setMaterial(1, mat)
domain.resizeCrossSectionModels(1)
domain.setCrossSection(1, cs)
domain.resizeBoundaryConditions(1)
domain.setBoundaryCondition(1, bc)
domain.resizeFunctions(1)
domain.setFunction(ltf.number, ltf)

vtkxmlModule.initialize() # (!)
############################################################

def vtkExport(i):
	problem.giveExportModuleManager().giveModule(1).doForcedOutput(problem.giveCurrentStep())

if __name__ == "__main__":
	problem.checkProblemConsistency();
	problem.init();
	problem.postInitialize();
	problem.setRenumberFlag();
	problem.solveYourself();
	problem.terminateAnalysis();
