#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

######################################################################
# This Python script is the coupling script a FEMxDEM simulation of a cutter blade (CB) test.
# The FEM and DEM parts of the simulation is defined in a separate Python scripts entitled vol1_oofem.py and vol1_yade.py.
# This user interface with this script is mainly the definition of the number of simulation steps, critical time step (equal to the 
#one defined in vol1_yade) and the interval of output request (nSteps,dt,output, respectively). Please do not change other sections.
######################################################################
#####################################################
#
# volume fem-dem coupling test
#
#####################################################

# import interfaces
import liboofem
import libyade
from demfemcoupling import OofemInterface,YadeInterface,OofemYadeMeshSurfaceMap,FemDemSurfaceCoupler

def vtkExport(i,fem,dem):
	"""Do VTK export"""
	fem.vtkExport(i)
	dem.vtkExport(i)

# initialize both domains
femName = 'vol1_oofem'
demName = 'vol1_yade'
fem = OofemInterface(femName,liboofem)
dem = YadeInterface(demName,libyade)

# create coupler object
femSurf = fem.toUnstructuredGrid().getSurface()
demSurf = dem.addSurface(femSurf)
femDemMeshMap = OofemYadeMeshSurfaceMap(fem,dem,femSurf,demSurf)
coupler = FemDemSurfaceCoupler(fem,dem,femDemMeshMap,vtkExport)

# run the simulation
nSteps,dt,output = 2000000,1e-7,50000
coupler.solve(nSteps,dt,output)
