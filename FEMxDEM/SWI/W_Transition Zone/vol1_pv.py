#Copyright 2023 Blueshift, LLC
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, #including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to #do so, subject to the following conditions:
         #The Software is subject to all use, distribution, modification, sales, and other restrictions applicable to the software-as-a-service product specified in the Agreement.
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND #NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR #IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################
# This Python script generates a synchronized FEMxDEM VTK output to run via paraview. No adjustments required
# here. Please see the instruction for information about running this script via paraview. 
######################################################################
from __future__ import print_function
import os
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

directory = "/tmp"

renderView = GetActiveViewOrCreate('RenderView')

# FEM
fem = [os.path.join(directory,f) for f in os.listdir(directory) if f.startswith("vol1_oofem.out.m1.") and f.endswith(".vtu")]
fem.sort(key=lambda f: int(os.path.splitext(f)[0].rpartition(".")[2]))
fem = XMLUnstructuredGridReader(FileName=fem)


Show(fem, renderView)
warpByVector1 = WarpByVector(Input=fem)
warpByVector1Display = Show(warpByVector1, renderView)
Hide(fem, renderView)
warpByVector1Display.SetRepresentationType('Surface With Edges')
ColorBy(warpByVector1Display, ('POINTS', 'IST_StressTensor', '8'))

iSTStressTensorLUT = GetColorTransferFunction('IST_StressTensor')
iSTStressTensorLUT.ApplyPreset('Warm to Cool', True)
iSTStressTensorLUT.RescaleTransferFunction(-300000.0, 300000.0)

iSTStressTensorPWF = GetOpacityTransferFunction('IST_StressTensor')
iSTStressTensorPWF.RescaleTransferFunction(-300000.0, 300000.0)

# DEM/facets
facets = [os.path.join(directory,f) for f in os.listdir(directory) if f.startswith("vol1_yade-facets-") and f.endswith(".vtk")]
facets.sort(key=lambda f: int(os.path.splitext(f)[0].rpartition("-")[2]))
facets = LegacyVTKReader(FileNames=facets)
Show(facets, renderView)
warpByVector2 = WarpByVector(Input=facets)
warpByVector2Display = Show(warpByVector2, renderView)
Hide(facets, renderView)
warpByVector2Display.SetRepresentationType('Surface With Edges')

# DEM/spheres
dem = [os.path.join(directory,f) for f in os.listdir(directory) if f.startswith("vol1_yade-spheres-") and f.endswith(".vtk")]
dem.sort(key=lambda f: int(os.path.splitext(f)[0].rpartition("-")[2]))
dem = LegacyVTKReader(FileNames=dem)

Show(dem, renderView)
glyph1 = Glyph(Input=dem, GlyphType='Sphere')
glyph1.ScaleMode = 'scalar'
glyph1.ScaleFactor = 1.0
glyph1.GlyphType.Radius = 1.0
glyph1.GlyphType.ThetaResolution = 16
glyph1.GlyphType.PhiResolution = 16
glyph1.GlyphMode = 'All Points'

glyph1Display = Show(glyph1, renderView)
ColorBy(glyph1Display, ('POINTS', 'dmg'))

dmgLUT = GetColorTransferFunction('dmg')
dmgLUT.ApplyPreset('Blue to Red Rainbow', True)
dmgLUT.RescaleTransferFunction(0,.5)

dmgPWF = GetOpacityTransferFunction('dmg')
dmgPWF.RescaleTransferFunction(0,.5)

#
animationScene1 = GetAnimationScene()
animationScene1.UpdateAnimationUsingDataTimeSteps()

renderView.ResetCamera()
renderView.OrientationAxesVisibility = False
renderView.Background = [1.0, 1.0, 1.0]
renderView.ViewSize = [989, 322]
renderView.CameraPosition = [0.0, 0.0665725152939558, -1.9729839408780194]
renderView.CameraFocalPoint = [0.0, 0.0665725152939558, 0.0]
renderView.CameraViewUp = [0.0, 1.0, 0.0]
renderView.CameraParallelProjection = True
renderView.CameraParallelScale = 0.169536833167

RenderAllViews()
out = '/tmp/vol1.png'
WriteAnimation(out)
print('animation saved to {}'.format(out))
print()
