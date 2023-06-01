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
"""
This script generates the grain as clumps and stand alone spheres.
Volume of each grain is calculated from the sum of its constitute sphere. This is used for porosity calculations

Grains size and shape distributions are obtained from the csv or txt files in  "InputData"

LMA1 grains are made of elipsoids that are filled by spheres; The spheres are then randomly removed; 
The removal percentage is based on the ratio of the grain cross section area (A_g) and the area of the ellipse circumbscribing the grain (A_e). 
; that is, merely the mean and standard deviation of (1- A_g/A_e) for each bin is used to find the percentage of the total spheres that should be removed from a grain 

GRC3 grains are made of polyhedra. c/a, b/a, and c/b in these polyhedra are specified based on their distribution that is obtained from Input Data File but each polyhedron is randomly formed. 

The number of spheres to generate each grain is defined in a way that a grain is made up of two or more spheres only if two spheres of the minimum size grain can fit into it.
For example for the case that minimum size grain is 90 microns, each grain is made up of one single sphere unless more than one sphere of size 90 microns can be fit into it. The GCSC parameter in the inputData file defines how many spheres can be fit into the intermeduate axis of grains. Increasing this value increases the numbe rof spheres used to generate grains, so increases
the grain shape complexity while increasing runtime. This also means, the spheres do not have the same size in different grains and thus they scale with the grain size. 

Output of this script is a txt file containing the spheres with the grain type, grain id, and volume of the grains that they belong to. 

"""

######################################################################
#Importing necessary, built-in functions to Yade
######################################################################
from __future__ import print_function
from yade import export,ymport
import random
from yade import polyhedra_utils
from yade import plot
from yade import pack
from yade import qt,timing
from builtins import range
from builtins import object
from yade.wrapper import *
from yade import utils,Matrix3,Vector3
from decimal import Decimal
from scipy.spatial import ConvexHull

######################################################################
#Importing user-defined variables from "InputData" script to this script like Particle Size Distribution
######################################################################
from InputData import * 
# This scripts contains the input data for soil grain generator and bond installer scripts


#####################################################################
Mat1=FrictMat(frictionAngle=FRICTShallow)
#####################################################################

#_______________________________________________________________________________________________________
def textExtAPriodic(filename, format='x_y_z_r', comment='',mask=-1,attrs=[],vol=[],GrainType=[],GrainSplit=[],GrainID=[]):
	O=Omega()
	
	try:
		out=open(filename,'w')
	except:
		raise RuntimeError("Problem to write into the file")
	
	count=0
	
	# TODO use output=[] instrad of ''???
	output = ''
	outputVel=''
	if (format!='liggghts_in'):
		output = '#format ' + format + '\n'
		if (comment):
			if format=='x_y_z_r_attrs':
				cmts = comment.split('\n')
				for cmt in cmts[:-1]:
					output += cmt
				output += '# x y z r ' + cmts[-1] + '\n'
			else:
				output += '# ' + comment + '\n'
	
	minCoord= Vector3.Zero
	maxCoord= Vector3.Zero
	maskNumber = []
	
	for b in O.bodies:
		try:
			if (isinstance(b.shape,Sphere) and ((mask<0) or ((mask&b.mask)>0))):
				if (format=='x_y_z_r'):
					output+=('%g\t%g\t%g\t%g\n'%(b.state.pos[0],b.state.pos[1],b.state.pos[2],b.shape.radius))
				elif (format=='x_y_z_r_matId'):
					output+=('%g\t%g\t%g\t%g\t%d\n'%(b.state.pos[0],b.state.pos[1],b.state.pos[2],b.shape.radius,b.material.id))
				elif (format=='x_y_z_r_attrs'):
					if b.isClumpMember==0:
						output+=('%g\t%g\t%g\t%g'%(b.state.pos[0],b.state.pos[1],b.state.pos[2],b.shape.radius))
						for cmd in attrs:
							v = eval(cmd)
							if isinstance(v,(int,float)):
								output+='\t%g'%v
							elif isinstance(v,Vector3):
								output+='\t%g\t%g\t%g'%tuple(v[i] for i in range(3))
							elif isinstance(v,Matrix3):
								output+='\t%g'%tuple(v[i] for i in range(9))
						output += '\n'
				elif (format=='x_y_z_r_attrs_vol_GrainType_GrainSplit_GrainID'):
					if b:#b.isClumpMember==0 and 
						output+=('%g\t%g\t%g\t%g'%(b.state.pos[0],b.state.pos[1],b.state.pos[2],b.shape.radius))
						for cmd in attrs:
							v = eval(cmd)
							if isinstance(v,(int,float)):
								output+='\t%g'%v
							elif isinstance(v,(str)):
								output+='\t%s'%v
							elif isinstance(v,Vector3):
								output+='\t%g\t%g\t%g'%tuple(v[i] for i in range(3))
							elif isinstance(v,Matrix3):
								output+='\t%g'%tuple(v[i] for i in range(9))
						#output += '\n'
						for cmd in vol:
							v = eval(cmd)
							if isinstance(v,(int,float)):
								output+='\t%g'%v
							elif isinstance(v,(str)):
								output+='\t%s'%v
							elif isinstance(v,Vector3):
								output+='\t%g\t%g\t%g'%tuple(v[i] for i in range(3))
							elif isinstance(v,Matrix3):
								output+='\t%g'%tuple(v[i] for i in range(9))
						for cmd in GrainType:
							v = eval(cmd)
							if isinstance(v,(int,float)):
								output+='\t%g'%v
							elif isinstance(v,(str)):
								output+='\t%s'%v
							elif isinstance(v,Vector3):
								output+='\t%g\t%g\t%g'%tuple(v[i] for i in range(3))
							elif isinstance(v,Matrix3):
								output+='\t%g'%tuple(v[i] for i in range(9))
						for cmd in GrainSplit:
							v = eval(cmd)
							if isinstance(v,(int,float)):
								output+='\t%g'%v
							elif isinstance(v,(str)):
								output+='\t%s'%v
							elif isinstance(v,Vector3):
								output+='\t%g\t%g\t%g'%tuple(v[i] for i in range(3))
							elif isinstance(v,Matrix3):
								output+='\t%g'%tuple(v[i] for i in range(9))
						for cmd in GrainID:
							v = eval(cmd)
							if isinstance(v,(int,float)):
								output+='\t%g'%v
							elif isinstance(v,(str)):
								output+='\t%s'%v
							elif isinstance(v,Vector3):
								output+='\t%g\t%g\t%g'%tuple(v[i] for i in range(3))
							elif isinstance(v,Matrix3):
								output+='\t%g'%tuple(v[i] for i in range(9))
						output += '\n'
				elif (format=='x_y_z_r_attrs'):
					if b.isClumpMember==0 and b:
						output+=('%g\t%g\t%g\t%g'%(O.cell.wrap(b.state.pos)[0],O.cell.wrap(b.state.pos)[1],O.cell.wrap(b.state.pos)[2],b.shape.radius))
						for cmd in attrs:
							v = eval(cmd)
							if isinstance(v,(int,float)):
								output+='\t%g'%v
							elif isinstance(v,(str)):
								output+='\t%s'%v
							elif isinstance(v,Vector3):
								output+='\t%g\t%g\t%g'%tuple(v[i] for i in range(3))
							elif isinstance(v,Matrix3):
								output+='\t%g'%tuple(v[i] for i in range(9))
						output += '\n'
				elif (format=='id_x_y_z_r_matId'):
					output+=('%d\t%g\t%g\t%g\t%g\t%d\n'%(b.id,b.state.pos[0],b.state.pos[1],b.state.pos[2],b.shape.radius,b.material.id))
				elif (format=='jointedPM'):
					output+=('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n'%(b.id,b.state.onJoint,b.state.joint,b.state.jointNormal1[0],b.state.jointNormal1[1],b.state.jointNormal1[2],b.state.jointNormal2[0],b.state.jointNormal2[1],b.state.jointNormal2[2],b.state.jointNormal3[0],b.state.jointNormal3[1],b.state.jointNormal3[2]))
				elif (format=='liggghts_in'):
					output+=('%g %g %g %g %g %g %g\n'%(count+1,b.mask,b.shape.radius,b.material.density,b.state.pos[0],b.state.pos[1],b.state.pos[2]))
					outputVel+=('%g %g %g %g %g %g %g\n'%(count+1,b.state.vel[0],b.state.vel[1],b.state.vel[2],b.state.angVel[0],b.state.angVel[1],b.state.angVel[2]))
				else:
					raise RuntimeError("Please, specify a correct format output!");
				count+=1
				if  (count==1):
					minCoord = b.state.pos - Vector3(b.shape.radius,b.shape.radius,b.shape.radius)
					maxCoord = b.state.pos + Vector3(b.shape.radius,b.shape.radius,b.shape.radius)
				else:
					minCoord = Vector3(min(minCoord[0], b.state.pos[0]-b.shape.radius),min(minCoord[1], b.state.pos[1]-b.shape.radius),min(minCoord[2], b.state.pos[2]-b.shape.radius))
					maxCoord = Vector3(max(maxCoord[0], b.state.pos[0]+b.shape.radius),max(maxCoord[1], b.state.pos[1]+b.shape.radius),max(minCoord[2], b.state.pos[2]+b.shape.radius))
				if b.mask not in maskNumber:
					maskNumber.append(b.mask)
		except AttributeError:
			pass
			
	out.write(output)
	out.close()
	return count

#_______________________________________________________________________________________________________

######################################################################
#Generating Polyhedra with specified shape index but random shape. 
#This finction receives a stl file that is later generated by other functions
#a, b, and c are the lenght in major, intermediate, and minor axes
#oz is the coordinate at which the polyhedra are generated; it is not critical to simulation and can be arbitrary  
#PolMat is defined here because it has no effect on the final grains
#This function returns the initial volume of the generated polyhedron. This volume is later scaled in "Filler" based on the size distribution
######################################################################
PolMat = PolyhedraMat()
PolMat.IsSplitable = True
PolMat.strength = 7.9e6
PolMat.density = 2678#kg/m^3
PolMat.young = 5.98e7 #Pa
PolMat.poisson =0.3
PolMat.frictionAngle = 0.5 #rad

def polyhedra2stl(stl,a,b,c,oz):
#create a polyhedron
#a, b, c are major, intermediate, and minor axes lengths (diameters)
  poly=polyhedra_utils.polyhedra(material=PolMat,size=(a,b,c))
  poly.shape.wire = True
  O.bodies.append(poly)
  PolVolume=poly.shape.GetVolume()
  poly.state.pos=(0,0,oz)
  trisCoords = []
  for b in O.bodies:
    if not isinstance(b.shape,Polyhedra):
      continue
    vs = [b.state.pos + b.state.ori*v for v in b.shape.v]
    trisIds = b.shape.GetSurfaceTriangulation()
    l = int(len(trisIds)/3)
    
    trisIds = [trisIds[3*i:3*i+3] for i in range(l)]
    for t in trisIds:
      trisCoords.append([vs[i] for i in t])
  #Turning polyhedra to stl files
  lines = ['solid polyhedra']
  for v1,v2,v3 in trisCoords:
    n = (v2-v1).cross(v3-v1)
    lines.append('facet normal {} {} {}'.format(n[0],n[1],n[2]))
    lines.append(' outer loop')
    for v in (v1,v2,v3):
      lines.append(' vertex {} {} {}'.format(v[0],v[1],v[2]))
    lines.append(' endloop')
    lines.append('endfacet')
  lines.append('endsolid polyhedra')
  #
  with open(stl,'w') as f:
    f.writelines(l+'\n' for l in lines)
  return PolVolume

import numpy as np
import math
# Function to find distance
def shortest_distance(x1, y1, z1, a, b, c, d):
     
    d = abs((a * x1 + b * y1 + c * z1 + d))
    e = (math.sqrt(a * a + b * b + c * c))
    #print("Perpendicular distance is", d/e)
    return(d/e)
#__________________________________________________________________________________________________
def ExternalGrainVolume(SpheresData):
 CenterPoints = []
 ClusterDataforOuterConvex = []
 SolidPartofInnerConvex = []

##Data need for the outer convexhull
 SpheresVolume=0
 for p2 in range (0,len(SpheresData)):
  SpherePoints = []
  X = SpheresData[p2,0]
  Y = SpheresData[p2,1]
  Z = SpheresData[p2,2]
  R = SpheresData[p2,3]
  SpheresVolume += 4.*math.pi/3.*R**3.
  Number=20

##Random placement of points on surface of sphere
  for nu in numpy.arange(0, (math.pi), (math.pi/(2*Number))):
    for phi in numpy.arange(0, (2*math.pi), (math.pi/Number)):
      x = R*numpy.sin(nu)*numpy.cos(phi)
      y = R*numpy.sin(nu)*numpy.sin(phi)
      z = R*numpy.cos(nu)
      ClusterDataforOuterConvex.append(numpy.array([float(X+x),float(Y+y),float(Z+z)]))
 ClusterDataforOuterConvex = numpy.array(ClusterDataforOuterConvex).reshape(-1,3)
 OuterConvexVolume = ConvexHull(ClusterDataforOuterConvex).volume

 return (OuterConvexVolume) # OuterConvexVolume
 
 
 
 
 ####################################################################
def InternalGrainVolume(SpheresData):
 #print("SpheresData",SpheresData)
 CenterPoints = []
 ClusterDataforOuterConvex = []
 SolidPartofInnerConvex = []

 #SpheresData=SpheresData*(1e9)

##Generate the inner convexhull from particles centerpoint
 for p1 in range (0,len(SpheresData)):
  CenterPoints.append(numpy.array([SpheresData[p1,0],SpheresData[p1,1],SpheresData[p1,2]]))
 if len(SpheresData)>3:
  CenterPoints = numpy.array(CenterPoints).reshape(-1,3)
  #print("this is the number of points",CenterPoints)
  _p1 = np.array([SpheresData[0][0],SpheresData[0][1],SpheresData[0][2]])
  _p2 = np.array([SpheresData[1][0],SpheresData[1][1],SpheresData[1][2]])
  _p3 = np.array([SpheresData[2][0],SpheresData[2][1],SpheresData[2][2]])
# These two vectors are in the plane
  _v1 = _p3 - _p1
  _v2 = _p2 - _p1

# the cross product is a vector normal to the plane
  _cp = np.cross(_v1, _v2)
  _a, _b, _c = _cp

# This evaluates a * x3 + b * y3 + c * z3 which equals d
  _d = np.dot(_cp, _p3)

  #print('The equation is {0}x + {1}y + {2}z = {3}'.format(_a, _b, _c, _d))
  DistList=[]; _Checker=0
  for io in range(len(SpheresData)):
   if io>2:
    _Distance=shortest_distance(SpheresData[io][0], SpheresData[io][1], SpheresData[io][2], _a, _b, _c, -_d)
    #print("This is the distance:",_Distance)
    DistList.append(_Distance)
  print("DistList",DistList,"SpheresData[0][3]",SpheresData[0][3])
  for idi in  DistList:
   if abs(idi)<2.*SpheresData[0][3]:
    _Checker+=1
  if len(DistList)-_Checker==0:
   TotalVolumeofGrain=ExternalGrainVolume(SpheresData)
   print("One grain is planar *************************************")
  else:
   InnerConvexVertices = ConvexHull(CenterPoints).vertices
   #print("this is the InnerConvexVertices",InnerConvexVertices)

   InnerConvexPointsofFacets = ConvexHull(CenterPoints).simplices
   InnerConvexFacets = ConvexHull(CenterPoints).equations
   InnerConvexVolume = ConvexHull(CenterPoints).volume

##Data need for the outer convexhull
   SpheresVolume=0
   for p2 in range (0,len(SpheresData)):
    SpherePoints = []
    X = SpheresData[p2,0]
    Y = SpheresData[p2,1]
    Z = SpheresData[p2,2]
    R = SpheresData[p2,3]
    SpheresVolume += 4.*math.pi/3.*R**3.
    Number=20

##Random placement of points on surface of sphere
    for nu in numpy.arange(0, (math.pi), (math.pi/(2*Number))):
      for phi in numpy.arange(0, (2*math.pi), (math.pi/Number)):
        x = R*numpy.sin(nu)*numpy.cos(phi)
        y = R*numpy.sin(nu)*numpy.sin(phi)
        z = R*numpy.cos(nu)
        ClusterDataforOuterConvex.append(numpy.array([float(X+x),float(Y+y),float(Z+z)]))
        SpherePoints.append(numpy.array([float(X+x),float(Y+y),float(Z+z)]))
    SpherePoints = numpy.array(SpherePoints).reshape(-1,3)

##Check if a point on the sphere surface is inside the inner convexhull
    InsidePoints =[]
    IntersectPoint=[]
    for p3 in range(0,len(SpherePoints)):
      SinglePoint = numpy.array([SpherePoints[p3,0],SpherePoints[p3,1],SpherePoints[p3,2]]).reshape(-1,3)
      NewInnerConvexVertices = ConvexHull(numpy.concatenate((CenterPoints, SinglePoint))).vertices
    #print(NewInnerConvexVertices)

      if numpy.array_equal(NewInnerConvexVertices, InnerConvexVertices):
        InsidePoints.append(SinglePoint)
    InsidePoints = numpy.array(InsidePoints).reshape(-1,3)

##Find the intersection points between sphere p2 and others 
    for p4 in range (0,len(SpheresData)):
      Distance = math.sqrt((SpheresData[p2,0]-SpheresData[p4,0])**2+(SpheresData[p2,1]-SpheresData[p4,1])**2+(SpheresData[p2,2]-SpheresData[p4,2])**2)
      Eps = R/100

      if abs(Distance-(SpheresData[p2,3]+SpheresData[p4,3]))<=Eps:
        IntersectPoint.append(numpy.array([(SpheresData[p2,0]+SpheresData[p4,0])/2,(SpheresData[p2,1]+SpheresData[p4,1])/2,(SpheresData[p2,2]+SpheresData[p4,2])/2]))
    IntersectPoint = numpy.array(IntersectPoint).reshape(-1,3)
    Center = numpy.array([SpheresData[p2,0],SpheresData[p2,1],SpheresData[p2,2]]).reshape(-1,3)
    EdgeandCenterPoints=numpy.concatenate((IntersectPoint, Center)).reshape(-1,3)
##Volume of solid parts inside inner convexhull
    SolidPartofInnerConvex.append(ConvexHull(numpy.concatenate((InsidePoints, EdgeandCenterPoints))).volume)

   SolidPartofInnerConvex = numpy.array(SolidPartofInnerConvex).reshape(-1)

##Points to generate the outer convexhull

   InnerVoidVolume = InnerConvexVolume-sum(SolidPartofInnerConvex)
   TotalVolumeofGrain = (InnerVoidVolume+SpheresVolume)
   TotalInnerVoidVolume = InnerVoidVolume
   print("Grain with Internal Volume Cacs is generated===========================")
 else:
  TotalVolumeofGrain=ExternalGrainVolume(SpheresData)
  print("Grain with less than 4 spheres +++++++++++++++++++++++++++++++++++++++++")
 return (TotalVolumeofGrain)


#In this section, the functions are called in a specific order
#_______________________________________________________________________________________________________
######################################################################
#This function receives the randomly generated polyhedron and fills it with spheres
#Sphere radius is defined by "TempRad" a user-defined value based on the cutoff size, which is a passed variable
#"stl" is the name of the stl file generated previously by "polyhedra2stl"
#PolVolume comes from polyhedra2stl
#The return is the polyhedron volume scaled with respect to the sum of sphere volumes within the polyhedron
#later, we scale the volume back based on the size distribution of each grain
# Another return is a tempelate for generating grains with the polyhedron shape
# This template will be used later to generate grains following this shape and the grain size distribution
######################################################################
def Filler(stl,TempRad,plc,plb,PolVolume):
 BAxis=max(plc,plb)
 global templatepl
 facets = ymport.stl(stl)
 surf = gts.Surface()
 for f in facets:
  vs = [f.state.pos + f.state.ori*v for v in f.shape.vertices]
  gtsvs = [gts.Vertex(v[0],v[1],v[2]) for v in vs]
  es = [gts.Edge(gtsvs[i],gtsvs[j]) for i,j in ((0,1),(1,2),(2,0))]
  face = gts.Face(es[0],es[1],es[2])
  surf.add(face)
 surf.cleanup(1.e-6)
 assert surf.is_closed()
 pred = pack.inGtsSurface(surf)
 
 for ib in range(99):
  Radif=TempRad*(100.-ib)/100.
  sphs=pack.regularHexa(pred,radius=Radif,gap=0.)
  if len(sphs)>0:
   break
   
 if len(sphs)<2:
  if plc/2.>1.e-4*1.e5:
    for ib in range(99):
     sphs=pack.regularHexa(pred,radius=(Radif*(100.-ib)/100.),gap=0.)
     if len(sphs)>1:
      break
   
 for ib2 in range(15000):
  if Radif*(100.+ib2)/100.<BAxis/GCSC/2. and len(sphs)>4:
     sphs=pack.regularHexa(pred,radius=(Radif*(100.+ib2)/100.),gap=0.)
     #print("len(sphs)",len(sphs),"BAxis/GCSC/2",BAxis/GCSC/2.,"Rad",(Radif*(100.+ib2)/100.),"TempRad",TempRad)
  if BAxis/GCSC/2.<Radif*(100.+ib2)/100. or len(sphs)<5:#:#(Rad*(100.+ib2)/100.):# Note plc is diameter
   break
 #print("len(sphs)",len(sphs))

 
 if len(sphs)>1:
  #ids=O.bodies.append(sphs)
  nn=0;ij=0
  plVolSp=0.
  relRadList2=[]
  relPosList2=[]
  j=0
  for i in sphs:
   relRadList2.append(i.shape.radius)
   relPosList2.append(i.state.pos)
   j+=1 
  templatepl= []
  templatepl.append(clumpTemplate(relRadii=relRadList2,relPositions=relPosList2))
#  for i in sphs:
#   O.bodies.erase(i.id)
  for b in O.bodies:
   if (isinstance(b.shape,Polyhedra)): 
    O.bodies.erase(b.id)
#    ij+=1
  plVolume=PolVolume

  return (1,plVolume)
 else:
  for b in O.bodies:
   if (isinstance(b.shape,Polyhedra)): 
    O.bodies.erase(b.id)
  return (2,69)#AgVolume

######################################################################
#This function receives size (a, b, c) and generate ellipsoids at location oz (dummy)
#The returns are the volume of the generated ellipsoid and a template for grain generation for Agglutinates
######################################################################

#_______________________________________________________________________________________________________
 
def ellipsoidtoAg(a,b,c,oz):
#a, b, c should be in radii
 #a= a/2.
 #b= b/2.
 #c= c/2.
 BAxis=max(c,b)
 pred=pack.inEllipsoid((0.,0.,oz),(a,max(c,b),min(c,b)))
 global templateAg
 TempRad=0.00009/2.
 #print("a,b,c,TempRad:",a,b,c,TempRad)
 for ib in range(99):
  Rad=TempRad*(100.-ib)/100.
  sphAg=pack.regularHexa(pred,radius=(Rad),gap=0.)
  if len(sphAg)>0:
   break
 if len(sphAg)<2:
  if c>1.e-4:
    for ib in range(99):
     Rad=Rad*(100.-ib)/100.
     sphAg=pack.regularHexa(pred,radius=(Rad),gap=0.)
     if len(sphAg)>1:
      break
 for ib in range(15000):
  if Rad*(100.+ib)/100.<BAxis/GCSC and len(sphAg)>7:
   sphAg=pack.regularHexa(pred,radius=(Rad*(100.+ib)/100.),gap=0.)
   #print("len(sphAg)",len(sphAg),"(Rad*(100.+ib)/100.)",(Rad*(100.+ib)/100.),"BAxis/GCSC",BAxis/GCSC,a,b,c)
  if Rad*(100.+ib)/100.>BAxis/GCSC or len(sphAg)<5:#len(sphAg)<50:#
   break



 #print("len(sphAg)",len(sphAg))
 if len(sphAg)>1:
  relRadList1=[]
  relPosList1=[]
  AgVolSp=0
  for i in sphAg:
   relRadList1.append(i.shape.radius)
   relPosList1.append(i.state.pos)
   AgVolSp+=4.*math.pi/3.*(i.shape.radius)**3.
  templateAg= []
  AgTem=templateAg.append(clumpTemplate(relRadii=relRadList1,relPositions=relPosList1))
  AgVolume=a*b*c*4.*math.pi/3.
  AgVolume=AgVolume/AgVolSp
 #print (AgVolume,AgVolSp)
  return (1,templateAg)#AgVolume
 if len(sphAg)<=1:
  return (2,sphAg)#AgVolume
 
#________________________________________________________________________________________________
######################################################################
#This function receives the generated ellipsoid and scale it based on the size distribution
#Some spheres are then removed based on the ratio of the solid ag to the volume of the elipsoid (Imported from the "InputData" script)
#Returns are an array for the generated spheres and the volume of the grain
######################################################################
import statistics 

#################################################

LMA1psd=[]
LMA1data=[]

GRC3psd=[]
GRC3data=[]
#################################################
infile = open(LMA1fileName,"r")
LMA1lines = infile.readlines()
LMA1Grains = [[] for _ in range(len(LMA1lines)-2)]
infile.close()

infile = open(GRC3fileName,"r")
GRC3lines = infile.readlines()
GRC3Grains = [[] for _ in range(len(GRC3lines)-2)]
infile.close()

#################################################
iaL=0
for line in LMA1lines:#line in lines:
	LMA1data = line.split()
	#print("LMA1data",LMA1data)
	if (LMA1data[0] == "Id" or LMA1data[2] == "µm") :
		continue
	LMA1Grains[iaL].append(LMA1data)
	iaL+=1
print("The number of LMA1 grains imported is:",iaL)

#################################################
ia=0
for line in GRC3lines:#line in lines:
	GRC3data = line.split()
	#print("GRC3data",GRC3data)
	if GRC3data==[]: continue
	if (GRC3data[0] == "Id"  or GRC3data[2] == "µm" ):
		continue
	GRC3Grains[ia].append(GRC3data)
	ia+=1
print("The number of GRC3 grains imported is:",ia)
TIGN=ia+iaL
print("Total number of imported grains  is:",ia+iaL)
#################################################
for gr in range(len(LMA1Grains)):
	for igr in range(len(LMA1Grains[0][0])):
		if not(LMA1Grains[gr][0][igr]=='Reject' or LMA1Grains[gr][0][igr]=='Accept') :
			LMA1Grains[gr][0][igr]=float(LMA1Grains[gr][0][igr])

for gr in range(len(GRC3Grains)):
	if GRC3Grains[gr]==[]:continue
	for igr in range(len(GRC3Grains[0][0])):
		if not(GRC3Grains[gr][0][igr]=='Reject' or GRC3Grains[gr][0][igr]=='Accept') :
			GRC3Grains[gr][0][igr]=float(GRC3Grains[gr][0][igr])
#################################################
#Filling PSD
for i in range(LMA1VFs):
	LMA1psd.append(LMA1MinSize+(LMA1MaxSize-LMA1MinSize)/(LMA1VFs-1)*i)

for i in range(GRC3VFs):
	GRC3psd.append(GRC3MinSize+(GRC3MaxSize-GRC3MinSize)/(GRC3VFs-1)*i)
#################################################
LMA1LenRatio = [[] for _ in range(LMA1VFs-1)]
LMA1LenRatioMean = [[] for _ in range(LMA1VFs-1)]
LMA1LenRatioSD = [[] for _ in range(LMA1VFs-1)]
LMA1EllipsoidalVoid = [[] for _ in range(LMA1VFs-1)]
LMA1EllipsoidalVoidMean = [[] for _ in range(LMA1VFs-1)]
LMA1EllipsoidalVoidSD = [[] for _ in range(LMA1VFs-1)]
LMA1VolFracAg = [[] for _ in range(LMA1VFs-1)]

GRC3LenRatio = [[] for _ in range(GRC3VFs-1)]
GRC3LenRatioMean = [[] for _ in range(GRC3VFs-1)]
GRC3LenRatioSD = [[] for _ in range(GRC3VFs-1)]
GRC3VolFrac = [[] for _ in range(GRC3VFs-1)]
#################################################
#Calculate Total Volume
LMA1TGV=0.
for grainLine in LMA1Grains:
	if grainLine[0][2]>LMA1MinSize and grainLine[0][2]<LMA1MaxSize:
		LMA1TGV+=grainLine[0][8]/1.e18#[0][8] is the grain Volume

GRC3TGV=0.
for grainLine in GRC3Grains:
	if grainLine==[]:continue
	if grainLine[0][2]>GRC3MinSize and grainLine[0][2]<GRC3MaxSize:
		GRC3TGV+=grainLine[0][8]/1.e18#[0][8] is the grain Volume

#################################################
def Gama():
 global Requested_LMA1_Volume; global Requested_GRC3_Volume
 if LMA1_Volume_To_Total_Volume_Percentage==1:
  Requested_LMA1_Volume=RequestedSampleVolume# for generation of exact input grains use: LMA1TGV#
  Requested_GRC3_Volume=0.# for generation of exact input grains use: GRC3TGV#
 else:
  Gama=LMA1_Volume_To_Total_Volume_Percentage/(1.-LMA1_Volume_To_Total_Volume_Percentage)*Density_GRC3/Density_LMA1
  Requested_LMA1_Volume=RequestedSampleVolume*Gama/(1.+Gama)# for generation of exact input grains use: LMA1TGV#
  Requested_GRC3_Volume=RequestedSampleVolume/(1.+Gama)# for generation of exact input grains use: GRC3TGV#

LMA1_Volume_To_Total_Volume_Percentage=LMA1_Mass_To_Total_Mass_Percentage
Gama()

#################################################
f= open("InputLMA1SphereDiameter.txt", 'w')
hi=0
for grainLine in LMA1Grains:
	hi+=1
	f.write(' %f' % grainLine[0][2])
	f.write('\n')
f.close()

f= open("InputGRC3SphereDiameter.txt", 'w')
hi=0
for grainLine in GRC3Grains:
	if grainLine==[]:continue
	hi+=1
	#print(grainLine[0][2])
	f.write(' %f' % grainLine[0][2])
	f.write('\n')
f.close()


#################################################

#Calculate Aspet ratio and volume fraction
for psdBin in range(LMA1VFs-1):
	LMA1GV=0.
	for grainLine in LMA1Grains:
		if grainLine[0][2]>=LMA1psd[psdBin] and grainLine[0][2]<LMA1psd[psdBin+1]:#[0][2] is the grain Da
			LMA1LenRatio[psdBin].append(grainLine[0][7]/grainLine[0][6])#[0][7] is the grain Ewidth
			LMA1EllipsoidalVoid[psdBin].append(grainLine[0][9]/grainLine[0][6]/grainLine[0][7]/math.pi*4.)#[0][7] is the grain Ewidth
			LMA1GV+=grainLine[0][8]/1.e18#[0][8] is the grain Volume
	LMA1VolFracAg[psdBin].append(LMA1GV/LMA1TGV)
for psdBin in range(LMA1VFs-1):
	if len(LMA1LenRatio[psdBin])>1:
		LMA1LenRatioMean[psdBin].append(statistics.mean(LMA1LenRatio[psdBin]))
		LMA1LenRatioSD[psdBin].append(statistics.stdev(LMA1LenRatio[psdBin]))
		LMA1EllipsoidalVoidMean[psdBin].append(statistics.mean(LMA1EllipsoidalVoid[psdBin]))
		LMA1EllipsoidalVoidSD[psdBin].append(statistics.stdev(LMA1EllipsoidalVoid[psdBin]))
	if len(LMA1LenRatio[psdBin])==1:
		LMA1LenRatioMean[psdBin].append(LMA1LenRatio[psdBin][0])
		LMA1LenRatioSD[psdBin].append(0)
		LMA1EllipsoidalVoidMean[psdBin].append(LMA1EllipsoidalVoid[psdBin][0])
		LMA1EllipsoidalVoidSD[psdBin].append(0)
ft=0.
for ih in LMA1VolFracAg:
	ft+=ih[0]
print("Sum of all Volume Fractions for LMA1 is",ft, "there is a red flag if this number is less than 0.99999 ")




for psdBin in range(GRC3VFs-1):
	GRC3GV=0.
	for grainLine in GRC3Grains:
		if grainLine==[]:continue
		if grainLine[0][2]>=GRC3psd[psdBin] and grainLine[0][2]<GRC3psd[psdBin+1]:#[0][2] is the grain Da
			GRC3LenRatio[psdBin].append(grainLine[0][7]/grainLine[0][6])#[0][7] is the grain Ewidth
			GRC3GV+=grainLine[0][8]/1.e18#[0][8] is the grain Volume
	GRC3VolFrac[psdBin].append(GRC3GV/GRC3TGV)
for psdBin in range(GRC3VFs-1):
	if len(GRC3LenRatio[psdBin])>1:
		GRC3LenRatioMean[psdBin].append(statistics.mean(GRC3LenRatio[psdBin]))
		GRC3LenRatioSD[psdBin].append(statistics.stdev(GRC3LenRatio[psdBin]))
	if len(GRC3LenRatio[psdBin])==1:
		GRC3LenRatioMean[psdBin].append(GRC3LenRatio[psdBin][0])
		GRC3LenRatioSD[psdBin].append(0)
		
ft=0.
for ih in GRC3VolFrac:
	ft+=ih[0]
print("Sum of all Volume Fractions for GRC3 is",ft, "there is a red flag if this number is less than 0.99999 ")
minAspectRatioGRC3=1e6
for ii in range(len(GRC3LenRatioMean)):
	if GRC3LenRatioMean[ii] and GRC3LenRatioSD[ii]:
		if minAspectRatioGRC3>min(minAspectRatioGRC3,GRC3LenRatioMean[ii][0]-GRC3LenRatioSD[ii][0]):
			minAspectRatioGRC3=min(minAspectRatioGRC3,GRC3LenRatioMean[ii][0]-GRC3LenRatioSD[ii][0])
			hx=ii
			
minAspectRatioLMA1=1e6
for ii in range(len(LMA1LenRatioMean)):
	if LMA1LenRatioMean[ii] and LMA1LenRatioSD[ii]:
		if minAspectRatioLMA1>min(minAspectRatioLMA1,LMA1LenRatioMean[ii][0]-LMA1LenRatioSD[ii][0]):
			minAspectRatioLMA1=min(minAspectRatioLMA1,LMA1LenRatioMean[ii][0]-LMA1LenRatioSD[ii][0])
			hx=ii
#minAspectRatio=min([ii[0] for ii in GRC3LenRatioMean if ii])

print(minAspectRatioLMA1,hx,LMA1LenRatioMean[hx][0],LMA1LenRatioSD[hx][0])
MinAspectRatio=min(minAspectRatioLMA1,minAspectRatioGRC3)

print(MinAspectRatio)

#################################################
#Defining the initial size of the compaction box
LMA1TGN=RequestedSampleVolume/( 4./3.*math.pi*(min(GRC3MinSize,LMA1MinSize)/2./1.e6)**3.)
vv=(LMA1psd[len(LMA1psd)-1]/1.e6)**3  *LMA1TGN *100
NoOfBlocksOnEachSide=int((LMA1TGN)**(1/3)*1.2)
BlockLength=max(GRC3MaxSize,LMA1MaxSize)/2./1.e6*1.2*  2./(MinAspectRatio)**0.5 ####LMA1psd[len(LMA1psd)-1]/1.e6#(4./3.*math.pi*(psd[0]/2.)**3.)**1/3
print('NoOfBlocksOnEachSide',NoOfBlocksOnEachSide)
OccupiedBlock=[]
def RandomLocation():
	_Checker=0
	while _Checker==0:
		Xrand=random.randint(1,NoOfBlocksOnEachSide) 
		Yrand=random.randint(1,NoOfBlocksOnEachSide) 
		Zrand=random.randint(1,NoOfBlocksOnEachSide) 
		if not [Xrand,Yrand,Zrand] in OccupiedBlock: 
			OccupiedBlock.append([Xrand,Yrand,Zrand])
			_Checker=1
			CenterX=Xrand*BlockLength-0.5*BlockLength
			CenterY=Yrand*BlockLength-0.5*BlockLength
			CenterZ=Zrand*BlockLength-0.5*BlockLength
	return(Vector3(CenterX,CenterY,CenterZ))
#####################################################################
# Roughness estimation
#(numpy.random.normal(LenRatioMean[ii],LenRatioSD[ii]))
#####################################################################
#Ag

VolSoFar=0.
COF=0.;spMannual=[]

def AgGenerator():
	l=0;
	global ccsAgSmall; global crsAgSmall; global TLAgVol; global Agco; idsAg=[];  ret=[];keyyss=[]; ijh=0;global VolSoFar;VolSoFar=0.;global Aga;global Agb; global spMannualLMA1; spMannualLMA1=[]; global AgStandAlone; AgStandAlone=[]
	TLAgVol=0;Agco=0
	for ii in range(len(LMA1psd)-1):
		if LMA1VolFracAg[ii][0]>0 and LMA1LenRatioMean[ii] and LMA1LenRatioSD[ii]:
			TLAgVol=0.
			while TLAgVol<LMA1VolFracAg[ii][0]*Requested_LMA1_Volume:
				#print("LMA1VolFracAg[ii][0]*Requested_LMA1_Volume",LMA1VolFracAg[ii][0]*Requested_LMA1_Volume)
				ijh+=1
				keyyss=[];_Volume=0.;SpheresData=[]
				Rr=(random.uniform(LMA1psd[ii]/2., LMA1psd[ii+1]/2.))/1.e6
				SPp=utils.sphere(RandomLocation(),Rr)
				spMannual.append(SPp)
				spMannualLMA1.append(SPp)
				print("The ",ijh,"th Ag Grain with Diameter ",Rr*1e6*2., " microns is being generated from the ",ii,"th psd","with Volume:",TLAgVol)
				V_known=4./3.*math.pi*Rr**3
				Agc=Rr#
				c_a= numpy.random.normal(LMA1LenRatioMean[ii][0],LMA1LenRatioSD[ii][0])
				if c_a< (LMA1LenRatioMean[ii][0]-LMA1LenRatioSD[ii][0]):
					c_a=(LMA1LenRatioMean[ii][0]-LMA1LenRatioSD[ii][0])
				if c_a> (LMA1LenRatioMean[ii][0]+LMA1LenRatioSD[ii][0]):
					c_a=(LMA1LenRatioMean[ii][0]+LMA1LenRatioSD[ii][0])			
				c_b= numpy.random.normal(LMA1LenRatioMean[ii][0],LMA1LenRatioSD[ii][0])
				if c_b< (LMA1LenRatioMean[ii][0]-LMA1LenRatioSD[ii][0]):
					c_b=(LMA1LenRatioMean[ii][0]-LMA1LenRatioSD[ii][0])
				if c_b> (LMA1LenRatioMean[ii][0]+LMA1LenRatioSD[ii][0]):
					c_b=(LMA1LenRatioMean[ii][0]+LMA1LenRatioSD[ii][0])	
				Aga = Agc/c_a
				Agb = Aga*c_b
				#print("LenRatioMean[ii][0],LenRatioSD[ii][0]",LMA1LenRatioMean[ii][0],LMA1LenRatioSD[ii][0])
				Status=0
				if GCSC!=0:
					(Status,Returned)=ellipsoidtoAg(Aga,Agb,Agc,-0.1)
				if Status==1 and GCSC!=0:
					SphIDAg=O.bodies.append([SPp])
					clpMemIDAg=O.bodies.replaceByClumps(templateAg,[1.])
					for cc,mm in clpMemIDAg: 
						a=len(mm)
						c = (1.-numpy.random.normal(LMA1EllipsoidalVoidMean[ii][0],LMA1EllipsoidalVoidSD[ii][0])) *a  #AgRanRemovPerc*a
						iij=0;te=0
						n=int(c)
						while iij<n:
							b=mm[random.randint(0,a-1)]
							bb=O.bodies[b]
							if bb: 
								O.bodies.erase(b)
								iij+=1
								te+=1
						#print("a",a,"c",c,"te",te)

					keyyss=list(O.bodies[clpMemIDAg[0][0]].shape.members.keys())
					for k in keyyss:
						if O.bodies[k]:
							_Volume+=4./3.*math.pi*(O.bodies[k].shape.radius)**3.

					keyyss=list(O.bodies[clpMemIDAg[0][0]].shape.members.keys())
					AgVolume=_Volume
					TLAgVol+=AgVolume
					#print("TLAgVol",TLAgVol,"LMA1VolFracAg[ii][0]*Requested_LMA1_Volume",LMA1VolFracAg[ii][0]*Requested_LMA1_Volume)
					ret.append([AgVolume,"LMA1",clpMemIDAg])			
				if Status==2 or GCSC==0:
					#if Rr<1.e-4:
						AgVolume=V_known
						TLAgVol+=AgVolume
						#print("TLAgVol",TLAgVol)
						AgStandAlone.append(SPp)				
			VolSoFar+=TLAgVol
	return VolSoFar,ret				
				

######################################################################
#This function calls "polyhedra2stl" and "Filler" and scale each template based on the size distribution for pl grains
#Returns are an array for the generated spheres and the volume of the grain
######################################################################
#### GRC3:

def GRC3Generator(VolSoFar):
	cr=0;ijk=0
	global crsplSmall; global ccsplSmall; global TLplVol; global plco;idspl=[]; global clpMemIDpl;clpMemIDpl=[];retpl=[]; SpheresData=[];keyyss=[]
	TLplVol=0;plco=0; global spMannualGRC3; spMannualGRC3=[]; global plStandAlone; plStandAlone=[]
	for ii in range(len(GRC3psd)-1):
		if GRC3VolFrac[ii][0]>0:
			TLplVol=0.
			while TLplVol<GRC3VolFrac[ii][0]*Requested_GRC3_Volume:
				SpheresData=[];keyyss=[];_Volume=0.
				ijk+=1
			#sp.makeCloud(minCorner=min_corner,maxCorner=max_corner,rMean=((psd[ii]+psd[ii+1])/4.),rRelFuzz=((psd[ii]-psd[ii+1])/4./((psd[ii]+psd[ii+1])/4.)),num=1)#
			#Rr=[r for c,r in sp]
				Rr=random.uniform(GRC3psd[ii]/2./1.e6, GRC3psd[ii+1]/2./1.e6)
				SPp=utils.sphere(RandomLocation(),Rr)
				spMannual.append(SPp)
				spMannualGRC3.append(SPp)
				V_known=4./3.*math.pi*Rr**3
			#idspl.append(O.bodies.append([utils.sphere(RandomLocation(),Rr[-1])]))
				print("The ",ijk,"th GRC3 Grain with Diameter ",Rr*1e6*2., " microns is being generated from the ",ii,"th psd")
				c_a= numpy.random.normal(GRC3LenRatioMean[ii][0],GRC3LenRatioSD[ii][0])
				if c_a< (GRC3LenRatioMean[ii][0]-GRC3LenRatioSD[ii][0]):
					c_a=(GRC3LenRatioMean[ii][0]-GRC3LenRatioSD[ii][0])
				if c_a> (GRC3LenRatioMean[ii][0]+GRC3LenRatioSD[ii][0]):
					c_a=(GRC3LenRatioMean[ii][0]+GRC3LenRatioSD[ii][0])			
				c_b= numpy.random.normal(GRC3LenRatioMean[ii][0],GRC3LenRatioSD[ii][0])
				if c_b< (GRC3LenRatioMean[ii][0]-GRC3LenRatioSD[ii][0]):
					c_b=(GRC3LenRatioMean[ii][0]-GRC3LenRatioSD[ii][0])
				if c_b> (GRC3LenRatioMean[ii][0]+GRC3LenRatioSD[ii][0]):
					c_b=(GRC3LenRatioMean[ii][0]+GRC3LenRatioSD[ii][0])
					
				plc=2.*Rr*1.e5 #should be diameter
				pla = plc/c_a
				plb = pla*c_b			
				TempRad=0.00009/2.*1.e5#min(plc,plb)/5.#4.*plCutoff*1.e4
				#print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++pla,plb,plc",pla,plb,plc,TempRad)
				Status=0
				if GCSC!=0:				
					PolVolume=polyhedra2stl('GRC3.stl',pla,max(plc,plb),min(plc,plb),0.3)
					(Status,Returned)=Filler('GRC3.stl',TempRad,plc,plb,PolVolume=PolVolume)
				if (Status==1 and GCSC!=0):
					_Volume=0.
					O.bodies.append(SPp)
					clpMemIDpl=O.bodies.replaceByClumps(templatepl,[1.])
					keyyss=list(O.bodies[clpMemIDpl[0][0]].shape.members.keys())
					for k in keyyss:
						_Volume+=4./3.*math.pi*(O.bodies[k].shape.radius)**3.		
					plVolume=_Volume
					#print("TLplVol",TLplVol,"GRC3VolFrac[ii][0]*Requested_GRC3_Volume",GRC3VolFrac[ii][0]*Requested_GRC3_Volume)					
					retpl.append([plVolume,"GRC3",clpMemIDpl])
					TLplVol+=plVolume
				if (Status==2 or GCSC==0):
					#if Rr<1.e-4:
						plVolume=V_known
						TLplVol+=plVolume
						#print("TLplVol",TLplVol,"GRC3VolFrac[ii][0]*Requested_GRC3_Volume",GRC3VolFrac[ii][0]*Requested_GRC3_Volume)					
						plStandAlone.append(SPp)				

			VolSoFar+=TLplVol
	return VolSoFar,retpl 

######################################################################
def SmallGrainGenerator():
 global MinSSRGRC3; global MaxSSRGRC3; global MinSSRLMA1; global MaxSSRLMA1; 
 MinSSRGRC3=1e6
 MaxSSRGRC3=-1e6
 MinSSRLMA1=1e6
 MaxSSRLMA1=-1e6
 for i in range(len(AgStandAlone)):
  AgStandAlone[i].GrainID =0#
  AgStandAlone[i].ClumpID =0#
  AgStandAlone[i].vol = 0.#
  AgStandAlone[i].GrainType = "LMA1" # tell each particle what is its grain type they belong to
  AgStandAlone[i].GrainSplit =0
  AP=O.bodies.append([AgStandAlone[i]])#
  AgStandAlone[i].GrainID =AP[0]#tell each particle who is its agglomerate
  AgStandAlone[i].ClumpID =AP[0]#tell each particle who is its agglomerate
  AgStandAlone[i].vol = 4./3.*math.pi*(O.bodies[AP[0]].shape.radius)**3.#
  print("LMA1 Small Grains",O.bodies[AP[0]].shape.radius)
  MinSSRLMA1=min(MinSSRLMA1,O.bodies[AP[0]].shape.radius)
  MaxSSRLMA1=max(MaxSSRLMA1,O.bodies[AP[0]].shape.radius) 
  
 for i in range(len(plStandAlone)):
  plStandAlone[i].vol = 0.#
  plStandAlone[i].GrainID =0#tell each particle who is its agglomerate
  plStandAlone[i].ClumpID =0#
  plStandAlone[i].GrainType = "GRC3" # tell each particle what is its grain type they belong to
  plStandAlone[i].GrainSplit =0
  AP=O.bodies.append([plStandAlone[i]])#
  plStandAlone[i].GrainID =AP[0]#tell each particle who is its agglomerate
  plStandAlone[i].ClumpID =AP[0]#
  plStandAlone[i].vol = 4./3.*math.pi*(O.bodies[AP[0]].shape.radius)**3.#
  print("GRC3 Small Grains",O.bodies[AP[0]].shape.radius)
  MinSSRGRC3=min(MinSSRGRC3,O.bodies[AP[0]].shape.radius)
  MaxSSRGRC3=max(MaxSSRGRC3,O.bodies[AP[0]].shape.radius) 
######################################################################
#This function replaces all clump members with spheres along with the volume of each grain, grain type, and grain Id that the spheres belong to
#We do not call this function until the sample is compacted because the goal is to benefit the reduced computation time when working with clumps
def releaseAmid(ret=[]):
	i=0;ImagesSign=0; 
	for i in range(len(ret)):
		for j in ret[i][2][0][1]:#first is row, second is clump id, third is members
			if O.bodies[j] and isinstance(O.bodies[j].shape,Sphere) and O.bodies[j].isClumpMember:
				ppo=O.bodies[j].state.pos;pra=O.bodies[j].shape.radius
				O.bodies.erase(O.bodies[j].id)
				SIds=O.bodies.append([utils.sphere(ppo,pra)])
				O.bodies[SIds[0]].GrainID = ret[i][2][0][0] # tell each particle what is its inital volume of the grain they belong to
				O.bodies[SIds[0]].ClumpID = O.bodies[SIds[0]].GrainID # tell each particle who is its agglomerate
				O.bodies[SIds[0]].vol = ret[i][0] # tell each particle what is its inital volume of the grain they belong to
				O.bodies[SIds[0]].GrainType = ret[i][1] # tell each particle what is its grain type they belong to
				O.bodies[SIds[0]].GrainSplit = 100 # We give a random value here but will update it later 
#we here mark the grains that are split by the PB
#this function informs all duplicated grains regarding their attributres

######################################################################
#Here Ag, pl, brA, brB grains are generated with size above the cutoff sizes are generated
######################################################################

VolSoFar,AAg=AgGenerator()



f= open("LMA1SphereDiameter.txt", 'w')
hi=0
for spm in spMannualLMA1:
	hi+=1
	Srad=spm.shape.radius*1.e6
	f.write(' %f' % Srad)
	f.write('\n')
f.close()




VolSoFar,AGRC3=GRC3Generator(VolSoFar) 

SmallGrainGenerator()

_Volume=0;Nin=0
for b in O.bodies:
	if b.isClump:#(b.shape,Sphere):
		Nin+=1
for b in O.bodies:
	if isinstance(b.shape,Sphere):
		#Nin+=1
		_Volume+=4./3.*math.pi*(b.shape.radius)**3.
print("Total Number of Grains is",Nin, "that should be",TIGN)
print("Total Volume of Grains is", _Volume, "while", RequestedSampleVolume," was requested")	

SpheresNo=0
ClumpNo=0
ClumpMemberNo=0
StandAlonesphere=0
MNM=0

for b in O.bodies:
	if isinstance(b.shape,Sphere):
		SpheresNo+=1

for b in O.bodies:
	if b.isClump:
		NM=0
		ClumpNo+=1
		keyyss=list(b.shape.members.keys())
		for k in keyyss:
			NM+=1
		MNM=max(NM,MNM)
	if b.isClumpMember:
		ClumpMemberNo+=1
for b in O.bodies:
	if not b.isClumpMember and not b.isClump:
		StandAlonesphere+=1
print('ClumpNo',ClumpNo,'ClumpMemberNo',ClumpMemberNo,'StandAlonesphere',StandAlonesphere,"MNM",MNM,"SpheresNo",SpheresNo)


maxR=max([b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]) 
minR=min([b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)]) 
RR=maxR/minR
print("maxR",maxR,"minR",minR,"RR",RR,"MinSSRGRC3,MaxSSRGRC3,MinSSRLMA1,MaxSSRLMA1",MinSSRGRC3,MaxSSRGRC3,MinSSRLMA1,MaxSSRLMA1)







#releasing the grains and storing their attributes to their memory
releaseAmid(AAg)
releaseAmid(AGRC3)




TGLIST=[];TGVolume=0.
for b in O.bodies:
	if isinstance(b.shape,Sphere):
		if b.GrainID in TGLIST:continue
		TGVolume+=float(b.vol)
		TGLIST.append(b.GrainID)		
print('Total Volume of Grains using clump s recorded volumes',TGVolume, "that should be", _Volume, "while", RequestedSampleVolume," was requested")	



##Here, the loose clouds of GRC-3 or GRC-3+LMA-1 are exported the final phases of sample preparation including alignment and relaxation
export.text("/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Typical Samples/Sample_GRC-3PlusLMA-1_Triaxial_OutwardPSD_CLumps_GCSC2p5_TestA_Initial.txt")
#textExtAPriodic('/home/ngoudarzi/Desktop/GRC-3 Calibration/Outward PSD/Triaxial/Sample Generation/RepCloud_Triaxial_PSD75-500_Spherical_n0337&n0365.txt','x_y_z_r_attrs_vol_GrainType_GrainSplit_GrainID',attrs=['b.GrainID'],vol=['b.vol'],GrainType=['b.GrainType'],GrainSplit=['b.GrainSplit'],GrainID=['b.GrainID'])

f= open("SphereDiameter.txt", 'w')
hi=0
for spm in spMannual:
	hi+=1
	f.write(' %f' % spm.shape.radius)
	f.write('\n')
f.close()



f= open("GRC3SphereDiameter.txt", 'w')
hi=0
for spm in spMannualGRC3:
	hi+=1
	f.write(' %f' % spm.shape.radius)
	f.write('\n')
f.close()
# ret




# #plotting



# import matplotlib.pyplot as plt

# # An "interface" to matplotlib.axes.Axes.hist() method
# n, bins, patches = plt.hist(x=d, bins='auto', color='#0504aa',
#                             alpha=0.7, rwidth=0.85)
# plt.grid(axis='y', alpha=0.75)
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# plt.title('My Very Own Histogram')
# plt.text(23, 45, r'$\mu=15, b=3$')
# maxfreq = n.max()
# # Set a clean upper y-axis limit.
# plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


