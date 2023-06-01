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
######################################################################
LMA1fileName="/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Initial Cloud Generation/LMA1StatiscalAnalysesTest.csv"
GRC3fileName="/home/ngoudarzi/Desktop/Triaxial Test/GRC-3&GRC-3A Sample Generation Paradigm/Initial Cloud Generation/GRC3.txt"
#################################################
# Number of bins for obtaining psd for each simulant type. This should be between 10 to 50. The higher values make the shape and size distribution closer to that of the
# input file as long as a large volume of sample is being requested. Small samples less than the volume of the input grains should have lower range of bin numbers.
# This should be 500 or above if generation of the exact input grain file is desired. Otherwise, a value between 20 to 50; otherwise, there might not be enough number of values
# in each bin to generate statisticaally meaningful distribution (size and shape) for each bin.
GRC3VFs=50
LMA1VFs=50
#################################################
# Min and Max size for grains
LMA1MinSize=90.
LMA1MaxSize=800.
GRC3MinSize=90.
GRC3MaxSize=800.
#################################################
# This is Grain Constitute Sphere Controller (GCSC) with recommanded value between 2.5 to 10.
# Increasing this value increases the runtime but reduces grain roughness and thus makes the mechanical response more brittle with higher peak strength
# The optimum number is 4 from parametric analyses in both runtime and triaxial peak strength.
# The fastest case can be obtained by assigning 2.5 to GCSC. GCSC=0 gives samples made of single-sphere grains
GCSC=2.5		
######################################################################
# This is the volume of the sample you would like to generate. 
RequestedSampleVolume=0.0000000504# Solid Volume related to the desired porosity for the selected sample sizes
#################################################
# This is the ratio of LMA1 mass to the total mass of the requested sample
LMA1_Mass_To_Total_Mass_Percentage=0.5
#################################################
Density_GRC3=2.63 # gr/cm3
Density_LMA1=2.9  # gr/cm3

# Grain microproperties
######################################################################
intR=1.5 # allows near neighbour interaction (can be adjusted for every packing); this is bond Generator Factor defining the gap between spheres in ordeer for them to be bonded. This is later adjusted to 1 to allow only interaction between grains that are touching
factor=intR
DENS=3100. # could be adapted to match material density: dens_DEM=dens_rock*(V_rock/V_particles)=dens_rock*1/(1-poro_DEM) -> packing porosity as to be computed?
YOUNG=2.e9
FRICTShallow=42.# in degrees
FRICTDeep=54.# in degrees
Poisson=0.1
TENS=1.e6
COH=1.e6
#Unbalance force limit:
UnBalLimit=5.e-5




