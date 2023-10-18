# -----------------------------------------------------------------------------
#
#  Gmsh Python extended tutorial 2
#
#  Mesh import, discrete entities, hybrid models, terrain meshing
#
# -----------------------------------------------------------------------------

import gmsh
import sys
import math
import random
import numpy as np

filepath = "/mnt/g/My Drive/Fenicsx/gmsh/stl/"
metrics_filename = "gmsh_terrain_surface2"
print_to_file = True
print_header = True

#case_name = "eggshell_slope"
#mesh_name = "slope_1mm_meshsize_0003"

#case_name = "node_thin_1mm"
case_name = "terrain_modified"
mesh_name = "x2_slopev1_2mm_meshsize_0004"
min_meshsize = 0.0004
max_meshsize = 0.0004
lcmin = min_meshsize
lcmax = max_meshsize
lc = max_meshsize

def export_metrics(get_metrics_file,case_name,mesh_name,Eq,h_min,h_max,h_avg,kf,min_meshsize,max_meshsize):

    get_metrics_file.write("%s,%s,%s,%.8g,%.8g,%.4e,%.4e,%.4e,%.8g\n" % (case_name,mesh_name,Eq,min_meshsize,max_meshsize,h_min,h_max,h_avg,kf))

def ah_metrics(z_coords,print_to_file,get_metrics_file,case_name,mesh_name,Eq,min_meshsize,max_meshsize):

    h_min = min(z_coords)
    h_max = max(z_coords)
    h_avg = np.average(z_coords)
    kf = np.power(h_avg,2)/12
    print("Min: %.4e\n" % h_min)
    print("Max: %.4e\n" % h_max)
    print("Avg: %.4e\n" % h_avg)
    print("Cubic law kf: %.4e\n" % kf)

    if print_to_file:
        export_metrics(get_metrics_file,case_name,mesh_name,Eq,h_min,h_max,h_avg,kf,min_meshsize,max_meshsize)
    

# The API can be used to import a mesh without reading it from a file, by
# creating nodes and elements on the fly and storing them in model
# entities. These model entities can be existing CAD entities, or can be
# discrete entities, entirely defined by the mesh.
#
# Discrete entities can be reparametrized (see `t13.py') so that they can be
# remeshed later on; and they can also be combined with built-in CAD entities to
# produce hybrid models.
#
# We combine all these features in this tutorial to perform terrain meshing,
# where the terrain is described by a discrete surface (that we then
# reparametrize) combined with a CAD representation of the underground.

gmsh.initialize()

#gmsh.model.add("x2")

# We will create the terrain surface mesh from N x N input data points:
M = 200
N = 100

# Helper function to return a node tag given two indices i and j:
def tag(i, j):
    return (N + 1) * i + j + 1

# The x, y, z coordinates of all the nodes:
coords = []

# The tags of the corresponding nodes:
nodes = []

# The connectivities of the triangle elements (3 node tags per triangle) on the
# terrain surface:
tris = []
z_coords = []
# The connectivities of the line elements on the 4 boundaries (2 node tags
# for each line element):
lin = [[], [], [], []]

# The connectivities of the point elements on the 4 corners (1 node tag for each
# point element):
pnt = [tag(0, 0), tag(M, 0), tag(M, N), tag(0, N)]

#2E-5*sin(50 * x[1] * pi)*sin(30 * x[0] * pi) - (x[0]/2E3) + 1E-4

random.seed(10)

#Eq = "5E-4 * math.sin(10 * float(i+j) / 10*N) + 0.001 - (random.randrange(1,5,1)/10000)"
#Eq = "3E-4 * math.sin(10 * float(i+j) / 10*N) + 0.0007"
#Eq = "5E-4 * math.sin(80*float(j) / M) * math.sin(15 * float(i) / N) + 0.0015"
#Eq = "1E-4 * math.sin(10 * float(i+j) / N) + 0.001"
Eq = "1E-4 * math.sin(10 * float(i+j) / 9*N) + 0.002 - ((10*float(i))/(10000*M))"

for i in range(M + 1):
    for j in range(N + 1):
        nodes.append(tag(i, j))
        #z_eq = 1E-4 * math.sin(10 * float(i+j) / 10*N) + 0.001 #- (random.randrange(1,5,1)/10000)
        #z_eq = 3E-4 * math.sin(10 * float(i+j) / 10*N) + 0.0007
        #z_eq = 5E-4 * math.sin(80*float(j) / M) * math.sin(15 * float(i) / N) + 0.0015 - ((5*float(i))/(10000*M))
        z_eq = 1E-4 * math.sin(5 * float(i+j) / 10*N) + 0.002 - ((10*float(i))/(10000*M))
        coords.extend([
            (float(i) / M)/20,
            (float(j) / N)/50, z_eq
        ])
        z_coords.append(z_eq)
        if i > 0 and j > 0:
            tris.extend([tag(i - 1, j - 1), tag(i, j - 1), tag(i - 1, j)])
            tris.extend([tag(i, j - 1), tag(i, j), tag(i - 1, j)])
        if (i == 0 or i == M) and j > 0:
            lin[3 if i == 0 else 1].extend([tag(i, j - 1), tag(i, j)])
        if (j == 0 or j == N) and i > 0:
            lin[0 if j == 0 else 2].extend([tag(i - 1, j), tag(i, j)])

'''
coords[3 * tag(0, 0) - 1] = 0
coords[3 * tag(N, 0) - 1] = -0.02720105554446849
coords[3 * tag(N, N) - 1] = 0.045647262536381385 
coords[3 * tag(0, N) - 1] = -0.02720105554446849
'''

# Create 4 discrete points for the 4 corners of the terrain surface:
for i in range(4):
    gmsh.model.addDiscreteEntity(0, i + 1)

# setCoordinates(tag, x, y, z)
# coords[3 * tag(0, 0) - 1] = 0
gmsh.model.setCoordinates(1, 0, 0, coords[3 * tag(0, 0) - 1])
# coords[3 * tag(N, 0) - 1] = -0.02720105554446849
gmsh.model.setCoordinates(2, 1/20, 0, coords[3 * tag(M, 0) - 1])
# coords[3 * tag(N, N) - 1] = 0.045647262536381385
gmsh.model.setCoordinates(3, 1/20, 1/50, coords[3 * tag(M, N) - 1])
# coords[3 * tag(0, N) - 1] = -0.02720105554446849
gmsh.model.setCoordinates(4, 0, 1/50, coords[3 * tag(0, N) - 1])


# Create 4 discrete bounding curves, with their boundary points:
for i in range(4):
    gmsh.model.addDiscreteEntity(1, i + 1, [i + 1, i + 2 if i < 3 else 1])

# Create one discrete surface, with its bounding curves:
gmsh.model.addDiscreteEntity(2, 1, [1, 2, -3, -4])

# Add all the nodes on the surface (for simplicity... see below):
gmsh.model.mesh.addNodes(2, 1, nodes, coords)

# Add point elements on the 4 points, line elements on the 4 curves, and
# triangle elements on the surface:
for i in range(4):
    # Type 15 for point elements:
    gmsh.model.mesh.addElementsByType(i + 1, 15, [], [pnt[i]])
    # Type 1 for 2-node line elements:
    gmsh.model.mesh.addElementsByType(i + 1, 1, [], lin[i])
# Type 2 for 3-node triangle elements:
gmsh.model.mesh.addElementsByType(1, 2, [], tris)

# Reclassify the nodes on the curves and the points (since we put them all on
# the surface before with `addNodes' for simplicity)
gmsh.model.mesh.reclassifyNodes()

# Create a geometry for the discrete curves and surfaces, so that we can remesh
# them later on:
gmsh.model.mesh.createGeometry()

# Note that for more complicated meshes, e.g. for on input unstructured STL
# mesh, we could use `classifySurfaces()' to automatically create the discrete
# entities and the topology; but we would then have to extract the boundaries
# afterwards.

# Create other build-in CAD entities to form one volume below the terrain
# surface. Beware that only built-in CAD entities can be hybrid, i.e. have
# discrete entities on their boundary: OpenCASCADE does not support this
# feature.

p1 = gmsh.model.geo.addPoint(0, 0, 0,lcmax)
p2 = gmsh.model.geo.addPoint(1/20, 0, 0,lcmin)
p3 = gmsh.model.geo.addPoint(1/20, 1/50, 0,lcmin)
p4 = gmsh.model.geo.addPoint(0, 1/50, 0,lcmax)

c1 = gmsh.model.geo.addLine(p1, p2)
c2 = gmsh.model.geo.addLine(p2, p3)
c3 = gmsh.model.geo.addLine(p3, p4)
c4 = gmsh.model.geo.addLine(p4, p1)
c10 = gmsh.model.geo.addLine(p1, 1)
c11 = gmsh.model.geo.addLine(p2, 2)
c12 = gmsh.model.geo.addLine(p3, 3)
c13 = gmsh.model.geo.addLine(p4, 4)

ll1 = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])
s1 = gmsh.model.geo.addPlaneSurface([ll1])
ll3 = gmsh.model.geo.addCurveLoop([c1, c11, -1, -c10])
s3 = gmsh.model.geo.addPlaneSurface([ll3])
ll4 = gmsh.model.geo.addCurveLoop([c2, c12, -2, -c11])
s4 = gmsh.model.geo.addPlaneSurface([ll4])
ll5 = gmsh.model.geo.addCurveLoop([c3, c13, 3, -c12])
s5 = gmsh.model.geo.addPlaneSurface([ll5])
ll6 = gmsh.model.geo.addCurveLoop([c4, c10, 4, -c13])
s6 = gmsh.model.geo.addPlaneSurface([ll6])

sl1 = gmsh.model.geo.addSurfaceLoop([s1, s3, s4, s5, s6, 1])

v1 = gmsh.model.geo.addVolume([sl1])
gmsh.model.geo.synchronize()

'''
# Set this to True to build a fully hex mesh:
#transfinite = True
transfinite = False
transfiniteAuto = False

if transfinite:
    NN = 30
    for c in gmsh.model.getEntities(1):
        gmsh.model.mesh.setTransfiniteCurve(c[1], NN)
    for s in gmsh.model.getEntities(2):
        gmsh.model.mesh.setTransfiniteSurface(s[1])
        gmsh.model.mesh.setRecombine(s[0], s[1])
        gmsh.model.mesh.setSmoothing(s[0], s[1], 100)
    gmsh.model.mesh.setTransfiniteVolume(v1)
elif transfiniteAuto:
    gmsh.option.setNumber('Mesh.MeshSizeMin', 0.5)
    gmsh.option.setNumber('Mesh.MeshSizeMax', 0.5)
    # setTransfiniteAutomatic() uses the sizing constraints to set the number
    # of points
    gmsh.model.mesh.setTransfiniteAutomatic()
else:
    #gmsh.option.setNumber('Mesh.MeshSizeMin', min_meshsize)
    gmsh.option.setNumber('Mesh.MeshSizeMax', max_meshsize)
'''
# Say we would like to obtain mesh elements with size lc/30 near curve 2 and
# point 5, and size lc elsewhere. To achieve this, we can use two fields:
# "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
# points 5 and on curve 2. This field returns the distance to point 5 and to
# (100 equidistant points on) curve 2.
gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "PointsList", [2,3])
gmsh.model.mesh.field.setNumbers(1, "CurvesList", [c2])
gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

# We then define a `Threshold' field, which uses the return value of the
# `Distance' field 1 in order to define a simple change in element size
# depending on the computed distances
#
# SizeMax -                     /------------------
#                              /
#                             /
#                            /
# SizeMin -o----------------/
#          |                |    |
#        Point         DistMin  DistMax

gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "InField", 1)
gmsh.model.mesh.field.setNumber(2, "SizeMin", lc/ 2)
gmsh.model.mesh.field.setNumber(2, "SizeMax", lc)
gmsh.model.mesh.field.setNumber(2, "DistMin", 0.005)
gmsh.model.mesh.field.setNumber(2, "DistMax", 0.01)

# Let's use the minimum of all the fields as the mesh size field:
gmsh.model.mesh.field.add("Min", 7)
gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2])

gmsh.model.mesh.field.setAsBackgroundMesh(7)

# The API also allows to set a global mesh size callback, which is called each
# time the mesh size is queried
def meshSizeCallback(dim, tag, x, y, z, lc):
    #return min(lc, 0.002 * x + 0.001)
    return min(lc, 0.0002 * x + 0.001)

#gmsh.model.mesh.setSizeCallback(meshSizeCallback)

# To determine the size of mesh elements, Gmsh locally computes the minimum of
#
# 1) the size of the model bounding box;
# 2) if `Mesh.MeshSizeFromPoints' is set, the mesh size specified at geometrical
#    points;
# 3) if `Mesh.MeshSizeFromCurvature' is positive, the mesh size based on
#    curvature (the value specifying the number of elements per 2 * pi rad);
# 4) the background mesh size field;
# 5) any per-entity mesh size constraint;
#
# The value can then be further modified by the mesh size callback, if any,
# before being constrained in the interval [`Mesh.MeshSizeMin',
# `Mesh.MeshSizeMax'] and multiplied by `Mesh.MeshSizeFactor'.  In addition,
# boundary mesh sizes are interpolated inside surfaces and/or volumes depending
# on the value of `Mesh.MeshSizeExtendFromBoundary' (which is set by default).
#
# When the element size is fully specified by a mesh size field (as it is in
# this example), it is thus often desirable to set

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.option.setNumber("Mesh.Algorithm", 5)
gmsh.model.mesh.optimize("Netgen")
gmsh.model.mesh.generate(3)
gmsh.write(filepath + mesh_name + '.msh')

# Launch the GUI to see the results:
'''if '-nopopup' not in sys.argv:
    gmsh.fltk.run()'''

gmsh.model.geo.synchronize()
gmsh.write(filepath + "x2.geo_unrolled")

gmsh.finalize()


if print_to_file:
    get_metrics_file = open(filepath + metrics_filename + ".txt", "a")
    
    if print_header:
        get_metrics_file.write("Case,Equation,mesh_name,min_meshsize,max_meshsize,min,max,Av_Height,kf\n")

    ah_metrics(z_coords,print_to_file,get_metrics_file,case_name,mesh_name,Eq,min_meshsize,max_meshsize) 
