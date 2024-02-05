# -----------------------------------------------------------------------------
#
#  Mesh generation with Gmsh 
#  Create artificial fracture based on terrain meshing - x2.py Python extended tutorial 2
#   
#  author: monique @mfdali
#  company: @lmmp-puc-rio
#  project: fractured medium
#
#  Mesh import, discrete entities, hybrid models, terrain meshing
#
# -----------------------------------------------------------------------------

# Load libraries
import gmsh
import sys
import math
import random
import numpy as np

# Surface functions
def slope(min_meshsize):
    # Periodic mesh with slope
    '''
    #Geometry
    x_ = 0.05 # Lenght
    y_ = 0.02 # Width
    h_ = 0.0009
    tol = 1E-6
    h_avg = 1.5000e-03
    '''
    # Define name of the surface type
    surface_type = 'slope'

    # Set the name of the msh file
    mesh_name = "x2_slopev1_10N_2mm_meshsize_" + str(min_meshsize)[2:]
    # Set the surface equation
    surface_equation = "1E-4 * math.sin(10 * float(i+j) / 10*N) + 0.002 - ((10*float(i))/(10000*M))"
    
    return surface_type, mesh_name, surface_equation


def eggshell(min_meshsize):
    # Periodic mesh
    '''
    #Geometry
    case_filename = "x2_reducedv1_1mm_meshsize_00025"
    x_ = 0.05 # Lenght
    y_ = 0.025 # Width
    h_ = 0.0009
    tol = 1E-6
    h_avg = 1.000e-03
    '''
    # Define name of the surface type
    surface_type = 'eggshell'
    # Set the name of the msh file
    mesh_name = "x2_reducedv2_1mm_meshsize_" + str(min_meshsize)[2:]
    # Set the surface equation
    surface_equation = "1E-4 * math.sin(10 * float(i+j) / 10*N) + 0.001"
    
    return surface_type, mesh_name, surface_equation

def export_metrics(get_metrics_file,case_name,mesh_name,Eq,h_min,h_max,h_avg,kf,min_meshsize,max_meshsize):
    # Save mesh properties to file
    get_metrics_file.write("%s,%s,%s,%.8g,%.8g,%.4e,%.4e,%.4e,%.8g\n" % (case_name,mesh_name,Eq,min_meshsize,max_meshsize,h_min,h_max,h_avg,kf))

def ah_metrics(z_coords,print_to_file,get_metrics_file,case_name,mesh_name,Eq,min_meshsize,max_meshsize):
    # Show mesh properties

    h_min = min(z_coords)
    h_max = max(z_coords)
    h_avg = np.average(z_coords)
    kf = np.power(h_avg,2)/12
    print("Min: %.4e\n" % h_min)
    print("Max: %.4e\n" % h_max)
    print("Avg: %.4e\n" % h_avg)
    print("Cubic law kf: %.4e\n" % kf)

    # Save mesh properties to file
    if print_to_file:
        export_metrics(get_metrics_file,case_name,mesh_name,Eq,h_min,h_max,h_avg,kf,min_meshsize,max_meshsize)
    

def create_mesh(set_meshsize, surface_mode, metrics_filename, case_name, print_to_file, print_header):

    # Set parameters
    min_meshsize = set_meshsize
    max_meshsize = min_meshsize + 0.00005
    lcmin = min_meshsize
    lcmax = max_meshsize
    lc = min_meshsize

    # Function to create the desired surface: periodic or periodic with slope 
    surface_type, mesh_name, surface_equation = surface_mode(min_meshsize)

    # Mesh generation
    gmsh.initialize()

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

    # Set seed to replicate surface
    random.seed(10)

    # create the points based on a given equation
    for i in range(M + 1):
        for j in range(N + 1):
            nodes.append(tag(i, j))
            if surface_type == 'slope':
                z_eq = 1E-4 * math.sin(10 * float(i+j) / 10*N) + 0.002 - ((10*float(i))/(10000*M))
            
            if surface_type == 'eggshell':
                z_eq = 1E-4 * math.sin(10 * float(i+j) / 10*N) + 0.001

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

    # Add corner points
    if surface_type == 'slope':
        p1 = gmsh.model.geo.addPoint(0, 0, 0,lcmax)
        p2 = gmsh.model.geo.addPoint(1/20, 0, 0,lcmin)
        p3 = gmsh.model.geo.addPoint(1/20, 1/50, 0,lcmin)
        p4 = gmsh.model.geo.addPoint(0, 1/50, 0,lcmax)

    if surface_type == 'eggshell':
        p1 = gmsh.model.geo.addPoint(0, 0, 0)
        p2 = gmsh.model.geo.addPoint(1/20, 0, 0)
        p3 = gmsh.model.geo.addPoint(1/20, 1/50, 0)
        p4 = gmsh.model.geo.addPoint(0, 1/50, 0)

    # Connect these points
    c1 = gmsh.model.geo.addLine(p1, p2)
    c2 = gmsh.model.geo.addLine(p2, p3)
    c3 = gmsh.model.geo.addLine(p3, p4)
    c4 = gmsh.model.geo.addLine(p4, p1)
    c10 = gmsh.model.geo.addLine(p1, 1)
    c11 = gmsh.model.geo.addLine(p2, 2)
    c12 = gmsh.model.geo.addLine(p3, 3)
    c13 = gmsh.model.geo.addLine(p4, 4)

    # Connect the lines except surface
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

    # Connect the boundaries to create a closed volume 
    sl1 = gmsh.model.geo.addSurfaceLoop([s1, s3, s4, s5, s6, 1])

    # Create volume
    v1 = gmsh.model.geo.addVolume([sl1])
    gmsh.model.geo.synchronize()

    # Get boundaries
    boundary = gmsh.model.getBoundary(gmsh.model.getEntities(3))
    boundary_ids = [b[1] for b in boundary]
    # Add a dict containing labels for each boundary
    #boundary_tags = ['surface','back','bottom','outlet','top','inlet']

    print(" - boundary entities " + str(boundary))
    # Set tags to each boundary
    for i in range(len(boundary_ids)):
        b = boundary_ids[i]
        gmsh.model.addPhysicalGroup(2, [b], tag=b)
        # FIXME: how to call the PhysicalName in FEniCS legacy 
        # Set name to each boundary
        #gmsh.model.setPhysicalName(2, b, boundary_tags[i])

    #Add physical tag 2 for the volume
    volume_entities = [model[1] for model in gmsh.model.getEntities(3)]
    gmsh.model.addPhysicalGroup(3, [v1], tag=0)
    gmsh.model.setPhysicalName(3, 2, "volume")

    if surface_type == 'slope':
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
        gmsh.model.mesh.field.setNumber(2, "SizeMin", lc / 1.5)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lcmax)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 0.002)
        gmsh.model.mesh.field.setNumber(2, "DistMax", 0.03)

        # Let's use the minimum of all the fields as the mesh size field:
        gmsh.model.mesh.field.add("Min", 7)
        gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2])

        gmsh.model.mesh.field.setAsBackgroundMesh(7)

    # The API also allows to set a global mesh size callback, which is called each
    # time the mesh size is queried
    def meshSizeCallback(dim, tag, x, y, z, lc):
        #return min(lc, 0.002 * x + 0.001)
        return min(lc, 0.0002 * x + 0.001)

    if surface_type == 'eggshell':
        gmsh.option.setNumber('Mesh.MeshSizeMin', min_meshsize)
        gmsh.option.setNumber('Mesh.MeshSizeMax', max_meshsize)

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
    gmsh.write(mesh_name + '.msh')

    gmsh.model.geo.synchronize()

    gmsh.finalize()

    # Save to file parameters of mesh creation
    if print_to_file:
        get_metrics_file = open(metrics_filename + ".txt", "a")
        # Print header
        if print_header:
            get_metrics_file.write("Case,Equation,mesh_name,min_meshsize,max_meshsize,min,max,Av_Height,kf\n")

        ah_metrics(z_coords,print_to_file,get_metrics_file,case_name,mesh_name,surface_equation,min_meshsize,max_meshsize) 

    # Create 3d array of surface points
    array_coords = np.asarray(coords)
    grid_coords = array_coords.reshape((int(len(coords)/3)),3)
    
    # Save surface coordinates to file
    np.savetxt(surface_type + '10N_coords.csv', grid_coords , fmt="%.8e", delimiter=",")  

############ Main ############

metrics_filename = "gmsh_terrain_surface2" # Filename to write parameters of mesh creation
print_to_file = True
print_header = True
set_meshsize = 0.0002 #meshsize for Gmsh
surface_mode = eggshell # eggshell or slope

case_name = "terrain_modified"

create_mesh(set_meshsize, surface_mode, metrics_filename, case_name, print_to_file, print_header)
