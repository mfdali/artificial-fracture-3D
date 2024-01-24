#------------------------------
# Convert mesh with meshio to hdf5,xdmf
# @author: monique @mfdali based on Stelio @steliohlopes
# project: fractured medium
# company: @lmmp-puc-rio
#------------------------------

# Load requirements
import meshio
import numpy as np

# Insert the name of the msh file 
msh_file = "x2_slopev1_10N_2mm_meshsize_0003" #Filename

# Read gmsh file
mesh = meshio.read(meshFile + ".msh")

# Create empty list
tetra_cells = None
tetra_data = None
triangle_cells = None
triangle_data = None

# Get physical entities in msh file
for key in mesh.cell_data_dict["gmsh:physical"].keys():
    if key == "triangle":
        triangle_data = mesh.cell_data_dict["gmsh:physical"][key]
    elif key == "tetra":
        print("______________ found: tetra_data")
        tetra_data = mesh.cell_data_dict["gmsh:physical"][key]

# Get tetras (from 3D mesh) and triangles
for cell in mesh.cells:
    if cell.type == "tetra":
        if tetra_cells is None:
            tetra_cells = cell.data
        else:
            tetra_cells = np.vstack([tetra_cells, cell.data])
    elif cell.type == "triangle":
        if triangle_cells is None:
            triangle_cells = cell.data
        else:
            triangle_cells = np.vstack([triangle_cells, cell.data])

# Convert 3d mesh
tetra_mesh = meshio.Mesh(points=mesh.points, cells={"tetra": tetra_cells},
                         cell_data={"name_to_read":[tetra_data]})
# Convert 2d boundaries
triangle_mesh =meshio.Mesh(points=mesh.points,
                           cells=[("triangle", triangle_cells)],
                           cell_data={"name_to_read":[triangle_data]})

#Save mesh
meshio.write(meshPath + meshFile +".xdmf", tetra_mesh)

#Save boundaries
meshio.write(meshPath + meshFile +"_boundaries.xdmf", triangle_mesh)