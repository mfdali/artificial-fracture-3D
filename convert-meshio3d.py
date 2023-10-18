import matplotlib.pyplot as plt
import meshio
import h5py
import numpy as np

#Path
output = "/mnt/g/My Drive/Fenicsx/gmsh/stl/"

msh_file = "x2_slopev1_2mm_meshsize_0004"

mesh_from_file = meshio.read(output + msh_file + ".msh")

'''
def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:, :2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={
                           "name_to_read": [cell_data]})
    return out_mesh

triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write(output + msh_file + "_boundaries.xdmf", triangle_mesh)

tetra_mesh = create_mesh(mesh_from_file, "tetra", prune_z=True)
meshio.write(output + msh_file + ".xdmf", tetra_mesh)


# Path to your Gmsh-generated .msh file
msh_file = "x2_random"

# Load the .msh mesh using meshio
mesh = meshio.read(output + msh_file + ".msh")

# Convert the mesh to the .xdmf format
xdmf_file = msh_file + ".xdmf"
with h5py.File(xdmf_file, "w") as f:
    # Create the mesh group
    mesh_group = f.create_group("Mesh")
    
    # Create the geometry dataset
    geo_dataset = mesh_group.create_dataset("Geometry", data=mesh.points, compression="gzip")

    # Create the topology dataset for cells
    topo_dataset = mesh_group.create_dataset("Topology", data=mesh.cells[0].data, compression="gzip")
    topo_dataset.attrs["TopologyType"] = mesh.cells[0].type

    # Create the datasets for the cell types
    cell_types = [t.name for t in mesh.cells]
    cell_types_dataset = mesh_group.create_dataset("CellTypes", data=cell_types, compression="gzip")
    
    # Add the attributes
    geo_dataset.attrs["Partitioned"] = 0  # Assuming it's not a partitioned mesh
    geo_dataset.attrs["GridType"] = "Curvilinear"
    geo_dataset.attrs["Number of Points"] = mesh.points.shape[0]
    topo_dataset.attrs["NumberOfElements"] = mesh.cells[0].data.shape[0]
'''
'''
meshio.write(msh_file + ".xdmf", meshio.Mesh(points=mesh_from_file.points, cells={"tetra": mesh_from_file.cells_dict["triangle"]}))
meshio.write("obstacles_facet_region.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
                                    cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
meshio.write("obstacles_physical_region.xdmf", meshio.Mesh(
    points=msh.points, cells={"tetra": msh.cells["tetra"]},
    cell_data={"tetra": {"name_to_read":
                            msh.cell_data["tetra"]["gmsh:physical"]}}))'''

tetra_cells = []
for cell in mesh_from_file.cells:
    if  cell.type == "tetra":
        if len(tetra_cells) == 0:
            tetra_cells = cell.data
        else:
            tetra_cells = np.vstack([tetra_cells, cell.data])

tetra_mesh = meshio.Mesh(points=mesh_from_file.points, cells={"tetra": tetra_cells})
meshio.write(output + msh_file + ".xdmf", tetra_mesh)