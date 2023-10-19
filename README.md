# artificial-fracture-3D
Python scripts to simulate flow through thin atificial fractures in 3D

# 1st
  Create mesh
  terrain.py
  Create desirable 3D geometry using Gmsh via Python API
  In this case, it was constructed an artificial fracture combining an harmonic surface at top with 1 milimiter of average height and flat walls to close the cuboid
  
# 2nd
  Convert mesh from .msh to .xdmf, .h5
  convert-meshio3d.py
  Use meshio lib to convert Gmsh file into readable mesh for FEniCS
  
# 3rd
  Solve flow
  stokes.py 
  Code to solve Stokes flow with FEniCS package through an artificial fracture 
