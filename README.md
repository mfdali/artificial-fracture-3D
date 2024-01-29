# Artificial fracture 3D

Python scripts to create and simulate flow through thin artificial fractures in 3D.
The goal of this project is to simulate flow through thin 3d box, $h << (x and y)$. These geometries mimic the behaviour of fractures in oil reservoirs.
We want to compare computational results for 3 different surface arrangements. For this work, we created 3 types of surface: flat, harmonic and harmonic slope.
Then, one flat surface is replaced by a harmonic one and after a slope is added.

### Modified surface

![Harmonic surface](https://github.com/mfdali/artificial-fracture-3D/blob/main/x2_reducedv3_1mm_surface.png)

### Surface manipulation

![sketch](https://github.com/mfdali/artificial-fracture-3D/blob/main/slope_sketch.png)


  1. ## Create mesh
  
  ```fracture_surface_3dmesh.py```
  
  Create desirable 3D geometry using Gmsh via Python API based on [terrain mesh demo](https://gitlab.onelab.info/gmsh/gmsh/blob/master/tutorials/python/x2.py)
  
  In this case, it was constructed an artificial fracture combining an harmonic surface at top with 1 milimiter of average height and flat walls to close the cuboid.

  ### Flat surface
  ![Box](https://github.com/mfdali/artificial-fracture-3D/blob/main/thin_box_2mm_meshsize_0007_mesh.png)

  ### Harmonic surface
  ![Harmonic pertubation](https://github.com/mfdali/artificial-fracture-3D/blob/main/x2_reducedv3_1mm_meshsize_0004_mesh.png)

  ### Surface slope
  ![Slope](https://github.com/mfdali/artificial-fracture-3D/blob/main/x2_slopev1_10N_2mm_meshsize_0003.png)
  
  2. ## Convert mesh from .msh to .xdmf, .h5
  
  ```convert-meshio3d.py```
  
  Use meshio lib to convert Gmsh file into readable mesh for FEniCS. It also creates facet region file containing boundary entities.

  ### Boundaries
  
  Harmonic
  ![Harmonic pertubation](https://github.com/mfdali/artificial-fracture-3D/blob/main/x2_reducedv3_1mm_meshsize_boundaries.png)
  
  Slope
  ![Harmonic slope](https://github.com/mfdali/artificial-fracture-3D/blob/main/x2_slopev1_10N_2mm_boundaries.png)
  

  3. ## Solve flow
  
  ```stokes.py```
  
  Code to solve Stokes flow with FEniCS package through an artificial fracture.

  ### Problem
  
  ![Boundary conditions](https://github.com/mfdali/artificial-fracture-3D/blob/main/boundary_conditions.png)
  
  4. ## Mesh test
     
  ```3D_fracture_simulation_analysis.ipynb```
  
  ![mesh_test](https://github.com/mfdali/artificial-fracture-3D/blob/main/mesh_test_periodic.png)

  5. ## Fracture analysis

     ```fracture_equation.ipynb```
     
     ### Slope
     - Aperture distribution
     
     ![Histogram](https://github.com/mfdali/artificial-fracture-3D/blob/main/slope_aperture_distribution.png)
     
     - Surface
     
     ![Slope_surface](https://github.com/mfdali/artificial-fracture-3D/blob/main/slope_surface.png)

     - Inlet & outlet curve
     
     ![inlet](https://github.com/mfdali/artificial-fracture-3D/blob/main/slope_inlet_outlet_curve.png)

     - Top curve
     
     ![top](https://github.com/mfdali/artificial-fracture-3D/blob/main/slope_top_curve.png)

     - Bottom curve
     
     ![bottom](https://github.com/mfdali/artificial-fracture-3D/blob/main/slope_bottom_curve.png)
  
  6. ## 3D-2D Comparison

     ```stokes-fracture-artificial-comparison.py```

     
