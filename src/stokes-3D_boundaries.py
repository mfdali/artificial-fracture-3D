#--------------------------
# Solve flow through artificial 3D fractures
# Stokes equations with an iterative solver 
# @author: Monique Dali @lmmp-puc-rio
# project: fractured medium
#------------------------------

# Load libraries
from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import time

## Functions to keep mesh parameters

def slope():
    # slope surface
    case_filename = "x2_slopev1_10N_2mm_meshsize_0002"

    #Geometry
    x_ = 0.05 # Lenght
    y_ = 0.02 # Width
    h_ = 0.0009 # Minimum
    h_avg = 1.5000e-03 # z average

    return case_filename,x_,y_,h_,h_avg


def eggshell():
    # Periodic surface
    case_filename = "x2_reducedv2_1mm_meshsize_00035"

    #Geometry
    x_ = 0.05 # Lenght
    y_ = 0.02 # Width
    h_ = 0.0009 # Minimum
    h_avg = 1.000e-03 # z average

    return case_filename,x_,y_,h_,h_avg

# Use perdiodic mesh
case_filename,x_,y_,h_,h_avg = eggshell()

# Use perdiodic mesh with slope
#case_filename,x_,y_,h_,h_avg = slope()

# Test for PETSc or Tpetra
if not has_linear_algebra_backend("PETSc") and not has_linear_algebra_backend("Tpetra"):
    info("DOLFIN has not been configured with Trilinos or PETSc. Exiting.")
    exit()

if not has_krylov_solver_preconditioner("amg"):
    info("Sorry, this demo is only available when DOLFIN is compiled with AMG "
         "preconditioner, Hypre or ML.")
    exit()

if has_krylov_solver_method("minres"):
    krylov_method = "minres"
    solver_method = 'minres'
elif has_krylov_solver_method("tfqmr"):
    krylov_method = "tfqmr"
    solver_method = 'tfqmr'
else:
    info("Default linear algebra backend was not compiled with MINRES or TFQMR "
         "Krylov subspace method. Terminating.")
    exit()

#############
## Path
path_to_folder = "meshes/"
output_folder = "results/"

# name of the project
project_name = "test_3d_meshes"

# Create dictionary for the boundaries
boundary_tags = {'surface':1,'back':2,'bottom':3,'outlet':4,'top':5,'inlet':6}

################# Mesh #######################
# Load mesh
mesh = Mesh()
with XDMFFile(path_to_folder + case_filename + ".xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()) 

# Load boundary mesh
with XDMFFile(path_to_folder + case_filename + "_boundaries.xdmf") as infile:
        infile.read(mvc, "name_to_read")

# Mark boundaries
boundaries = cpp.mesh.MeshFunctionSizet(mesh, mvc)

# Mark mesh domain
domains = MeshFunction("size_t", mesh, mesh.topology().dim())

dx = Measure("dx",domain=mesh)#, subdomain_data=domains) #or mf
ds = Measure('ds',domain=mesh, subdomain_data = boundaries)

# Save boundaries to check in paraview
file = File(output_folder + case_filename + "_boundaries.pvd")
file << boundaries

################################
# Fluid properties

mu = Constant(0.001) #Water viscosity [Pa.s]
mu_ = 7.5*mu #R.C. Givler, S.A. Altobelli, J. Fluid Mech. 258 (1994) 355.
pin = Constant(2*6894.76)#*x_ #Imposed pressure at the entrance [Pa]
pout = Constant(6894.76)#*x_ #Imposed pressure on output [Pa]

####### Numerical Solution #######
start_time = time.time()

#Function space over the mesh
V = VectorElement('CG',mesh.ufl_cell(),2)
Q = FiniteElement('CG',mesh.ufl_cell(),1)
Element = V*Q
W = FunctionSpace(mesh,Element)

#Define variational problem
(u,p) = TrialFunctions(W)
(v,q) = TestFunctions(W)

# Print degrees of freedom
info("Num DOFs {}".format(W.dim()))
##############################################################
# Parameters

# Define expressions used in variational forms

#pressure gradient
dp = pin-pout

# No-slip boundary condition for velocity
noslip = Constant((0.0,0.0,0.0)) #No-slip condition for velocity, u=0 at y=h
f = Constant((0.0,0.0,0.0)) #External force
n = FacetNormal(mesh) #Normal vector to mesh

# Inflow boundary condition for velocity
inflow = Expression(("sin((x[1]*pi)/w)", "0.0", "0.0"),w = y_, degree=2)

##############################################################
 
# Define Dirichlet boundary conditions
bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,boundary_tags['inlet'])
bc2 = DirichletBC(W.sub(0),noslip,boundaries,boundary_tags['top'])
bc3 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,boundary_tags['outlet'])
bc4 = DirichletBC(W.sub(0),noslip,boundaries,boundary_tags['bottom'])
bc5 = DirichletBC(W.sub(0),noslip,boundaries,boundary_tags['surface'])
bc6 = DirichletBC(W.sub(0),noslip,boundaries,boundary_tags['back'])

# Collect boundary conditions
bcs = [bc1,bc2,bc3,bc4,bc5,bc6]

##############################################################

# Weak form for Stokes flow
F = (mu*inner(grad(u),grad(v))*dx - div(v)*p*dx - div(u)*q*dx) - inner(f,v)*dx + pin*dot(v,n)*ds(boundary_tags['inlet']) + pout*dot(v,n)*ds(boundary_tags['outlet']) 

a = lhs(F)
L = rhs(F)

#Compute solution
U = Function(W, name = "field")

# Form for use in constructing preconditioner matrix
b = mu*inner(grad(u), grad(v))*dx + p*q*dx

# Assemble system
A, bb = assemble_system(a, L, bcs)

# Assemble preconditioner system
P, btmp = assemble_system(b, L, bcs)

# Create Krylov solver and AMG preconditioner
solver = KrylovSolver(krylov_method, "amg")

# Associate operator (A) and preconditioner matrix (P)
solver.set_operators(A, P)

#Solve
solver.solve(U.vector(),bb)

##############################################################

#Get sub functions
(u, p) = U.split()

######### Calculate metrics #######

# Calculate inlet flow
form = -dot(n, u) * ds(boundary_tags['inlet'])
inflow = assemble(form)

# Calculate outlet flow
form = dot(n, u) * ds(boundary_tags['outlet'])
outflow = assemble(form)
deviation = ((inflow-outflow)/inflow)*100
print("Flow: %.10g\t %.10g\t deviation: %e\n" % (inflow,outflow,deviation))

# Get execution time
exec_time = time.time() - start_time
print(exec_time)

# Save solution in VTK format
file = File(output_folder + case_filename + "-u.pvd")
file << u
file = File(output_folder + case_filename + "-p.pvd")
file << p

# Measure height
# Define the Z-coordinate function
z = Expression('x[2]', degree=1)

# Calculate the total area of the surface
Volume = float(assemble(1 * ds(boundary_tags['surface'])))  # Surface area of the mesh

# Calculate the integral of the Z-coordinate over the surface mesh
z_integral = assemble(z * ds(boundary_tags['surface']))

# Calculate the average Z-coordinate
int_H2 = z_integral/Volume
print("a_h_average = %e" % int_H2)

# Average from lubrification theory
a_H = assemble((1/(h_avg**3))*dx)
#Area = float(assemble(1 * dx))
#int_H = pow((Area/a_H),(1/3))
#print("a_h_average_lub = %e" % int_H)

# mesh size
triangles = mesh.num_cells()
print("%g\n" % triangles)

# degrees of freedom
dofs = W.dim()

# Calculate outflow area
outflow_area = float(assemble(1*ds(boundary_tags['outlet'])))

# Calculate homogeneous permeability
k_eq = (outflow * mu * x_)/(outflow_area*dp)
# Convert to darcy
k_D = k_eq*1013250000000
# Calculate permeability with the cubic law
k_f = pow(int_H2,2)/12

print("Mesh: %s\n" % case_filename)
print("keq(m²):%.5g\t keq(D):%g\n" % (k_eq,k_D))

# Calculate permeability deviation betweeen cubic law and numerical simulation
kf_dev = (abs(k_f-k_eq)/k_f)*100
print("kf(m²):%.5g\t deviation: %e\n" % (k_f,kf_dev))

### Save output ###

# Save simulation parameters and results to file
# open file
print_to_file = open(output_folder + project_name + ".txt", "a")
# header
print_to_file.write("mesh,triangles,dofs,solver,avg_ah,x,y,z,inflow,outflow,deviation,k_eq,k_f,exec_time,Area,ah_lub\n")
# results
print_to_file.write("%s,%g,%10.1g,%s,%e,%e,%e,%e,%.10g,%.10g,%e,%.5g,%.5g,%g,%g,%e\n" % (case_filename,triangles,dofs,solver_method,int_H2,x_,y_,h_,inflow,outflow,deviation,k_eq,k_f,exec_time,Volume,a_H))
# close file
print_to_file.close()