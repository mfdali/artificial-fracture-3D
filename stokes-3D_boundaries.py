from dolfin import *
#from fenics_adjoint import *
import numpy as np
import math
import matplotlib.pyplot as plt
import time

#class BorderBC(SubDomain):
def inside(x, on_boundary):
    r = 1E-4 * math.sin(10 * float(x[0]+x[1]) / 10*y_) + 0.002 - ((10*float(x[0]))/(10000*x_))
    return on_boundary and (x[2]> r - 1E-4)

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
elif has_krylov_solver_method("tfqmr"):
    krylov_method = "tfqmr"
else:
    info("Default linear algebra backend was not compiled with MINRES or TFQMR "
         "Krylov subspace method. Terminating.")
    exit()

#############
## Path
path_to_folder = "/home/monique/fem-fenics/3D/"
#case_filename = "slope_1mm_meshsize_0005" 0.001 1.5006E-03

project_name = "test_3d_meshes"
solver_method = 'minres'
'''
case_filename = "thin_box_1mm_meshsize_0004"

#Geometry
x_ = 0.05 # Lenght
y_ = 0.02 # Width
h_ = 0.001
tol = 1E-6
h_avg = 0.001
'''

case_filename = "x2_slopev1_2mm_meshsize_0003"

#Geometry
x_ = 0.05 # Lenght
y_ = 0.02 # Width
h_ = 0.0009
tol = 1E-6
h_avg = 1.5000e-03

'''
#Geometry
case_filename = "x2_reducedv1_1mm_meshsize_00025"
x_ = 0.05 # Lenght
y_ = 0.025 # Width
h_ = 0.0009
tol = 1E-6
h_avg = 1.000e-03
'''
#################
mesh = Mesh()
with XDMFFile(path_to_folder + case_filename + ".xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim()) 
'''with XDMFFile(path_to_folder +"obstacles_boundaries.xdmf") as infile:
    infile.read(mvc, "name_to_read")
'''
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
'''
plt.figure()
plot(mesh)
#plot(mvc)
plt.axis('on')
plt.show()
plt.savefig(path_to_folder + case_filename + ".png", dpi=300)
print(mesh.num_cells())
'''
domains = MeshFunction("size_t", mesh, mesh.topology().dim())

dx = Measure("dx",domain=mesh, subdomain_data=domains) #or mf

# Define boundary condition
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
# Smooth surface
#front = AutoSubDomain(lambda x: x[2] > h_-tol)
# Rigid surface
#front = AutoSubDomain(lambda x: near(x[2], h_))

left = AutoSubDomain(lambda x: near(x[0], 0.0))
right = AutoSubDomain(lambda x: near(x[0], x_) )
bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
top = AutoSubDomain(lambda x: near(x[1], y_))
back = AutoSubDomain(lambda x: near(x[2], 0.0))


# Define boundary markers
#front.mark(boundaries, 6)
#borderbc = BorderBC()
#borderbc.mark(boundaries,0)
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
back.mark(boundaries, 5)

ds = Measure('ds',domain=mesh, subdomain_data = boundaries)

file = File(path_to_folder + case_filename + "_boundaries.pvd")
file << boundaries

################################
#Material properties

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

info("Num DOFs {}".format(W.dim()))
##############################################################
# Parameters

# Define expressions used in variational forms
dp = pin-pout
#g = Expression("b-(a/l)*y[0]" , degree=1 ,a = Constant(dp), b = Constant(pin), l = 0.027) #g = div(u)
#u_in = Constant((0.0)) #Initial velocity in x [m/s]
noslip = Constant((0.0,0.0,0.0)) #No-slip condition for velocity, u=0 at y=h
f = Constant((0.0,0.0,0.0)) #External force
n = FacetNormal(mesh) #Normal vector to mesh

# Inflow boundary condition for velocity
inflow = Expression(("sin((x[1]*pi)/w)", "0.0", "0.0"),w = y_, degree=2)

##############################################################

#Define Dirichlet boundary conditions
#bc1 = DirichletBC(W.sub(0),inflow,boundaries,1)
bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,1)
bc2 = DirichletBC(W.sub(0),noslip,boundaries,2)
bc3 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,3)
bc4 = DirichletBC(W.sub(0),noslip,boundaries,4)
bc5 = DirichletBC(W.sub(0),noslip,boundaries,5)
bc6 = DirichletBC(W.sub(0),noslip,inside)

bcs = [bc1,bc2,bc3,bc4,bc5,bc6]
'''
#Define measures associated with the boundaries and holes
ds = Measure('ds',domain=mesh, subdomain_data = boundaries)
dx = Measure('dx',domain=mesh)
'''

##############################################################
# Compensate the drag at glass walls (3D)
#Drag_force = -(12/(h_**2))*mu*inner(u,v)*dx

# Weak form for Stokes flow
F = (mu*inner(grad(u),grad(v))*dx - div(v)*p*dx - div(u)*q*dx) - inner(f,v)*dx + pin*dot(v,n)*ds(1) + pout*dot(v,n)*ds(3)  #- Drag_force
#(a, L) = system(F)
a = lhs(F)
L = rhs(F)
'''
a = (mu*inner(grad(u),grad(v))*dx + div(v)*p*dx + div(u)*q*dx) + pin*dot(v,n)*ds(1) + pout*dot(v,n)*ds(3) 
L = inner(f,v)*dx 
'''
#Compute solution
U = Function(W, name = "field")
#solve(a,U.vector(),L)

#problem = LinearVariationalProblem(a, L, U, bcs)
#solver = LinearVariationalSolver(problem)
#solver.solve()
'''
A = assemble(a)
b = assemble(L)
[bc.apply(A,b) for bc in bcs]
#[bc.apply(U.vector()) for bc in bcs]
solve(A,U.vector(),b,solver_method)
'''
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
#ux, uy = u.split(deepcopy=True)

# inlet flow
form = -dot(n, u) * ds(1)
inflow = assemble(form)

# outlet flow
form = dot(n, u) * ds(3)
outflow = assemble(form)
deviation = ((inflow-outflow)/inflow)*100
print("Flow: %.10g\t %.10g\t deviation: %e\n" % (inflow,outflow,deviation))

exec_time = time.time() - start_time

print(exec_time)
'''
#Save data to file
#data_file.write("%s,%s,%s,%g,%g,%g,%g,%2.10g,%2.10g,%g\n" % (mode,mesh_type,sm,triangles,exec_time,k,ps,inflow,outflow,k_medium))

'''
# Save solution in VTK format
file = File(path_to_folder + case_filename + "_morep-u.pvd")
file << u
file = File(path_to_folder + case_filename + "_morep-p.pvd")
file << p

# Measure height
# Define the Z-coordinate function
z = Expression('x[2]', degree=1)

# Calculate the total volume of the mesh
Volume = float(assemble(1 * ds(0)))  # Total volume of the mesh

# Calculate the integral of the Z-coordinate over the mesh
z_integral = assemble(z * ds(0))

# Calculate the average Z-coordinate
int_H2 = z_integral/Volume
print("a_h_average = %e" % int_H2)

a_H = assemble((1/(h_avg**3))*dx)
#Area = float(assemble(1 * dx))
#int_H = pow((Area/a_H),(1/3))
#print("a_h_average_lub = %e" % int_H)

#mesh size
triangles = mesh.num_cells()
print("%g\n" % triangles)
#info("Num DOFs {}".format(W.dim()))
dofs = W.dim()

outflow_area = float(assemble(1*ds(3)))

k_eq = (outflow * mu * x_)/(outflow_area*dp)
k_D = k_eq*1013250000000

k_f = pow(int_H2,2)/12
print("keq(m²):%.5g\t keq(D):%g\n" % (k_eq,k_D))
kf_dev = (abs(k_f-k_eq)/k_f)*100
print("kf(m²):%.5g\t deviation: %e\n" % (k_f,kf_dev))

print_to_file = open(path_to_folder + project_name + ".txt", "a")
#print_to_file.write("mesh,triangles,dofs,solver,avg_ah,x,y,z,inflow,outflow,deviation,k_eq,k_f,exec_time,Area,ah_lub\n")
print_to_file.write("%s,%g,%10.1g,%s,%e,%e,%e,%e,%.10g,%.10g,%e,%.5g,%.5g,%g,%g,%e\n" % (case_filename,triangles,dofs,solver_method,int_H2,x_,y_,h_,inflow,outflow,deviation,k_eq,k_f,exec_time,Volume,a_H))
print_to_file.close()