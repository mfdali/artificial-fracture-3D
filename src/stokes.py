from dolfin import *
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import time

#############
## Path
path_to_folder = "/home/monique/fem-fenics/3D/"
#case_filename = "slope_1mm_meshsize_0005" 0.001 1.5006E-03
case_filename = "x2_reduced_1mm_meshsize_0003"
project_name = "test_3d_meshes"
solver_method = "mumps"

#Geometry
x_ = 0.05 # Lenght
y_ = 0.025 # Width
h_ = 0.0009
tol = 1E-6
h_avg = 0.001

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

dx = Measure("dx",domain=mesh, subdomain_data=mf) #or domains
print(assemble(Constant(1)*dx))

# Define boundary condition
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left = AutoSubDomain(lambda x: near(x[0], 0.0))
right = AutoSubDomain(lambda x: near(x[0], x_) )
bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
top = AutoSubDomain(lambda x: near(x[1], y_))
back = AutoSubDomain(lambda x: near(x[2], 0.0))
# Smooth surface
front = AutoSubDomain(lambda x: x[2] > h_-tol)
# Rigid surface
#front = AutoSubDomain(lambda x: near(x[2], h_))


# Define boundary markers
left.mark(boundaries, 1)
top.mark(boundaries, 2)
right.mark(boundaries, 3)
bottom.mark(boundaries, 4)
back.mark(boundaries, 5)
front.mark(boundaries, 6)

ds = Measure('ds',domain=mesh, subdomain_data = boundaries)

################################
#Material properties

mu = Constant(0.001) #Water viscosity [Pa.s]
mu_ = 7.5*mu #R.C. Givler, S.A. Altobelli, J. Fluid Mech. 258 (1994) 355.
pin = Constant(2*6894.76) #Imposed pressure at the entrance [Pa]
pout = Constant(6894.76) #Imposed pressure on output [Pa]

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

##############################################################

#Define Dirichlet boundary conditions
bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,1)
bc2 = DirichletBC(W.sub(0),noslip,boundaries,2)
bc3 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,3)
bc4 = DirichletBC(W.sub(0),noslip,boundaries,4)
bc5 = DirichletBC(W.sub(0),noslip,boundaries,5)
bc6 = DirichletBC(W.sub(0),noslip,boundaries,6)

bcs = [bc1,bc2,bc3,bc4,bc5,bc6]
'''
#Define measures associated with the boundaries and holes
ds = Measure('ds',domain=mesh, subdomain_data = boundaries)
dx = Measure('dx',domain=mesh)
'''

##############################################################
# Compensate the drag at glass walls (3D)
#Drag_force = -(12/(h_**2))*mu*inner(u,v)*dx

#Define variational form for Stokes
#a = (mu_*inner(grad(u),grad(v))*dx(1) + (mu/k)*inner(u,v)*dx(0) - div(v)*p*dx(1) - div(v)*p*dx(0)-div(u)*q*dx(1) - div(u)*q*dx(0))
#L = (inner(f,v)*dx(1) + inner(f,v)*dx(0) - pin*dot(v,n)*ds(1) - pout*dot(v,n)*ds(3))

F = (mu*inner(grad(u),grad(v))*dx - div(v)*p*dx - div(u)*q*dx) + pin*dot(v,n)*ds(1) + pout*dot(v,n)*ds(3) - inner(f,v)*dx #- Drag_force
(a, L) = system(F)

#Compute solution
U = Function(W, name = "field")
#solve(a,U.vector(),L)
'''
problem = LinearVariationalProblem(a, L, U, bcs)
solver = LinearVariationalSolver(problem)
solver.solve()
'''
A  = assemble(a)
b = assemble(L)
[bc.apply(A,b) for bc in bcs]
#[bc.apply(U.vector()) for bc in bcs]
solve(A,U.vector(),b,solver_method)

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

'''print(exec_time)
#Save data to file
#data_file.write("%s,%s,%s,%g,%g,%g,%g,%2.10g,%2.10g,%g\n" % (mode,mesh_type,sm,triangles,exec_time,k,ps,inflow,outflow,k_medium))
'''

# Save solution in VTK format
file = File(path_to_folder + case_filename + "_" + str(h_) + "-u.pvd")
file << u
file = File(path_to_folder + case_filename + "_" + str(h_) + "-p.pvd")
file << p
'''
# Plot field and save figure
plt.figure()
plot(u, mode = "glyphs")'''
'''
plot(u,
    wireframe = True,              # use wireframe rendering
    scalarbar = False,             # hide the color mapping bar
    scale = 2.0,                    # scale the warping/glyphs
    title = "Fancy plot"           # Set your own title
    )'''
'''
plt.xlabel('u')
plt.savefig(path_to_folder + case_filename + "_" + str(h_) + "-u.png")
plt.close()
'''
# Measure height
# Define the Z-coordinate function
z = Expression('x[2]', degree=1)

# Calculate the total volume of the mesh
Volume = float(assemble(1 * ds(6)))  # Total volume of the mesh

# Calculate the integral of the Z-coordinate over the mesh
z_integral = assemble(z * ds(6))

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

k_f = (int_H2*2)/12
print("keq(m²):%.5g\t keq(D):%g\n" % (k_eq,k_D))
print("kf(m²):%.5g\n" % (k_f))

print_to_file = open(path_to_folder + project_name + ".txt", "a")
#print_to_file.write("mesh,triangles,dofs,solver,avg_ah,x,y,z,inflow,outflow,deviation,k_eq,k_f,exec_time,Area,ah_lub\n")
print_to_file.write("%s,%g,%10.1g,%s,%e,%e,%e,%e,%.10g,%.10g,%e,%.5g,%.5g,%g,%g,%e\n" % (case_filename,triangles,dofs,solver_method,int_H2,x_,y_,h_,inflow,outflow,deviation,k_eq,k_f,exec_time,Volume,a_H))
print_to_file.close()