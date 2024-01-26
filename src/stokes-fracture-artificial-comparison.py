import csv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import distance
from dolfin import *
from mshr import *
import random #pseudorandom
import numpy as np
#import matplotlib.pyplot as pl
from math import pi, sqrt, floor, ceil,pow
#import time
from datetime import datetime
#from secrets import randbelow #not cryptographically secure
import time
import matplotlib.pyplot as plt
from ufl.operators import min_value
from scipy.interpolate import LinearNDInterpolator,NearestNDInterpolator

class Params:
  def __init__(self):
    self.muw = 0.001
    self.Swc = 0
    self.nw = 2
    self.muo = 0.001
    self.rhoo = 1000
    self.rhow = 1000
    self.kromax = 1
    self.krwmax = 1
    
    self.Sor = 0
    self.no = 2
    self.Cw = 0
    self.Co = 0
    self.kpm = 1E-11
    self.pin = 80#5*6894.76
    self.pout = 10#6894.76

class MeshHeight(UserExpression):

    def __init__(self, **kwargs):
        #random.seed(2 + MPI.rank(MPI.comm_world))
        #random.seed(datetime.now())
        self.null_value = 1E-15#int(0)
        self.tol = 1E-10
        super().__init__(**kwargs)
        
    def set_h_values(self,X,Y,Z):
        self.X, self.Y, self.Z = X, Y, Z
        self.interpolator = LinearNDInterpolator((self.X, self.Y), self.Z)
        #self.interpolator = NearestNDInterpolator((self.X, self.Y), self.Z)

    def eval(self, values, x):

        interpolator_ = self.interpolator([x[0], x[1]])
        if (interpolator_[0] <= self.tol):
            values[0] = self.null_value #np.array([0.0])
        else:
            values[0] = interpolator_

    def value_shape(self):
        return ()

#class LoadAperture():

def load_data(filename):
    x = []
    y = []
    z = []

    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        #next(csvreader)  # Skip the header row if present
        for row in csvreader:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]))

    return np.array(x), np.array(y), np.array(z)

def _calculate_z(x, y, z, x_query, y_query):
    points = np.column_stack((x, y))
    values = z
    z_query = griddata(points, values, (x_query, y_query), method='nearest')
    return z_query

def matrixZ(X,Y,x,y,z):
    Z = _calculate_z(x, y, z, X, Y)
    return Z

def calculate_distance(Zup, Zdown):
    distance = np.sqrt((Zup - Zdown) ** 2)
    return distance

def get_distance(x, y):
  min_distances = []
  i = 20
  point1 = (x[i], y[i])  # Coordinates of the current point
  min_distance = float("inf")  # Initialize with a large value
  for j in range(len(x)):
      if i != j:  # Skip the current point
          point2 = (x[j], y[i])  # Coordinates of another point # Keep x change y
          dist = distance.euclidean(point1, point2)
          if dist > 0:
            min_distance = min(min_distance, dist)
      if min_distance > 0 or min_distance < float('inf'):
        min_distances.append(min_distance)
  
  dist = np.unique(min_distances)
  scientific_notation="{:.1e}".format(dist[0])
  x_dist = float(scientific_notation)

  return x_dist

#class NumModel():

def Dolfin_Mesh(resx,resy,a,b):

    # Domain size
    larea = a
    harea = b
    xsteps = int(round(a/resx))+1
    ysteps = int(round(b/resy))+1
    print(xsteps,ysteps)
    '''# Domain construction
    rectangle = Rectangle(Point(0., 0.), Point(larea, harea))
    domain = rectangle

    # Mesh generation
    mesh = generate_mesh(domain,res)'''
    mesh = RectangleMesh(Point(0., 0.), Point(larea, harea),xsteps,ysteps)

    # Define boundary condition
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    left = AutoSubDomain(lambda x: near(x[0], 0.0))
    right = AutoSubDomain(lambda x: near(x[0], larea))
    bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
    top = AutoSubDomain(lambda x: near(x[1], harea))

    # Define boundary markers
    left.mark(boundaries, 1)
    top.mark(boundaries, 2)
    right.mark(boundaries, 3)
    bottom.mark(boundaries, 4)
    '''
    # 2D plot
    plot(mesh)
    #plt.axis('off')
    plt.show()
    plt.savefig("mesh-savefig.png")'''
    return mesh, boundaries

def Max(a, b): return (a+b+abs(a-b))/Constant(2)
def Min(a, b): return (a+b-abs(a-b))/Constant(2)


def thickness(xah, yah, zah,x_,y_,case_filename):

    from matplotlib import cm
    from matplotlib.ticker import LinearLocator

    eq = case_filename

    #grid_pos = grid_position()
    
    #Aperture
    h = MeshHeight()
    h.set_h_values(xah,yah,zah)
    
    tol = 0.001 # avoid hitting points outside the domain
    x = np.linspace(0, x_, 1001)
    y = np.linspace(0, y_, 1001)
    X,Y = np.meshgrid(x,y)

    
    Z = matrixZ(X,Y, xah, yah,zah)

    fig = plt.figure(figsize=(9, 4))

    # Plot the surface.
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis, rstride=4, cstride=4, 
                        linewidth=0.2)#, antialiased=False)
    ax.view_init(30, 70)
    fig.colorbar(surf, shrink=0.5, aspect=10, location='bottom',orientation='horizontal')
    
    ax = fig.add_subplot(1, 2, 2)
    CS = ax.contour(X, Y, Z)
    ax.clabel(CS, inline=True, fontsize=10)
    ax.grid(True)
    fig.savefig(output_folder + case_filename + '-all.png',dpi=120)

    plt.close()
    
    return h,eq


def Num_Solution(x_,y_,xah, yah, zah,x_dist,y_dist,case_filename,output_folder):
    
    ### Parameters ###
    resx = x_dist
    resy = y_dist
    start_time = time.time()

    # Load mesh
    mesh, boundaries = Dolfin_Mesh(resx,resy,x_,y_)

    #Material properties
    Prop = Params()
    mu = Constant(Prop.muw) #Water viscosity [Pa.s]
    mu_ = 7.5*mu #R.C. Givler, S.A. Altobelli, J. Fluid Mech. 258 (1994) 355.
    pin = Constant(Prop.pin) #Imposed pressure at the entrance [Pa]
    pout = Constant(Prop.pout) #Imposed pressure on output [Pa]
    thick,eq = thickness(xah, yah, zah, x_, y_, case_filename) #Varing the height

    ####### Numerical Solution #######

    #Function space over the mesh
    V = VectorElement('CG',mesh.ufl_cell(),2)
    Q = FiniteElement('CG',mesh.ufl_cell(),1)
    Element = V*Q
    W = FunctionSpace(mesh,Element)

    Y = FunctionSpace(mesh, "Lagrange", 1)
    h_test = interpolate(thick,Y)
    h_min = h_test.vector().min()
    h_max = h_test.vector().max()
    h_avg = np.average(h_test.vector())
    h_vec = h_test.vector()
    h_avg2 = h_vec.sum()/h_vec.size()

    #Define variational problem
    (u,p) = TrialFunctions(W)
    (v,q) = TestFunctions(W)
                    
    ##############################################################
    # Parameters

    # Define expressions used in variational forms
    dp = pin-pout
    #g = Expression("b-(a/l)*y[0]" , degree=1 ,a = Constant(dp), b = Constant(pin), l = 0.027) #g = div(u)
    #u_in = Constant((0.0)) #Initial velocity in x [m/s]
    noslip = Constant((0.0,0.0)) #No-slip condition for velocity, u=0 at y=h
    f = Constant((0.0,0.0)) #External force
    n = FacetNormal(mesh) #Normal vector to mesh

    ##############################################################

    #Define Dirichlet boundary conditions
    bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,1)
    bc2 = DirichletBC(W.sub(0),noslip,boundaries,2)
    bc3 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,3)
    bc4 = DirichletBC(W.sub(0),noslip,boundaries,4)

    bcs = [bc1,bc2,bc3,bc4]

    '''mesh coordinates'''
    D = mesh.topology().dim()
        
    # get vertex coordinates
    coords = mesh.coordinates()
    
    # get elements at BC1
    mks = bc3.markers()

    #Define measures associated with the boundaries and holes
    ds = Measure('ds',domain=mesh, subdomain_data = boundaries)
    dx = Measure('dx',domain=mesh)

    ##############################################################
    # Compensate the drag at glass walls (3D)
    Drag_force = -(12/(thick**2))*mu*inner(u,v)*dx

    #Define variational form for Stokes
    #a = (mu_*inner(grad(u),grad(v))*dx(1) + (mu/k)*inner(u,v)*dx(0) - div(v)*p*dx(1) - div(v)*p*dx(0)-div(u)*q*dx(1) - div(u)*q*dx(0))
    #L = (inner(f,v)*dx(1) + inner(f,v)*dx(0) - pin*dot(v,n)*ds(1) - pout*dot(v,n)*ds(3))

    F = (mu*inner(grad(u),grad(v))*dx - div(v)*p*dx - div(u)*q*dx) + pin*dot(v,n)*ds(1) + pout*dot(v,n)*ds(3) - Drag_force - inner(f,v)*dx
    (a, L) = system(F)

    #Compute solution
    U = Function(W, name = "field")
    A  = assemble(a)
    b = assemble(L)
    [bc.apply(A,b) for bc in bcs]
    #[bc.apply(U.vector()) for bc in bcs]
    solve(A,U.vector(),b,"mumps")

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

    
    # Save solution in VTK format
    file = File(output_folder + case_filename + "-u.pvd")
    file << u
    file = File(output_folder + case_filename + "-p.pvd")
    file << p
    file = File(output_folder + case_filename + "-ah.pvd")
    file << h_test

    # Plot field and save figure
    plt.figure()
    plot(u, mode = "glyphs")
    
    plt.xlabel('u')
    plt.savefig(output_folder + case_filename + "-u.png")
    plt.close()

    print("\nMin: %.5e"%h_min)
    print("\nMax: %.5e"%h_max)
    print("\nAverage: %.5e"%h_avg)
    print("\nAverage: %.5e"%h_avg2)

    # Calculate aperture average using mesh
    a_H2 = assemble(thick*dx) 
    Area = float(assemble(1 * dx))
    int_H2 = a_H2/Area
    print("a_h_average2 = %e" % int_H2)

    try:
        a_H = assemble((1/(thick**3))*dx)
        Area = float(assemble(1 * dx))
        int_H = pow((Area/a_H),(1/3))
        print("int dA/hy^3 = %e" % a_H)
        print("a_h_average = %e" % int_H)

    except:
        a_H = 0
        int_H = 0

    #mesh size
    triangles = mesh.num_cells()
    print("%g\n" % triangles)
    info("Num DOFs {}".format(W.dim()))
    dofs = W.dim()

    k_eq = (outflow * mu * x_)/(y_*dp)
    k_D = k_eq*1013250000000

    k_f = pow(int_H2,2)/12
    print("keq(m²):%.5g\t keq(D):%g\n" % (k_eq,k_D))
    kf_dev = (abs(k_f-k_eq)/k_f)*100
    print("kf(m²):%.5g\t deviation: %.5g\n" % (k_f,kf_dev))
    
    print_to_file.write("%s,%s,%g,%e,%e,%.10g,%.10g,%e,%g,%g,%g,%g,%g,%g,%g,%g\n" % (case_filename,eq,int_H2,h_min,h_max,inflow,outflow,deviation,k_eq,exec_time,Area,a_H,int_H,triangles,dofs,exec_time))

def _geometry_params(aperture):
    x, y, z = load_data(aperture) 
    a1 = np.min(x)
    b1 = np.min(y)
    a2 = np.max(x)
    b2 = np.max(y)

    x_ = a2-a1 # Lenght
    y_ = b2-b1 # Width

    xah = x - a1
    yah = y - b1
    zh_ = z
    
    return x_,y_,xah, yah, zh_

'''Artificial'''
def slope(output_folder):
    
    aperture = output_folder + 'slope10N_coords.csv'
    project_name = '3d-artificial-fracture'
    case = ['2.5D-slope-10N']
    h_ = [0]

    x, y, z = load_data(aperture) 
    a1 = np.min(x)
    b1 = np.min(y)
    a2 = np.max(x)
    b2 = np.max(y)

    x_ = a2-a1 # Lenght
    y_ = b2-b1 # Width

    xah = x - a1
    yah = y - b1
    zh_ = z
    
    return x_,y_,xah, yah, zh_, aperture, project_name, case, h_

def eggshell(output_folder):
    
    aperture = 'eggshell10N_coords.csv'
    project_name = '3d-artificial-fracture'
    case = ['2.5D-eggshell-10N']
    h_ = [0]

    x_,y_,xah, yah, zh_ =_geometry_params(aperture)
    
    return x_,y_,xah, yah, zh_, aperture, project_name, case, h_

def box(output_folder):
    
    aperture = output_folder + 'thinbox_1mm_0002_coords.csv'
    project_name = '3d-artificial-fracture'
    case = ['2.5D-box-1mm-0003']
    h_ = [0]
    
    return x_,y_,xah, yah, zh_, aperture, project_name, case, h_

output_folder = 'fracture_2D/'

x_,y_,xah, yah, zh_, aperture, project_name, case, h_ = eggshell(output_folder)
#x_,y_,xah, yah, zh_, aperture, project_name, case, h_ = slope(output_folder)
#x_,y_,xah, yah, zh_, aperture, project_name, case, h_ = box(output_folder)


x_dist = get_distance(xah, yah)
y_dist = get_distance(yah, xah)
print(x_dist,y_dist)

print_to_file = open(output_folder + project_name + ".txt", "a")
print_to_file.write("Case,Equation,Av_Height,np.min,np.max,Inflow,Outflow,%Deviation,K_eq,Exec_time,Area,Height_3,Av_lub,cells,dof,time\n")

for step in range(len(h_)):
    case_filename = project_name + case[step]
    zah = zh_-h_[step]
    zah[zah<(1e-10)] = int(0)
    Num_Solution(x_,y_,xah,yah,zah,x_dist,y_dist,case_filename,output_folder)

print_to_file.close()
