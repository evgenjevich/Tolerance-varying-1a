
# coding: utf-8
# 1_edit.py file in /data/aem1/new1a/benchmarks
# In[ ]:

import numpy as np
import sympy
import fipy as fp
from fipy import numerix as nmx
import matplotlib.pyplot as plt
import os
import json
import sys


jsonfile = sys.argv[1]


if jsonfile:
    with open(jsonfile, 'rb') as ff:
        params = json.load(ff)
        
else:
    params = dict()
    
print 'my params:', params

#extract the parameters
N = params.get('N', 20)  
total_steps = params.get('steps', 2)
sumatra_label = params.get('sumatra_label', '')
total_sweeps = params.get ('sweeps', 2)
duration = params.get('duration', 500)
problem = params.get('problem', 'a')
tolerance = params.get('tolerance',0.1)


c, rho_s, c_alpha, c_beta = sympy.symbols("c_var rho_s c_alpha c_beta")
f_0 = rho_s * (c - c_alpha)**2 * (c_beta - c)**2

sympy.diff(f_0, c, 2)

# command format:
# python 1_edit.py a 2.0 100
# where nx = 100
# and dx = 2.0
# for benchmark problem 1a
# or
# python 1_edit.py d 100 2500
# for a sphere, radius 100 and containing 2500 cells

if (problem == "a"):
    print "Creating mesh for problem a"
    mesh = fp.PeriodicGrid2D(nx=N, ny=N, dx=200.0/float(N), dy=200.0/float(N))
    
elif (problem == "b"):
    print "Creating mesh for problem b"
    mesh = fp.Grid2D(nx=N, ny=N, dx=200.0/float(N), dy=200.0/float(N))
    
elif (problem == "c"):
    print "Creating mesh for problem c"
#    mesh = fp.Grid2D(dx=0.5, dy=0.5, nx=40, ny=200) + (fp.Grid2D(dx=0.5, dy=0.5, nx=200, ny=40) + [[-40],[100]])
    mesh = fp.Grid2D(Lx=20., Ly=100.0, nx=N / 5, ny=N) + (fp.Grid2D(Ly=20.0, Lx=100.0, nx=N, ny=N / 5) + [[-40],[100]])
       
elif (problem == "d"):
    print "Creating mesh for problem da"
    r = float(sys.argv[2])
    numCellsDesired = int(sys.argv[3])
    epsilon = 0.05
    cellSize = 16 * np.pi * r**2 / (nmx.sqrt(3.) * numCellsDesired)
    cellSize = nmx.sqrt(cellSize)
    
    substring1 = '''
    radius = {0};
    cellSize = {1};
    '''.format(r, round(cellSize, 6))

    mesh = fp.Gmsh2DIn3DSpace(substring1 + '''

    // create inner 1/8 shell
    Point(1) = {0, 0, 0, cellSize};
    Point(2) = {-radius, 0, 0, cellSize};
    Point(3) = {0, radius, 0, cellSize};
    Point(4) = {0, 0, radius, cellSize};
    Circle(1) = {2, 1, 3};
    Circle(2) = {4, 1, 2};
    Circle(3) = {4, 1, 3};
    Line Loop(1) = {1, -3, 2};
    Ruled Surface(1) = {1};

    // create remaining 7/8 inner shells
    t1[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{1};}};
    t2[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{1};}};
    t3[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{1};}};
    t4[] = Rotate {{0,1,0},{0,0,0},-Pi/2} {Duplicata{Surface{1};}};
    t5[] = Rotate {{0,0,1},{0,0,0},Pi/2} {Duplicata{Surface{t4[0]};}};
    t6[] = Rotate {{0,0,1},{0,0,0},Pi} {Duplicata{Surface{t4[0]};}};
    t7[] = Rotate {{0,0,1},{0,0,0},Pi*3/2} {Duplicata{Surface{t4[0]};}};

    // create entire inner and outer shell
    Surface Loop(100)={1, t1[0],t2[0],t3[0],t7[0],t4[0],t5[0],t6[0]};
    ''', order=2.0).extrude(extrudeFunc=lambda r: 1.01*r)

c_alpha = 0.3
c_beta = 0.7
kappa = 2.0
M = 5.0
c_0 = 0.5
epsilon = 0.01
rho_s = 5.0



# solution variable
c_var = fp.CellVariable(mesh=mesh, name=r"$c$", hasOld=True)

# array of sample c-values: used in f versus c plot
vals = np.linspace(-.1, 1.1, 1000)

if (problem == 'a' or 'b' or 'c'):
    # 2D mesh coordinates
    x, y = np.array(mesh.x), np.array(mesh.y)
    # initial value for square and T domains
    c_var[:] = c_0 + epsilon * (np.cos(0.105 * x) * np.cos(0.11 * y) + (np.cos(0.13 * x) * np.cos(0.087 * y))**2 + np.cos(0.025 * x - 0.15 * y) * np.cos(0.07 * x - 0.02 * y))
if (problem == 'd'):
    print "number of cells: " , mesh.numberOfCells
    # 3D mesh coordinates
    x, y, z = np.array(mesh.x), np.array(mesh.y), np.array(mesh.z)
    
    # convert from rectangular to spherical coordinates
    theta = fp.CellVariable(name=r"$\theta$", mesh=mesh)
    theta = nmx.arctan2(z, nmx.sqrt(x**2 + y**2))
    phi = fp.CellVariable(name=r"$\phi$", mesh=mesh)
    phi = nmx.arctan2(y, x)
     
    # initial value for spherical domain
    c_var[:]  = c_0 + epsilon * ((np.cos(8*theta))*(np.cos(15*phi)) + ((np.cos(12*theta))*(np.cos(10*phi)))**2 + ((np.cos(2.5*theta - 1.5*phi))*(np.cos(7*theta - 2*phi))))

# bulk free energy density
def f_0(c):
    return rho_s*((c - c_alpha)**2)*((c_beta-c)**2)
def f_0_var(c_var):
    return 2*rho_s*((c_alpha - c_var)**2 + 4*(c_alpha - c_var)*(c_beta - c_var) + (c_beta - c_var)**2)
# free energy
def f(c):
    return (f_0(c)+ .5*kappa*(c.grad.mag)**2)

#Method for making sure we save .mpz.npz at specified dump_times
def calc_dt(elapsed_time, dt, dt_old, dump_to_file, dump_times, filename):
    if dump_to_file: #if this is true, we have alreay saved the necessary .mpz.npz file
        dt = dt_old #reset back to normal dt
        dt = dt * 1.1 #continue as normal
        dump_to_file = False
        # filename = '1a_{0}_step{1}_data_time-{2:.2f}.npz'.format(N, str(steps).rjust(6, '0'), elapsed+dt)
    else:
        dt_old = dt
        dt = dt * 1.1
        # import ipdb; ipdb.set_trace()

        if len(dump_times) > 0:
            if (elapsed_time + dt) >= dump_times[0]:
                dt = dump_times[0] - elapsed_time
                dump_to_file = True
                
                #dump_time files will have .mpz.npz extension
                # filename = '1a_{0}_step{1}_data_time-{2:.2f}.mpz.npz'.format(N, str(steps).rjust(6, '0'), elapsed+dt)
                del dump_times[0]

    return dt, dt_old, dump_times, dump_to_file, filename
    
# save elapsed time and free energy at each data point
time_data = []
cvar_data = []
f_data = []
# checks whether a folder for the pickles from this simulation exists
# if not, creates one in the home directory
filepath = os.path.join('/data/aem1/new1a/benchmarks/Data', sumatra_label)


# solver equation    
eqn = fp.TransientTerm(coeff=1.) == fp.DiffusionTerm(M * f_0_var(c_var)) - fp.DiffusionTerm((M, kappa))

elapsed = 0.0
steps = 0
dt = 0.01
#tolerance = 1e-1
dump_times = [1.0, 5.0, 10.0, 20.0, 100.0, 200, 500, 1000, 2000, 3000, 10000]
dt_old = dt
dump_to_file = False
filename = 'anushka'


c_var.updateOld()
from fipy.solvers.pysparse import LinearPCGSolver as Solver
solver = Solver()
print "Starting Solver."
while (steps <= total_steps) and elapsed <= duration:
    res0 = eqn.sweep(c_var, dt=dt, solver=solver)
    #record the volume integral of the free energy 
    # equivalent to the average value of the free energy for any cell,
    # multiplied by the number of cells and the area of each cell
    # (since this is a 2D domain)
    
    for sweeps in range(total_sweeps):
        res = eqn.sweep(c_var, dt=dt, solver=solver)
                    
    if res < res0 * tolerance:  
        # anything in this loop will only be executed every 10 steps
        if dump_to_file or (steps%10==0):
            print steps
            print elapsed 
            print "Saving data"
            if (problem == 'a') or (problem == 'b') or (problem == 'c'):
                ##print "saving for a or b or c!!"
                ##save_data(elapsed, c_var, f(c_var).cellVolumeAverage*mesh.numberOfCells*dx*dx)
                if dump_to_file: filename = '1{3}_{0}_step{1}_data_time-{2:.2f}.mpz.npz'.format(N, str(steps).rjust(6, '0'), elapsed, problem)
                else: filename = '1{3}_{0}_step{1}_data_time-{2:.2f}.npz'.format(N, str(steps).rjust(6, '0'), elapsed. problem)
                #Saving necessary Data
                np.savez(os.path.join(filepath,filename),
                     c_var_array=np.array(c_var),
                     dt=dt,
                     elapsed=elapsed,
                     steps=steps,
                     dx=c_var.mesh.dx,
                     dy=c_var.mesh.dy,
                     nx=c_var.mesh.nx,
                     ny=c_var.mesh.ny,
                     sweeps = total_sweeps)
                
            elif (problem == 'd'):
                print "saving for d!!"
                if dump_to_file: filename = '1d_{0}_step{1}_data_time-{2:.2f}.mpz.npz'.format(N, str(steps).rjust(6, '0'), elapsed)
                else: filename = '1d_{0}_step{1}_data_time-{2:.2f}.npz'.format(N, str(steps).rjust(6, '0'), elapsed)
                save_data(elapsed, c_var, )

                np.savez(os.path.join(filepath,filename),
                     c_var_array=np.array(c_var),
                     dt=dt,
                     elapsed=elapsed,
                     steps=steps,
                     dx=c_var.mesh.dx,
                     dy=c_var.mesh.dy,
                     nx=c_var.mesh.nx,
                     ny=c_var.mesh.ny,
                     sweeps = total_sweeps)

        dt, dt_old, dump_times, dump_to_file, filename = calc_dt(elapsed, dt, dt_old, dump_to_file, dump_times, filename)
        steps += 1
        elapsed += dt
        # dt *= 1.1 #this should probably not be here since it is inculded in the calc_dt() method
        c_var.updateOld()
    else:
        dt *= 0.8
        c_var[:] = c_var.old

# simulation ends
print 'steps reached:', steps


