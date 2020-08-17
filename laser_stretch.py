import espressomd
import object_in_fluid as oif

from espressomd import lb
from espressomd import lbboundaries
from espressomd import shapes
from espressomd import interactions

import numpy as np
import os, sys, time

# this is a script to try stretch sphere with laser

#####################################
#INPUT DATA
######################################

sim_no = sys.argv[1]

# CELL
radius = 3.13
noNodes = 642
ks = 0.1
kb = 0.0054
kal = 0.02
kag = 0.7
kv = 0.9
kvisc = 0.0
mass = 100.0

#GEOMETRY
boxX = 30.0
boxY = 10.0
boxZ = 20.0

originX = boxX/2
originY = boxY/2
originZ = boxZ/2
partMass = mass/noNodes
stretch = radius

#INPUT FILES
fileNodes = "input/sphere" + str(noNodes) + "nodes.dat"
fileTriangles = "input/sphere" + str(noNodes) + "triangles.dat"

#LBFLUID
grid = 1.0
dens = 1.0
visc = 1.5

#FRICTION
# surface is the refence sphere with r = 4

surface = 4*np.pi*4**2
surface_sph = 4*np.pi*radius**2
friction =(393.0/noNodes)*np.sqrt(surface_sph/surface)*((5.6-1.82)/(6.0/1.025-1.5)*(visc-1.5)+(10-1.82)/(6-1.025)*(dens-1.025)+1.82)

#print ("surface_sph:" + str(surface_sph))
#print ("friction_sph:" + str(friction))

#ITERATION PARAMETERS
maxCycle = 500000
integr_steps = 100
time_step = 0.1
noLoops = 200000
counter = 0
noVtk = 100

accuracy = 0.001
previousL = -1.0
previousA = -1.0

#temp variables
tempVolume = 0.0 #I want to know the maximum deviation of volume and surface
tempSurface = 0.0
tempCycle = 0
cycle = 0

directory = "output/sim"+str(sim_no)
os.makedirs(directory)

#####################################
#initialization
######################################

system = espressomd.System(box_l=[boxX, boxY, boxZ])
system.cell_system.skin = 0.2
system.time_step = time_step

# creating the template for cells
cellType = oif.OifCellType(nodes_file=fileNodes, triangles_file=fileTriangles, system = system, ks=ks, kb=kb, kal=kal, kag=kag, kv=kv, kvisc=kvisc,\
                       resize =[stretch,stretch,stretch], normal=True)

# creating the RBC
cell = oif.OifCell(cell_type=cellType, particle_type=0, origin=[originX,originY,originZ], rotate=[np.pi/2.0, 0.0, 0.0], particle_mass=partMass)

# fluid
lbf = lb.LBFluid(agrid=grid, dens=dens, visc=visc, tau=time_step, fric=friction)
system.actors.add(lbf)

# laser stretching - forces applied to points       
laser = [1.0, 0.0, 0.0]
laser_length = oif.norm(laser)
sigma_0 = 0.00147
        
# for p in cell.mesh.points:
    # neighbors = cell.mesh.neighbors[p.id]
    # normal = neighbors.outer_normal()
    # lengthNormal = oif.norm(normal)
    # posI = p.get_pos()[0]
    # if (posI>originX):
        # laser = [1.0, 0.0, 0.0]
    # if (posI<originX):
        # laser = [-1.0, 0.0, 0.0]
    # angle = np.arccos(np.dot(laser,normal)/(laser_length * lengthNormal))
    # sigma = sigma_0*(np.cos(angle))**2
    # AppliedForce = np.dot(sigma,normal)
    # p.set_force(AppliedForce)

###################################################################
                        # MAIN LOOP
###################################################################

originalVolume = cell.volume()
originalSurface = cell.surface()

cell.output_vtk_pos_folded(file_name="output/sim" + str(sim_no) + "/cell_0.vtk")

stopSimulation = 0

while (counter < noLoops):
    startTime = time.time()
    system.integrator.run(steps=integr_steps)
    #print counter
    counter += 1
    cycle = counter*integr_steps
    # prepocitanie sil natahovacich
    for p in cell.mesh.points:
        neighbors = cell.mesh.neighbors[p.id]
        normal = neighbors.outer_normal()
        lengthNormal = oif.norm(normal)
        posI = p.get_pos()[0]
        if (posI>originX):
            laser = [1.0, 0.0, 0.0]
        if (posI<originX):
            laser = [-1.0, 0.0, 0.0]
        angle = np.arccos(np.dot(laser,normal)/(laser_length * lengthNormal))
        sigma = sigma_0*(np.cos(angle))**2
        AppliedForce = np.dot(sigma,normal)
        p.set_force(AppliedForce)
    # print "time: ", str(counter*time_step*integr_steps)

    # vtk file, normal I write only last vtk for given force and combination of coefficients
    if (counter % noVtk == 0):
        cell.output_vtk_pos_folded(file_name="output/sim" + str(sim_no) + "/cell_" + str(counter) + ".vtk")


    xMin = cell.pos_bounds()[0]
    xMax = cell.pos_bounds()[1]
    length = xMax - xMin

    zMin = cell.pos_bounds()[4]
    zMax = cell.pos_bounds()[5]
    thickness = zMax - zMin
    
    volume = cell.volume()
    surface = cell.surface()
    devVolume = (volume - originalVolume)*100/originalVolume
    devSurface = (surface - originalSurface)*100/originalSurface

    if (counter % noVtk  == 0):
        fout_MS = open("output/sim" + str(sim_no)+"outDatMS.dat", "w")
        fout_MS.write ("f l a l - previous_l previous_a - a devSurface devVolume ""\n")
        if (abs(devSurface) > 3 or abs(devVolume) > 3):
            fout_MS.write("Error: devSurface = "+str(devSurface)+", devVolume =  " +str(devVolume))
            stopSimulation = 1
        else:
            fout_MS.write(str(length)+" "+str(thickness)+" "+str(length-previousL)+" "+str(previousA-thickness)+" "+str(devSurface)+" "+str(devVolume))
        fout_MS.close()
        if ((length-previousL)/previousL < 0.0001 and counter > 100):
            stopSimulation = 1
        previousA = thickness
        previousL = length


    ### save the maximum of volume and surface deviation
    if ( abs(devVolume) > tempVolume):
        tempVolume = abs(devVolume)
        tempCycle = cycle

    if (abs(devSurface) > tempSurface):
        tempSurface = abs(devSurface)
        tempCycle = cycle

    if (stopSimulation == 1):
        counter = noLoops

# last vtk for this stretching
# cell.output_vtk_pos("output/sim" + str(sim_no) + "_force" + str(force) + ".vtk")
elapsedTime=time.time()-startTime
print "Simulation completed in: " +str(elapsedTime)
foutFinalData = open("output/sim" + str(sim_no)+"_finalData.dat","a")
foutFinalData.write(str(length)+" "+str(thickness)+" "+str(volume)+" "+str(surface)+" "+str(cycle)+" "+str(tempVolume)+" "+str(tempSurface)+" "\
                    +str(tempCycle)+" "+str(time))
        
# for i in range(1,maxCycle):
    # system.integrator.run(steps=integr_steps)
    # vtk file for cell
    # cell.output_vtk_pos_folded(file_name="output/sim" + str(sim_no) + "/cell_" + str(i) + ".vtk")
    # prepocitanie sil natahovacich
    # for p in cell.mesh.points:
        # neighbors = cell.mesh.neighbors[p.id]
        # normal = neighbors.outer_normal()
        # lengthNormal = oif.norm(normal)
        # posI = p.get_pos()[0]
        # if (posI>originX):
            # laser = [1.0, 0.0, 0.0]
        # if (posI<originX):
            # laser = [-1.0, 0.0, 0.0]
        # angle = np.arccos(np.dot(laser,normal)/(laser_length * lengthNormal))
        # sigma = sigma_0*(np.cos(angle))**2
        # AppliedForce = np.dot(sigma,normal)
        # p.set_force(AppliedForce)
    # print "time: ", str(i*time_step*integr_steps)

