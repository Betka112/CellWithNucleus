import espressomd
import object_in_fluid as oif

import Boundaries_15
from Boundaries_15 import *

from espressomd import lb
from espressomd import lbboundaries
from espressomd import shapes
from espressomd import interactions
from espressomd.interactions import HarmonicBond

import numpy as np
import os, sys
import random

# this is a script to try the nucleus with bonded interactions

#####################################
#INPUT DATA
######################################
large_number = 10000000.0

# cell
membrane_radius = float(sys.argv[1])
membrane_noNodes = int(sys.argv[2])
membrane_ks = float(sys.argv[3])
membrane_kb = float(sys.argv[4])
membrane_kal = float(sys.argv[5])
membrane_kag = float(sys.argv[6])
membrane_kv = float(sys.argv[7])
membrane_kvisc = float(sys.argv[8])
membrane_mass = float(sys.argv[9])
# nucleus
nucleus_radius = float(sys.argv[10])
nucleus_noNodes = int(sys.argv[11])
nucleus_ks = float(sys.argv[12])
nucleus_kb = float(sys.argv[13])
nucleus_kal = float(sys.argv[14])
nucleus_kag = float(sys.argv[15])
nucleus_kv = float(sys.argv[16])
nucleus_kvisc = float(sys.argv[17])
nucleus_mass = float(sys.argv[18])
# fluid force
fluid_force = float(sys.argv[19])
# M_force is the force applied to each point of the cell (M as membrane)
M_force = float(sys.argv[20])
sim_no = sys.argv[21]
# harmonic bond parameter
par_k = float(sys.argv[22])
angle_par = float(sys.argv[23])
# parameters of the constriction
par_ax = float(sys.argv[24])
par_ay = float(sys.argv[25])
par_az = float(sys.argv[26])

#potrebujes mat pricinok output a v nom priecinok data a vtk, do data sa ti ukladaju datove subory, do vtk sa ti ukladaju vtkacka
directory_data = "output/data/sim"+str(sim_no)
os.makedirs(directory_data)
directory_vtk = "output/vtk/sim"+str(sim_no)
os.makedirs(directory_vtk)
vtk_directory = directory_vtk+"/vtk"
os.makedirs(vtk_directory)

#CELL
originX = 3*membrane_radius/2
originY = boxY/2
originZ = boxZ/2
membrane_partMass = membrane_mass/membrane_noNodes
membrane_stretch = float(membrane_radius)
#NUCLEUS
nucleus_partMass = nucleus_mass/nucleus_noNodes
nucleus_stretch = float(nucleus_radius)

#INPUT FILES
membraneNodes = "input/sphere" + str(membrane_noNodes) + "nodes.dat"
membraneTriangles = "input/sphere" + str(membrane_noNodes) + "triangles.dat"
nucleusNodes = "input/sphere" + str(nucleus_noNodes) + "nodes.dat"
nucleusTriangles = "input/sphere" + str(nucleus_noNodes) + "triangles.dat"

#LBFLUID
LBgrid = 1.0
density = 1.025
viscosity = 1.5
time_step = 0.1

#FRICTION
# surface is the refence sphere with r = 4
surface = 4*np.pi*4**2
surface_sph = 4*np.pi*membrane_radius**2
friction =(393.0/membrane_noNodes)*np.sqrt(surface_sph/surface)*((5.6-1.82)/(6.0/1.025-1.5)*(viscosity-1.5)+(10-1.82)/(6-1.025)*(density-1.025)+1.82)
#print ("surface_sph:" + str(surface_sph))
#print ("friction_sph:" + str(friction))

system = espressomd.System(box_l=[boxX, boxY, boxZ])
system.cell_system.skin = 0.2
system.time_step = time_step

#ITERATION PARAMETERS
cycle = 0
maxCycle = 500000000000
integr_steps = 100

#####################################
#initialization
######################################

# creating templates
membrane_type = oif.OifCellType(nodes_file=membraneNodes, triangles_file=membraneTriangles, check_orientation=False,
                           system=system, ks=membrane_ks, kb=membrane_kb, kal=membrane_kal, kag=membrane_kag,
                           kv=membrane_kv, kvisc=membrane_kvisc, resize=[membrane_stretch, membrane_stretch, membrane_stretch], normal=True)
nucleus_type = oif.OifCellType(nodes_file=nucleusNodes, triangles_file=nucleusTriangles, check_orientation=False,
                          system=system, ks=nucleus_ks, kb=nucleus_kb, kal=nucleus_kal, kag=nucleus_kag,
                          kv=nucleus_kv, kvisc=nucleus_kvisc, resize=[nucleus_stretch, nucleus_stretch, nucleus_stretch], normal=True)

# creating the cell
membrane = oif.OifCell(cell_type = membrane_type, particle_type=0, origin=[originX,originY,originZ], particle_mass=membrane_partMass)
nucleus = oif.OifCell(cell_type = nucleus_type, particle_type=1, origin=[originX,originY,originZ], particle_mass=nucleus_partMass)

#info about the cell membrane and nucleus membrane
surface_sph = membrane.mesh.surface()
min_hrana_sph = membrane.mesh.min_edge_length()
max_hrana_sph = membrane.mesh.max_edge_length()
aver_hrana_sph = membrane.mesh.aver_edge_length()
out_file = open(directory_data + "/sphere_info_" + str(sim_no) + ".dat", "a")
out_file.write ("surface min_edge max_edge ave_edge \n")
out_file.write(str(surface_sph)+" "+str(min_hrana_sph)+" "+str(max_hrana_sph)+" "+str(aver_hrana_sph))
out_file.close()

surface_nuc = nucleus.mesh.surface()
min_hrana_nuc = nucleus.mesh.min_edge_length()
max_hrana_nuc = nucleus.mesh.max_edge_length()
aver_hrana_nuc = nucleus.mesh.aver_edge_length()
out_file = open(directory_data + "/nucleus_info_" + str(sim_no) + ".dat", "a")
out_file.write ("surface min_edge max_edge ave_edge \n")
out_file.write(str(surface_nuc)+" "+str(min_hrana_nuc)+" "+str(max_hrana_nuc)+" "+str(aver_hrana_nuc))
out_file.close()


small_angle = np.pi/angle_par
countbonds = 0
average_dlzka = 0
max_dlzka = - large_number
min_dlzka = large_number

lines = []
pairs = []
for p in nucleus.mesh.points:
    neighbors = nucleus.mesh.neighbors[p.id]
    normal = neighbors.outer_normal()
    lengthNormal = oif.norm(normal)
    pPos = p.get_pos()
    for pMembrane in membrane.mesh.points:
        pMembranePos = pMembrane.get_pos()
        tmpConnection = pMembranePos - pPos
        tmpDist = oif.norm(tmpConnection)
        angle = np.arccos(np.dot(tmpConnection,normal)/(tmpDist * lengthNormal))
        if angle < small_angle:
            harmonicInter = HarmonicBond(k=par_k, r_0=tmpDist)
            system.bonded_inter.add(harmonicInter)
            p.part.add_bond((harmonicInter, pMembrane.part_id))
            tmplist = [pPos[0], pPos[1], pPos[2], pMembranePos[0], pMembranePos[1], pMembranePos[2]]
            tmppair = [p.id, pMembrane.id]
            lines.append(tmplist)
            pairs.append(tmppair)
            countbonds += 1

#initial info about the bonds , number, min_length, max_length, average
for pair in pairs:
    pPos = nucleus.mesh.points[pair[0]].get_pos()
    pMembranePos = membrane.mesh.points[pair[1]].get_pos()
    tmpConnection = pMembranePos - pPos
    dlzka = oif.norm(tmpConnection)
    average_dlzka = average_dlzka + dlzka
    if dlzka > max_dlzka:
        max_dlzka = dlzka
    if dlzka < min_dlzka:
        min_dlzka = dlzka
average_dlzka = average_dlzka / countbonds
out_file = open(directory_data + "/bonds_initial_info_" + str(sim_no) + ".dat", "a")
out_file.write ("NoBonds min_edge max_edge ave_edge \n")
out_file.write(str(countbonds)+" "+str(min_dlzka)+" "+str(max_dlzka)+" "+str(average_dlzka))
out_file.close()
average_dlzka = 0
max_dlzka = - large_number
min_dlzka = large_number


# cell-wall interactions
system.non_bonded_inter[0,10].soft_sphere.set_params(a = 0.0001, n = 1.2, cutoff = 0.1, offset = 0.0)
system.non_bonded_inter[1,10].soft_sphere.set_params(a = 0.0001, n = 1.2, cutoff = 0.1, offset = 0.0)


# fluid
#lbf = espressomd.lb.LBFluid(agrid = LBgrid, dens = density, visc = viscosity, tau = time_step, fric = friction , ext_force_density = [fluid_force, 0.0, 0.0])
#system.actors.add(lbf)

# fluid
lb_params = {'agrid': LBgrid, 'dens': density, 'visc': viscosity, 'tau': time_step,
             'ext_force_density': [fluid_force, 0, 0]}
lbf = espressomd.lb.LBFluid(**lb_params)
system.actors.add(lbf)

system.thermostat.set_lb(LB_fluid=lbf, gamma=friction)

boundaries = list()
FillBoundaries(boundaries, vtk_directory, par_ax, par_ay, par_az)
for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape=boundary))
    system.constraints.add(shape=boundary, particle_type=10, penetrable=False)

#force on individual mesh points
for it in membrane.mesh.points:
    it.set_force((M_force,0.0,0.0))

out_file_sph = open(directory_data + "/sphere_" + str(sim_no) + ".dat", "a")
out_file_sph.write ("cycle time velX velY velZ posX posY posZ sizeX sizeY sizeZ defI")
out_file_sph.write("\n")
out_file_sph.close()

out_file_nuc = open(directory_data + "/nucleus_" + str(sim_no) + ".dat", "a")
out_file_nuc.write ("cycle time velX velY velZ posX posY posZ sizeX sizeY sizeZ defI")
out_file_nuc.write("\n")
out_file_nuc.close()

out_file_diff = open(directory_data + "/difference_" + str(sim_no) + ".dat", "a")
out_file_diff.write ("posX diffMaxX defI")
out_file_diff.write("\n")
out_file_diff.close()

out_file = open(directory_data + "/bonds_info_" + str(sim_no) + ".dat", "a")
out_file.write ("cycle min_edge max_edge ave_edge \n")
out_file.close()
    
###################################################################
                       #MAIN LOOP
###################################################################

membrane.output_vtk_pos_folded(file_name=directory_vtk + "/membrane_0.vtk")
nucleus.output_vtk_pos_folded(file_name=directory_vtk + "/nucleus_0.vtk")
oif.output_vtk_lines(lines=lines, out_file=directory_vtk + "/lines_0.vtk" )


#for i in range(1,maxCycle):
while cycle < maxCycle:
    system.integrator.run(steps=integr_steps)
    # vtk file for membrane
    membrane.output_vtk_pos_folded(file_name=directory_vtk + "/membrane_" + str(cycle//integr_steps) + ".vtk")
    # dat file for membrane
    out_file_sph = open(directory_data + "/sphere_" + str(sim_no) + ".dat", "a")
    sph_velocity = membrane.get_velocity()
    sph_position = membrane.get_origin()
    sph_pos_bounds = membrane.pos_bounds()
    sph_x_size = sph_pos_bounds[1] - sph_pos_bounds[0]
    sph_y_size = sph_pos_bounds[3] - sph_pos_bounds[2]
    sph_z_size = sph_pos_bounds[5] - sph_pos_bounds[4]
    # out_file_sph.write ("cycle time velX velY velZ posX posY posZ sizeX sizeY sizeZ")
    out_file_sph.write(
        str(cycle) + " " + str(cycle*time_step) + " " + str(sph_velocity) + " " + str(sph_position) + " " + str(
            sph_x_size) + " " + str(sph_y_size) + " " + str(sph_z_size) + " " + str(sph_x_size/sph_y_size))
    out_file_sph.write("\n")
    out_file_sph.close()
    # vtk file for nucleus
    nucleus.output_vtk_pos_folded(file_name=directory_vtk + "/nucleus_"  + str(cycle//integr_steps) + ".vtk")
    # dat file for nucleus
    out_file_nuc = open(directory_data + "/nucleus_" + str(sim_no) + ".dat", "a")
    nuc_velocity = nucleus.get_velocity()
    nuc_position = nucleus.get_origin()
    nuc_pos_bounds = nucleus.pos_bounds()
    nuc_x_size = nuc_pos_bounds[1] - nuc_pos_bounds[0]
    nuc_y_size = nuc_pos_bounds[3] - nuc_pos_bounds[2]
    nuc_z_size = nuc_pos_bounds[5] - nuc_pos_bounds[4]
    out_file_nuc.write(
        str(cycle) + " " + str(cycle*time_step) + " " + str(nuc_velocity) + " " + str(nuc_position) + " " + str(
            nuc_x_size) + " " + str(nuc_y_size) + " " + str(nuc_z_size) + " " + str(nuc_x_size/nuc_y_size))
    out_file_nuc.write("\n")
    out_file_nuc.close()
    lines = []
    out_file_diff = open(directory_data + "/difference_" + str(sim_no) + ".dat", "a")
    out_file_diff.write(str(sph_position[0])+ " " + str(sph_pos_bounds[1]-nuc_pos_bounds[1]) + " " + str(sph_x_size/sph_y_size))
    out_file_diff.write("\n")
    out_file_diff.close()
    for pair in pairs:
        pPos = nucleus.mesh.points[pair[0]].get_pos()
        pMembranePos = membrane.mesh.points[pair[1]].get_pos()
        tmpConnection = pMembranePos - pPos
        dlzka = oif.norm(tmpConnection)
        average_dlzka = average_dlzka+dlzka
        if dlzka > max_dlzka:
            max_dlzka = dlzka
        if dlzka < min_dlzka:
            min_dlzka = dlzka
        tmplist = [pPos[0], pPos[1], pPos[2], pMembranePos[0], pMembranePos[1], pMembranePos[2]]
        lines.append(tmplist)
    average_dlzka = average_dlzka/countbonds
    oif.output_vtk_lines(lines=lines, out_file=directory_vtk + "/lines_" + str(cycle//integr_steps) + ".vtk")
    out_file = open(directory_data + "/bonds_info_" + str(sim_no) + ".dat", "a")
    out_file.write(str(cycle) + " " + str(min_dlzka) + " " + str(max_dlzka) + " " + str(average_dlzka))
    out_file.write("\n")
    out_file.close()
    if sph_pos_bounds[1] > 2 * boxX:
        cycle = maxCycle
    cycle = cycle + integr_steps
    max_dlzka = - large_number
    min_dlzka = large_number
    average_dlzka = 0
    #   print "time: ", str(i*time_step*integr_steps)
