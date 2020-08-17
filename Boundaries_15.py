import espressomd
import object_in_fluid as oif

from espressomd import shapes

# GEOMETRY
# rectangular channel with size boxX*boxY*boxZ, the with of the borders is set in width
# and narrow rectangle part in the middle with cross-section ay*az and length ax
boxX = 80.0
boxY = 40.0
boxZ = 40.0
width = 2.0
# ax = par_ax
# ay = par_ay
# az = par_az

# center of the channel
# cx = boxX/2
# cy = boxY/2
# cz = boxZ/2

# size of obstacles that create the narrow part
# horizontal
# ah = ax
# bh = 2*width + boxY
# ch = width + (boxZ - az)/2
# vertical
# av = ax
# bv = width + (boxY - ay)/2
# cv = 2*width + boxZ

# corners for obstacles
# A
# A_x = cx - ax/2
# A_y = -width
# A_z = -width
# B
# B_x = cx - ax/2
# B_y = (boxY + ay)/2
# B_z = -width
# C
# C_x = cx - ax/2
# C_y = -width
# C_z = (boxZ + az)/2

def FillBoundaries(boundaries, vtk_directory, ax, ay, az):
    # GEOMETRY
    # rectangular channel with size boxX*boxY*boxZ, the with of the borders is set in width
    # and narrow rectangle part in the middle with cross-section ay*az and length ax
    # boxX = 80.0
    # boxY = 40.0
    # boxZ = 40.0
    # width = 2.0
    # ax = par_ax
    # ay = par_ay
    # az = par_az

    # center of the channel
    cx = boxX / 2
    cy = boxY / 2
    cz = boxZ / 2

    # size of obstacles that create the narrow part
    # horizontal
    ah = ax
    bh = 2 * width + boxY
    ch = width + (boxZ - az) / 2
    # vertical
    av = ax
    bv = width + (boxY - ay) / 2
    cv = 2 * width + boxZ

    # corners for obstacles
    # A
    A_x = cx - ax / 2
    A_y = -width
    A_z = -width
    # B
    B_x = cx - ax / 2
    B_y = (boxY + ay) / 2
    B_z = -width
    # C
    C_x = cx - ax / 2
    C_y = -width
    C_z = (boxZ + az) / 2

    wallBottom = shapes.Rhomboid(corner=[-width, -width, -width], a=[boxX+2*width, 0.0, 0.0], b=[0.0, boxY+2*width, 0.0], c=[0.0, 0.0, width], direction=1)
    wallTop = shapes.Rhomboid(corner=[-width, -width, boxZ], a=[boxX+2*width, 0.0, 0.0], b=[0.0, boxY+2*width, 0.0], c=[0.0, 0.0, width], direction=1)
    wallBack = shapes.Rhomboid(corner=[-width, boxY, -width], a=[boxX+2*width, 0.0, 0.0], b=[0.0, width, 0.0], c=[0.0, 0.0, boxZ+2*width], direction=1)
    wallFront = shapes.Rhomboid(corner=[-width, -width, -width], a=[boxX+2*width, 0.0, 0.0], b=[0.0, width, 0.0], c=[0.0, 0.0, boxZ+2*width], direction=1)
    horizBottom = shapes.Rhomboid(corner=[A_x, A_y, A_z], a=[ah, 0, 0],    b=[0, bh, 0], c=[0, 0, ch], direction = 1)
    horizTop = shapes.Rhomboid(corner=[C_x, C_y, C_z], a=[ah, 0, 0],   b=[0, bh, 0], c=[0, 0, ch], direction = 1)
    verticLeft = shapes.Rhomboid(corner=[A_x, A_y, A_z], a=[av, 0, 0], b=[0, bv, 0], c=[0, 0, cv], direction = 1)
    verticRight = shapes.Rhomboid(corner=[B_x, B_y, B_z], a=[av, 0, 0],   b=[0, bv, 0], c=[0, 0, cv], direction = 1)
   
    oif.output_vtk_rhomboid(rhom_shape=wallBottom, out_file=vtk_directory + "/wallBottom.vtk")
    oif.output_vtk_rhomboid(rhom_shape=wallTop, out_file=vtk_directory + "/wallTop.vtk")
    oif.output_vtk_rhomboid(rhom_shape=wallBack, out_file=vtk_directory + "/wallBack.vtk")
    oif.output_vtk_rhomboid(rhom_shape=wallFront, out_file=vtk_directory + "/wallFront.vtk")
    oif.output_vtk_rhomboid(rhom_shape=horizBottom, out_file=vtk_directory + "/horizBottom.vtk")
    oif.output_vtk_rhomboid(rhom_shape=horizTop, out_file=vtk_directory + "/horizTop.vtk")
    oif.output_vtk_rhomboid(rhom_shape=verticLeft, out_file=vtk_directory + "/verticLeft.vtk")
    oif.output_vtk_rhomboid(rhom_shape=verticRight, out_file=vtk_directory + "/verticRight.vtk")
    
    # oif.output_vtk_rhomboid(rhom_shape=wallBottom, out_file="output/vtk/wallBottom.vtk")
    # oif.output_vtk_rhomboid(rhom_shape=wallTop, out_file="output/vtk/wallTop.vtk")
    # oif.output_vtk_rhomboid(rhom_shape=wallBack, out_file="output/vtk/wallBack.vtk")
    # oif.output_vtk_rhomboid(rhom_shape=wallFront, out_file="output/vtk/wallFront.vtk")
    # oif.output_vtk_rhomboid(rhom_shape=horizBottom, out_file="output/vtk/horizBottom.vtk")
    # oif.output_vtk_rhomboid(rhom_shape=horizTop, out_file="output/vtk/horizTop.vtk")
    # oif.output_vtk_rhomboid(rhom_shape=verticLeft, out_file="output/vtk/verticLeft.vtk")
    # oif.output_vtk_rhomboid(rhom_shape=verticRight, out_file="output/vtk/verticRight.vtk")
    
    boundaries.append(wallBottom)
    boundaries.append(wallTop)
    boundaries.append(wallBack)
    boundaries.append(wallFront)
    boundaries.append(horizBottom)
    boundaries.append(horizTop)
    boundaries.append(verticLeft)
    boundaries.append(verticRight)
        
    return 0
