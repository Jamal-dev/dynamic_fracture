# Import modules:
import gmsh
import sys
import numpy as np
def cosd(x):
    return np.cos(x*np.pi/180)
def sind(x):
    return np.sin(x*np.pi/180)
def angle_by_points(x1, y1, x2, y2):
    return np.arctan2(y2-y1, x2-x1)*180/np.pi
def hor_intersect(x1,y1,x2,y2,yc):
    m = (y2-y1)/(x2-x1)
    def line(y):
        return (y-y1)/m + x1
    xc = line(yc)
    return (xc, yc)
    
# Initialize gmsh:
gmsh.initialize()



# square points:

lc = 1
eps = 0.05
width = 10.0
height = 10.0
radius = 1.6
notch_vert_dist = 1.0
notch_height = 3.0
circ_x = width/2.0
circ_y = height/2.0 - 0.5
notch_center_x = width - notch_height/2
notch_center_y = height - notch_vert_dist - notch_height/2
def point_circle(theta):
    return circ_x + radius*cosd(theta), circ_y + radius*sind(theta)
def add_point_circle(theta):
    x, y = point_circle(theta)
    print(f'theta = {theta}: Adding point ({x} , {y})')
    return gmsh.model.geo.addPoint(x, y, 0, lc)

corner_pnt1 = gmsh.model.geo.add_point(0.0, height, 0, lc)
mid_pnt1 = gmsh.model.geo.add_point(width/2, height, 0, lc)
top_notch_pnt1 = gmsh.model.geo.add_point(width - notch_height/2, height , 0, lc)
corner_pnt2 = gmsh.model.geo.add_point(width, height, 0, lc)
notch_pnt1 = gmsh.model.geo.add_point(width, height - notch_vert_dist, 0, lc)
notch_pnt2 = gmsh.model.geo.add_point(notch_center_x, notch_center_y, 0, lc)
notch_pnt3 = gmsh.model.geo.add_point(width , height - notch_vert_dist - notch_height, 0, lc)
mid_pnt2 = gmsh.model.geo.add_point(width, circ_y, 0, lc)
corner_pnt3 = gmsh.model.geo.add_point(width, 0.0, 0, lc)
mid_pnt3 = gmsh.model.geo.add_point(height/2, 0.0, 0, lc)
corner_pnt4 = gmsh.model.geo.add_point(0.0, 0.0, 0, lc)
mid_pnt4 = gmsh.model.geo.add_point(0.0, circ_y, 0, lc)
center_pnt = gmsh.model.geo.add_point(circ_x, circ_y , 0, lc)

# xc,yc = hor_intersect(circ_x, circ_y, width - notch_height/2, height,notch_center_y)
# virtual_pnt1 = gmsh.model.geo.add_point(xc,yc, 0, lc)
# xc,yc = hor_intersect(circ_x, circ_y, width, height - notch_vert_dist,height - notch_vert_dist)
# virtual_pnt2 = gmsh.model.geo.add_point(xc,yc, 0, lc)

ang1 = angle_by_points(circ_x, circ_y, notch_center_x, notch_center_y)
ang2 = angle_by_points(circ_x, circ_y, width - notch_height/2, height)
theta3 = ang2 - ang1
theta = 90
pnt1 = add_point_circle(theta)

theta1 = (90 - 2 * theta3 )/2
theta -= theta1
pnt2 = add_point_circle(theta)
theta -= theta3
pnt3 = add_point_circle(theta)
theta -= theta3
pnt4 = add_point_circle(theta)
theta -= theta1
pnt5 = add_point_circle(theta)
theta = -45
pnt6 = add_point_circle(theta)
theta = -90
pnt7 = add_point_circle(theta)
theta = -135
pnt8 = add_point_circle(theta)
theta = -180
pnt9 = add_point_circle(theta)
theta = -225
pnt10 = add_point_circle(theta)

# Edge of circle
circ1 = gmsh.model.geo.addCircleArc(pnt1, center_pnt, pnt2) # 90-60
circ2 = gmsh.model.geo.addCircleArc(pnt2, center_pnt, pnt3) # 60-45
circ3 = gmsh.model.geo.addCircleArc(pnt3, center_pnt, pnt4) # 45-30
circ4 = gmsh.model.geo.addCircleArc(pnt4, center_pnt, pnt5) # 30-0
circ5 = gmsh.model.geo.addCircleArc(pnt5, center_pnt, pnt6)
circ6 = gmsh.model.geo.addCircleArc(pnt6, center_pnt, pnt7)
circ7 = gmsh.model.geo.addCircleArc(pnt7, center_pnt, pnt8)
circ8 = gmsh.model.geo.addCircleArc(pnt8, center_pnt, pnt9)
circ9 = gmsh.model.geo.addCircleArc(pnt9, center_pnt, pnt10)
circ10 = gmsh.model.geo.addCircleArc(pnt10, center_pnt, pnt1)

# Edge of square:
line1 = gmsh.model.geo.add_line(corner_pnt1, mid_pnt1)
line2 = gmsh.model.geo.add_line(mid_pnt1, top_notch_pnt1)
line3 = gmsh.model.geo.add_line(top_notch_pnt1, corner_pnt2)
line4 = gmsh.model.geo.add_line(corner_pnt2, notch_pnt1)
line5 = gmsh.model.geo.add_line(notch_pnt1, notch_pnt2)
line6 = gmsh.model.geo.add_line(notch_pnt2, notch_pnt3)
line7 = gmsh.model.geo.add_line(notch_pnt3, mid_pnt2)
line8 = gmsh.model.geo.add_line(mid_pnt2, corner_pnt3)
line9 = gmsh.model.geo.add_line(corner_pnt3, mid_pnt3)
line10 = gmsh.model.geo.add_line(mid_pnt3, corner_pnt4)
line11 = gmsh.model.geo.add_line(corner_pnt4, mid_pnt4)
line12 = gmsh.model.geo.add_line(mid_pnt4, corner_pnt1)

# Curve loop:
# square:
curve_loop1 = gmsh.model.geo.add_curve_loop([line1, line2, \
                                            line3, line4, \
                                            line5, line6, \
                                            line7, line8, \
                                            line9, line10, \
                                            line11, line12])
# circle:
curve_loop2 = gmsh.model.geo.add_curve_loop([circ1, circ2, \
                                            circ3, circ4, \
                                            circ5, circ6, \
                                            circ7, circ8, \
                                            circ9, circ10])
# plane surface:
plane_surface1 = gmsh.model.geo.add_plane_surface([curve_loop1, curve_loop2])

# top edge 
top_edge = [line1, line2, line3]
# right edge
right_edge = [line4, line5, line6, line7, line8]
# bottom edge
bottom_edge = [line9, line10]
# left edge
left_edge = [line11,  line12]
# notch
# notch_edge = [line5, line6]
# circle edge
circle_edge = [circ1, circ2, circ3, circ4, circ5, circ6, circ7, circ8, circ9, circ10]


# adding physical groups:
gmsh.model.add_physical_group(1,left_edge, 1)
gmsh.model.add_physical_group(1,right_edge, 0)
gmsh.model.add_physical_group(1,top_edge, 3)
gmsh.model.add_physical_group(1,bottom_edge, 2)
# gmsh.model.add_physical_group(1,notch_edge, 4)
gmsh.model.add_physical_group(1,circle_edge, 5)
gmsh.model.add_physical_group(2,[plane_surface1], 6)

# from Gmsh model.
gmsh.model.geo.synchronize()
# Generate mesh:

gmsh.option.setNumber('Mesh.Algorithm', 8)
gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 2)
gmsh.option.setNumber('Mesh.RecombineAll', 1)
gmsh.option.setNumber('Mesh.MshFileVersion', 4)

gmsh.model.mesh.generate()
# Write mesh data:
gmsh.write("p_notched_cavity.msh")
# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()
# It finalize t
gmsh.finalize()
