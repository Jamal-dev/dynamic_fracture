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
def in2m(inches):
    return inches * 0.0254
    
# Initialize gmsh:
gmsh.initialize()



# square points:

lc = in2m(0.1)
geo_eps = in2m(0.005)
width = in2m(20.0)
height = in2m(8.0)
radius = in2m(0.5)
a = in2m(1)
num_holes = 3
b = in2m(6)
dist_bw_holes = in2m(2)
dist_first_hole = in2m(1.5)
dist_hole_x = width/2 - b + in2m(a)
dist_a_x = width/2 - b
# both supposrts
left_supp_x = in2m(1)
right_supp_x = in2m(19)

circ_x = [dist_hole_x for _ in range(num_holes)]
circ_y = [height - dist_first_hole - i*dist_bw_holes for i in range(num_holes)]

def point_circle(theta, circ_num):
    return circ_x[circ_num] + radius*cosd(theta), circ_y[circ_num] + radius*sind(theta)
def add_point_circle(theta,circ_num):
    x, y = point_circle(theta,circ_num)
    print(f'theta = {theta}: Adding point ({x} , {y})')
    return gmsh.model.geo.addPoint(x, y, 0, lc)


corner_pnt1 = gmsh.model.geo.add_point(0.0, height, 0, lc)
a_pnt1 = gmsh.model.geo.add_point(dist_a_x, height, 0, lc)
if num_holes:
    top_circle_pnt1 = gmsh.model.geo.add_point(dist_hole_x, height, 0, lc)
mid_pnt1 = gmsh.model.geo.add_point(width/2, height, 0, lc)
corner_pnt2 = gmsh.model.geo.add_point(width, height, 0, lc)
right_circle_pnt1 = gmsh.model.geo.add_point(width, height - dist_first_hole, 0, lc)
right_circle_pnt2 = gmsh.model.geo.add_point(width, height - dist_first_hole - dist_bw_holes, 0, lc)
right_circle_pnt3 = gmsh.model.geo.add_point(width, height - dist_first_hole - 2*dist_bw_holes, 0, lc)
a_pnt2 = gmsh.model.geo.add_point(width , a, 0, lc)
corner_pnt3 = gmsh.model.geo.add_point(width, 0.0, 0, lc)
mid_pnt3 = gmsh.model.geo.add_point(width/2, 0.0, 0, lc)
if num_holes:
    bottom_circle_pnt1 = gmsh.model.geo.add_point(dist_hole_x, 0.0, 0, lc)

a_pnt3 = gmsh.model.geo.add_point(dist_a_x+geo_eps, 0.0, 0, lc)
a_pnt4 = gmsh.model.geo.add_point(dist_a_x, a, 0, lc)
a_pnt5 = gmsh.model.geo.add_point(dist_a_x-geo_eps, 0.0, 0, lc)

corner_pnt4 = gmsh.model.geo.add_point(0.0, 0.0, 0, lc)
a_pnt6 = gmsh.model.geo.add_point(0.0, a, 0, lc)
left_circle_pnt1 = gmsh.model.geo.add_point(0.0, height - dist_first_hole - 2*dist_bw_holes, 0, lc)
left_circle_pnt2 = gmsh.model.geo.add_point(0.0, height - dist_first_hole - dist_bw_holes, 0, lc)
left_circle_pnt3 = gmsh.model.geo.add_point(0.0, height - dist_first_hole, 0, lc)
center_pnts = [gmsh.model.geo.add_point(x, y , 0, lc) for x,y in zip(circ_x, circ_y)]

def add_mesh_help_pnts(x, num_holes):
    pnts = [gmsh.model.geo.add_point(x, height - dist_first_hole - i*dist_bw_holes, 0, lc) for i in range(num_holes)]
    return pnts
add_mesh_help_pnts(dist_a_x, num_holes=3)
add_mesh_help_pnts(width/2, num_holes=3)
gmsh.model.geo.add_point(width/2, a, 0, lc)

def add_arc(num_pnts=4):
    theta = [i*360/num_pnts for i in range(num_pnts)]
    circ_edges = []
    circ_loops = []
    for circ_num in range(num_holes):
        pnts = [add_point_circle(t,circ_num) for t in theta]
        pnts_temp = pnts.copy()
        pnts_temp.append(pnts[0])
        arc_edges = []
        for pn1, pn2 in zip(pnts, pnts_temp[1:]):
            arc_edges.append( gmsh.model.geo.addCircleArc(pn1, center_pnts[circ_num], pn2))
        circ_edges.append(arc_edges)
        circ_loops.append(gmsh.model.geo.add_curve_loop(arc_edges))

    return circ_edges, circ_loops

# xc,yc = hor_intersect(circ_x, circ_y, width - notch_height/2, height,notch_center_y)
# virtual_pnt1 = gmsh.model.geo.add_point(xc,yc, 0, lc)
# xc,yc = hor_intersect(circ_x, circ_y, width, height - notch_vert_dist,height - notch_vert_dist)
# virtual_pnt2 = gmsh.model.geo.add_point(xc,yc, 0, lc)

# adding circle arcs
num_pnts = 4
circ_edges, circ_loops = add_arc(num_pnts)
print(f'circ_edges = {circ_edges}')

# Edge of square:
line1 = gmsh.model.geo.add_line(corner_pnt1, a_pnt1)
if num_holes:
    line2 = gmsh.model.geo.add_line(a_pnt1, top_circle_pnt1)
    line3 = gmsh.model.geo.add_line(top_circle_pnt1, mid_pnt1)
else:
    line2 = gmsh.model.geo.add_line(a_pnt1, mid_pnt1)
line4 = gmsh.model.geo.add_line(mid_pnt1, corner_pnt2)
line5 = gmsh.model.geo.add_line(corner_pnt2, right_circle_pnt1)
line6 = gmsh.model.geo.add_line(right_circle_pnt1, right_circle_pnt2)
line7 = gmsh.model.geo.add_line(right_circle_pnt2, right_circle_pnt3)
line8 = gmsh.model.geo.add_line(right_circle_pnt3, a_pnt2)
line9 = gmsh.model.geo.add_line(a_pnt2, corner_pnt3)

line10 = gmsh.model.geo.add_line(corner_pnt3, mid_pnt3)
if num_holes:
    line11 = gmsh.model.geo.add_line(mid_pnt3, bottom_circle_pnt1)
    line12 = gmsh.model.geo.add_line(bottom_circle_pnt1, a_pnt3)
else:
    line11 = gmsh.model.geo.add_line(mid_pnt3, a_pnt3)
line13 = gmsh.model.geo.add_line(a_pnt3, a_pnt4)
line14 = gmsh.model.geo.add_line(a_pnt4, a_pnt5)
line15 = gmsh.model.geo.add_line(a_pnt5, corner_pnt4)

line16 = gmsh.model.geo.add_line(corner_pnt4, a_pnt6)
line17 = gmsh.model.geo.add_line(a_pnt6, left_circle_pnt1)
line18 = gmsh.model.geo.add_line(left_circle_pnt1, left_circle_pnt2)
line19 = gmsh.model.geo.add_line(left_circle_pnt2, left_circle_pnt3)
line20 = gmsh.model.geo.add_line(left_circle_pnt3, corner_pnt1)


# Curve loop:
# square:
if num_holes:
    curve_loop1 = gmsh.model.geo.add_curve_loop([line1, line2, \
                                                line3, line4, \
                                                line5, line6, \
                                                line7, line8, \
                                                line9, line10, \
                                                line11, line12, \
                                                line13, line14, \
                                                line15, line16, \
                                                line17, line18, \
                                                line19, line20])
else:
    curve_loop1 = gmsh.model.geo.add_curve_loop([line1, line2, \
                                                 line4, \
                                                line5, line6, \
                                                line7, line8, \
                                                line9, line10, \
                                                line11,  \
                                                line13, line14, \
                                                line15, line16, \
                                                line17, line18, \
                                                line19, line20])

# plane surface:
circ_loops.insert(0, curve_loop1)
plane_surface1 = gmsh.model.geo.add_plane_surface(circ_loops)

# top edge 
if num_holes:
    top_edge = [line1, line2, line3, line4]
else:
    top_edge = [line1, line2, line4]
# right edge
right_edge = [line5, line6, line7, line8, line9]
# bottom edge
if num_holes:
    bottom_edge = [line10, line11, line12,line15]
else:
    bottom_edge = [line10, line11, line15]
# slit edge
slit_edge = [line13, line14]
# left edge
left_edge = [line16,  line17, line18, line19, line20]


# adding physical groups:
gmsh.model.add_physical_group(1,left_edge, 1)
gmsh.model.add_physical_group(1,right_edge, 0)
gmsh.model.add_physical_group(1,top_edge, 3)
gmsh.model.add_physical_group(1,bottom_edge, 2)
gmsh.model.add_physical_group(1,slit_edge, 4)
circ_phy_group = [gmsh.model.add_physical_group(1,c_e, 5+i) for i,c_e in enumerate(circ_edges)]
gmsh.model.add_physical_group(2,[plane_surface1], 20)

# from Gmsh model.
gmsh.model.geo.synchronize()
# Generate mesh:

gmsh.option.setNumber('Mesh.Algorithm', 8)
gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 2)
gmsh.option.setNumber('Mesh.RecombineAll', 1)
gmsh.option.setNumber('Mesh.MshFileVersion', 4)

gmsh.model.mesh.generate()
# Write mesh data:
gmsh.write(f"p_compression_test_num_holes_{num_holes}.msh")
# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()
# It finalize t
gmsh.finalize()
