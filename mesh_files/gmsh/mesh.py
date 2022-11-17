# Import modules:
import gmsh
import sys


# Initialize gmsh:
gmsh.initialize()



# square points:
lc = 1e-1
eps = 0.05
width = 10.0
height = 10.0
pnt1 = gmsh.model.geo.add_point(0.0, height, 0, lc)
pnt2 = gmsh.model.geo.add_point(width/2, height, 0, lc)
pnt3 = gmsh.model.geo.add_point(width, height, 0, lc)
pnt4 = gmsh.model.geo.add_point(width, height/2, 0, lc)
pnt5 = gmsh.model.geo.add_point(width, 0.0, 0, lc)
pnt6 = gmsh.model.geo.add_point(height/2, 0.0, 0, lc)
pnt7 = gmsh.model.geo.add_point(0.0, 0.0, 0, lc)
pnt8 = gmsh.model.geo.add_point(0.0, height/2 - eps, 0, lc)
pnt9 = gmsh.model.geo.add_point(width/2, height/2 , 0, lc)
pnt10 = gmsh.model.geo.add_point(0.0, height/2 + eps, 0, lc)

# Edge of square:
line1 = gmsh.model.geo.add_line(pnt1, pnt2)
line2 = gmsh.model.geo.add_line(pnt2, pnt3)
line3 = gmsh.model.geo.add_line(pnt3, pnt4)
line4 = gmsh.model.geo.add_line(pnt4, pnt5)
line5 = gmsh.model.geo.add_line(pnt5, pnt6)
line6 = gmsh.model.geo.add_line(pnt6, pnt7)
line7 = gmsh.model.geo.add_line(pnt7, pnt8)
line8 = gmsh.model.geo.add_line(pnt8, pnt9)
line9 = gmsh.model.geo.add_line(pnt9, pnt10)
line10 = gmsh.model.geo.add_line(pnt10, pnt1)

# Curve loop:
curve_loop1 = gmsh.model.geo.add_curve_loop([line1, line2, \
                                            line3, line4, \
                                            line5, line6, \
                                            line7, line8, \
                                            line9, line10])
# plane surface:
plane_surface1 = gmsh.model.geo.add_plane_surface([curve_loop1])

# top edge 
top_edge = [line1, line2]
# right edge
right_edge = [line3, line4]
# bottom edge
bottom_edge = [line5, line6]
# left edge
left_edge = [line7,  line10]
# crack bottom
crack_bottom = [line8]
# crack top
crack_top = [line9]

# adding physical groups:
gmsh.model.add_physical_group(1,left_edge, 1)
gmsh.model.add_physical_group(1,right_edge, 0)
gmsh.model.add_physical_group(1,top_edge, 3)
gmsh.model.add_physical_group(1,bottom_edge, 2)
gmsh.model.add_physical_group(1,crack_bottom, 4)
gmsh.model.add_physical_group(1,crack_top, 5)
gmsh.model.add_physical_group(2,[plane_surface1], 6)

# from Gmsh model.
gmsh.model.geo.synchronize()
# Generate mesh:

gmsh.option.setNumber('Mesh.Algorithm', 8)
gmsh.option.setNumber('Mesh.RecombinationAlgorithm', 1)
gmsh.option.setNumber('Mesh.RecombineAll', 1)
gmsh.option.setNumber('Mesh.MshFileVersion', 4)

gmsh.model.mesh.generate()
# Write mesh data:
gmsh.write("mesh.msh")
# Creates  graphical user interface
if 'close' not in sys.argv:
    gmsh.fltk.run()
# It finalize t
gmsh.finalize()