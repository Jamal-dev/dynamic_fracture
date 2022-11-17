Mesh.SaveParametric = 1;
Mesh.Algorithm = 8;

Point(1) = {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e-01};
Point(2) = {0.000000e+00, 4.950000e+00, 0.000000e+00, 1.000000e-01};
Point(3) = {5.000000e+00, 5.000000e+00, 0.000000e+00, 1.000000e-01};
Point(4) = {0.000000e+00, 5.050000e+00, 0.000000e+00, 1.000000e-01};
Point(5) = {0.000000e+00, 1.000000e+01, 0.000000e+00, 1.000000e-01};
Point(6) = {1.000000e+01, 1.000000e+01, 0.000000e+00, 1.000000e-01};
Point(7) = {1.000000e+01, 0.000000e+00, 0.000000e+00, 1.000000e-01};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
Line Loop(1) = {5, 6, 7, 1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Surface(1) = {1};


Recombine Surface{1};
// right edge
Physical Curve(0) = {6};
// left edge
Physical Curve( 1) = {1, 4};
// bottom edge
Physical Curve( 2) = {7};
// top edge
Physical Curve( 3) = {5};
// bottom_crack
Physical Curve( 4) = {2};
// bottom_crack
Physical Curve( 5) = {3};


