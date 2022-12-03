Mesh.SaveParametric = 1;
Mesh.Algorithm = 8;

Point(1) = {0, 0, 0, 1e+22};
Point(2) = {10, 0, 0, 1e+22};
Point(3) = {10, 6, 0, 1e+22};
Point(4) = {10, 9, 0, 1e+22};
Point(5) = {10, 10, 0, 1e+22};
Point(7) = {0, 10, 0, 1e+22};
Point(6) = {8.5, 7.5, 0, 1e+22};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 6};
//+
Line(4) = {6, 4};
//+
Line(5) = {4, 5};
//+
Line(6) = {5, 7};
//+
Line(7) = {7, 1};

//+
Point(8) = {5, 4.5, 0, 1.0};
//+
Point(9) = {5, 3, 0, 1.0};
//+
Point(10) = {5, 6, 0, 1.0};
//+
Point(11) = {6.5, 4.5, 0, 1.0};
//+
Point(12) = {3.5, 4.5, 0, 1.0};
//+
Circle(8) = {10, 8, 12};
//+
Circle(9) = {11, 8, 10};
//+
Circle(10) = {12, 8, 9};
//+
Circle(11) = {9, 8, 11};
//+
Line Loop(1) = {6, 7, 1, 2, 3, 4, 5};
//+
Line Loop(2) = {9, 8, 10, 11};
//+
Plane Surface(1) = {1, 2};
//+
Physical Surface("boundary") = {1};
//+

Recombine Surface{1};

