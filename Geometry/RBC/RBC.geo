// Gmsh project created on Fri Jun 12 09:41:24 2020
// Gmsh project created on Mon Apr 27 08:53:20 2020
//+

sizemesh = 0.3;

l = 0.5;
L = 1;

coeff = 6;

plate_thickness = 0.02;

Point(1) = {-plate_thickness, l, l, sizemesh/coeff};
//+
Point(2) = {-plate_thickness, -l, l, sizemesh/coeff};
//+
Point(3) = {-plate_thickness, -l, -l, sizemesh/coeff};
//+
Point(4) = {-plate_thickness, +l, -l, sizemesh/coeff};
//+

//+
Point(21) = {0, l, l, sizemesh/coeff};
//+
Point(22) = {0, -l, l, sizemesh/coeff};
//+
Point(23) = {0, -l, -l, sizemesh/coeff};
//+
Point(24) = {0, l, -l, sizemesh/coeff};
//+


//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+

Line(21) = {24, 23};
//+
Line(22) = {23, 22};
//+
Line(23) = {22, 21};
//+
Line(24) = {21, 24};



//+
Point(25) = {0, L, L, sizemesh};
//+
Point(26) = {0, L, -L, sizemesh};
//+
Point(27) = {0, -L, -L, sizemesh};
//+
Point(28) = {0, -L, L, sizemesh};
//+
Point(29) = {2*L, L, L, sizemesh};
//+
Point(30) = {2*L, -L, L, sizemesh};
//+
Point(31) = {2*L, -L, -L, sizemesh};
//+
Point(32) = {2*L, L, -L, sizemesh};
//+
Line(25) = {26, 25};
//+
Line(26) = {25, 29};
//+
Line(27) = {29, 32};
//+
Line(28) = {32, 26};
//+
Line(29) = {26, 27};
//+
Line(30) = {27, 28};
//+
Line(31) = {28, 25};
//+
Line(32) = {28, 30};
//+
Line(33) = {30, 29};
//+
Line(34) = {30, 31};
//+
Line(35) = {31, 32};
//+
Line(36) = {31, 27};
//+
Line(37) = {21, 1};
//+
Line(38) = {24, 4};
//+
Line(39) = {23, 3};
//+
Line(40) = {22, 2};
//+
Line Loop(1) = {24, 21, 22, 23};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {2, 3, 4, 1};
//+
Plane Surface(2) = {2};
//+
Line Loop(4) = {24, 38, 4, -37};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {21, 39, 3, -38};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {22, 40, 2, -39};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {23, 37, 1, -40};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {25, 26, 27, 28};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {30, 32, 34, 36};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {34, 35, -27, -33};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {31, 26, -33, -32};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {28, 29, -36, 35};
//+
Plane Surface(12) = {12};
//+
Surface Loop(1) = {1, 4, 5, 6, 7, 2};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {3, 12, 8, 11, 10, 9, 2};
//+
Line Loop(13) = {25, -31, -30, -29};
//+
Plane Surface(13) = {1, 13};
//+
Surface Loop(3) = {13, 12, 8, 11, 10, 9, 1};
//+
Volume(2) = {3};
//+
Physical Surface("RBC_boundary",6) = {8, 11, 10, 9, 12};
//+
Physical Surface("coupling",3) = {1};
//+
Physical Surface("ext_pressure",4) = {2};
//+
Physical Surface("embedding",5) = {4, 5, 6, 7};
//+
Physical Volume("plate",2) = {1};
//+
Physical Volume("acoustic",1) = {2};
