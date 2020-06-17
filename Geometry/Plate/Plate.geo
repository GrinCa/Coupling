// Gmsh project created on Mon Apr 27 08:53:20 2020
//+

sizemesh = 0.6;

l = 0.5;
L = 1;
plate_thickness = 0.02;

coeff = 7;


Point(1) = {plate_thickness, l, l, sizemesh/coeff};
//+
Point(2) = {plate_thickness, -l, l, sizemesh/coeff};
//+
Point(3) = {plate_thickness, -l, -l, sizemesh/coeff};
//+
Point(4) = {plate_thickness, +l, -l, sizemesh/coeff};
//+
Point(13) = {plate_thickness, L, L, sizemesh};
//+
Point(14) = {plate_thickness, -L, L, sizemesh};
//+
Point(15) = {plate_thickness, -L, -L, sizemesh};
//+
Point(16) = {plate_thickness, L, -L, sizemesh};

Point(17) = {2*L, L, L, sizemesh};
//+
Point(18) = {2*L, -L, L, sizemesh};
//+
Point(19) = {2*L, -L, -L, sizemesh};
//+
Point(20) = {2*L, L, -L, sizemesh};
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
Line(5) = {13, 14};
//+
Line(6) = {14, 15};
//+
Line(7) = {15, 16};
//+
Line(8) = {16, 13};
//+

Line(9) = {13, 17};
//+
Line(10) = {17, 20};
//+
Line(11) = {20, 16};
//+
Line(12) = {14, 18};
//+
Line(13) = {18, 19};
//+
Line(14) = {19, 15};
//+
Line(15) = {18, 17};
//+
Line(16) = {19, 20};
//+
Line(17) = {1, 21};
//+
Line(18) = {2, 22};
//+
Line(19) = {3, 23};
//+
Line(20) = {4, 24};
//+
Line(21) = {24, 23};
//+
Line(22) = {23, 22};
//+
Line(23) = {22, 21};
//+
Line(24) = {21, 24};


//+
Line Loop(1) = {9, 10, 11, 8};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {5, 12, 15, -9};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {13, 16, -10, -15};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {11, -7, -14, 16};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {5, 6, 7, 8};
//+
Line Loop(6) = {1, 2, 3, 4};
//+
Plane Surface(5) = {5, 6};
//+
Line Loop(7) = {12, 13, 14, -6};
//+
Plane Surface(6) = {7};
//+
Plane Surface(7) = {6};
//+
Line Loop(8) = {22, 23, 24, 21};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {17, -23, -18, -1};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {18, -22, -19, -2};
//+
Plane Surface(10) = {10};
//+
Line Loop(11) = {19, -21, -20, -3};
//+
Plane Surface(11) = {11};
//+
Line Loop(12) = {4, 17, 24, -20};
//+
Plane Surface(12) = {12};
//+
Surface Loop(1) = {1, 2, 5, 6, 3, 4, 7};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {8, 10, 9, 12, 11, 7};
//+
Volume(2) = {2};
//+

//+
Physical Volume("cavity",1) = {1};
//+
Physical Volume("plate",2) = {2};
//+
Physical Surface("embedding",3) = {9, 12, 11, 10};
//+
Physical Surface("coupling",4) = {7};
//+
Physical Surface("ext_pressure", 5) = {8};
