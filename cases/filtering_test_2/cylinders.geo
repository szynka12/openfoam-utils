// Gmsh project created on Tue Feb  9 13:15:36 2021
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 4, 1, 0};
//+
Circle(5) = {0.5, 1/3, 0, 0.1, 0, 2*Pi};
//+
Circle(6) = {0.5, 2/3, 0, 0.1, 0, 2*Pi};
//+
Circle(7) = {0.8, 2/3, 0, 0.1, 0, 2*Pi};
//+
Circle(8) = {0.8, 1/3, 0, 0.1, 0, 2*Pi};

//+
Curve Loop(2) = {6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {5};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {8};
//+
Plane Surface(5) = {5};
//+
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Surface{3}; Surface{5}; Surface{4}; Delete; }
//+
Transfinite Curve {6, 7, 5, 8} = 75 Using Progression 1;
//+
Extrude {0, 0, 1} {
  Surface{1}; Layers{1}; Recombine;
}
//+
Physical Surface("walls") = {5, 2};
//+
Physical Surface("inlet") = {3};
//+
Physical Surface("frontAndBack") = {10, 1};
//+
Physical Surface("cylinders") = {9, 8, 6, 7};
//+
Physical Surface("outlet") = {4};
//+
Physical Volume("internal") = {1};
