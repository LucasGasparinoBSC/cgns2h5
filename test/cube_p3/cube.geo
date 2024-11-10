//pi = 3.14159265;
pi = 3.0;
n = 1;
h = 2*pi/n;

//+ Defiine the geometry
Point(1) = {0.0,0.0,0.0,h};
Point(2) = {2*pi,0.0,0.0,h};
Point(3) = {2*pi,2*pi,0.0,h};
Point(4) = {0.0,2*pi,0.0,h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1};

Extrude {0, 0, 2*pi} {
  Surface{1}; Layers{n}; Recombine;
}
Transfinite Curve {9, 8, 7, 6, 20, 11, 4, 16, 3, 12, 2, 1} = 2 Using Progression 1;

//+ Set every surface to periodic, using surf(1) as master
//Physical Surface("Periodic") = {21,13,17,25,26,1};

//+ Set the fluid region
Physical Volume("fluid") = {1};

//+ Set the output mesh file version
Mesh.MshFileVersion = 2.2;

//+ Options controlling mesh generation
Mesh.ElementOrder = 3; //+ Set desired element order beetween 3 and 8 (1,2 not supported in SOD)
Mesh 3;                //+ Volumetric mesh
//+
