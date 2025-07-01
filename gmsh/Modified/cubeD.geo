Lc = 0.5;
Point(1) = {0, 0, 0, Lc};
Point(2) = {1, 0, 0, Lc};
Point(3) = {0, 1, 0, Lc};
Point(4) = {1, 1, 0, Lc};

Line(1) = {1, 3};
Line(2) = {3, 4};
Line(3) = {4, 2};
Line(4) = {2, 1};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};

Extrude {0, 0, 1} {
  Surface{6};
}

// Additional Lines needed for construction

// Creation considering the loops possible
Line Loop(7)  = {22, 11, -13, -2};
Line Loop(8)  = {13, 8, -14, -3};
Line Loop(9)  = {14, 9, -18, -4};
Line Loop(10) = {18, 10, -22, -1};
Line Loop(11) = {-8, -9, -10, -11};

// Definition of surfaces based on loops
Plane Surface(12) = {7};
Plane Surface(13) = {8};
Plane Surface(14) = {9};
Plane Surface(16) = {10};
Plane Surface(17) = {11};

// Corrected orientation defined
Surface Loop(29) = {6, 12, 13, 14, 16, 17};

// Definition of volume
Volume(30) = {29};

// Definition of physical surface
Physical Surface(50) = {6, 12, 13, 14, 16, 17};

// Definition of volume
// Volume(30) = {29};



// Earlier creation of surface
// Surface Loop(29) = {23, 6, 15, 19, 28, 27};

//Volume(30) = {29};

//Physical Volume(40) = {1};

// Physical Surface(50) = {19, 28, 27, 6, 15, 23};

// Physical Surface(44) = {15, 27, 28};
// Physical Surface(45) = {23, 6, 19};

// Program additional changes
Coherence;
