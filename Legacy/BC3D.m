function [uhd,qhn]=BC3D(ux, uy, uz, gx, gy, gz, T, k, formula)

% [uhd,qhn]=BC3d(ux, uy, uz, gx, gy, gz, T, k, formula)
%
% Input:
%    ux,uy,uz: Dirichlet data; vectorized function of three variables
%    gx,gy,gz: Neumann data (corresponds to kappa*grad(u))
%              vectorized functions of three variables
%           T: Full tetrahedrization
%           k: polynomial degree
%     formula: quadrature formula in 2d (N x 4 matrix)
%
% Output:
%         uhd: d2 x Ndir,   Ndir, number of Dirichlet faces
%         qhn: d2 x Nneu,   Nneu, number of Neumann faces
%
% Last modified: March 21, 2013

% General definitions
Ndir   = size(T.dirichlet,1);
d2     = nchoosek(k+2,2);
Nnodes = size(formula,1);

% Separation on components of coordinates
x=T.coordinates(:,1);
y=T.coordinates(:,2);
z=T.coordinates(:,3);

% Dubiner evaluation on reference element
db = dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k); % Nnodes x d2
wD = bsxfun(@times,formula(:,4),db);                 % Nnodes x d2

% ------------------------ Dirichlet BC ------------------------ %

% Dirichlet elements
x12=x(T.dirichlet(:,2))-x(T.dirichlet(:,1)); %x2-x1
y12=y(T.dirichlet(:,2))-y(T.dirichlet(:,1)); %y2-y1
z12=z(T.dirichlet(:,2))-z(T.dirichlet(:,1)); %z2-z1
x13=x(T.dirichlet(:,3))-x(T.dirichlet(:,1)); %x3-x1
y13=y(T.dirichlet(:,3))-y(T.dirichlet(:,1)); %y3-y1
z13=z(T.dirichlet(:,3))-z(T.dirichlet(:,1)); %z3-z1

% Dirichlet normals: Decomposition of components
Dir_normals=0.5*[y12.*z13-z12.*y13,...
                 z12.*x13-x12.*z13,...
                 x12.*y13-x13.*y12];    % Nneu x 3,  ||n^e||=|e|

Nx_D = Dir_normals(:, 1);       % Dirichlet Normals x-component
Ny_D = Dir_normals(:, 2);       % Dirichlet Normals y-component      
Nz_D = Dir_normals(:, 3);       % Dirichlet Normals z-component

% ------- Redefinition of space: dirichlet conditions ------- %

% Let's remember the normals result after the transformation:
% N = (P2 - P1) x (P3 - P1), transform matrix asociated it's.
A_0 = 0.5 * [x12, y12, z12]; A_1 = 0.5 * [x13, y13, z13];

% Dimentional array of all possible transformations per element
A = zeros([3, 2, Ndir]);
A(:, 1, :) = A_1'; A(:, 2, :) = A_0';

% Now we define the new db space for the Dirichlet conditions
% - Dimentional of the basis
db_resized = zeros([2, 2*d2 * Nnodes]);
db_resized(1, 1:d2*Nnodes) = reshape(db, [1, d2 * Nnodes]);
db_resized(2, (d2*Nnodes+1):end) = reshape(db, [1, d2 * Nnodes]);

% - Repetition of the definitition per element
rep_db = repmat(db_resized, Ndir, 1);
rep_db = reshape(rep_db, [2, Ndir, 2*d2*Nnodes]);

% - Dimentional correction to have per each element the 
%   Transposed matrix
rep_db = permute(rep_db, [1, 3, 2]);

% Finally the polynomial is defined as:
db_dimentional = pagemtimes(A, rep_db);

% - Also for each coordiante we can define:
db_dim_x = [reshape(db_dimentional(1, 1:Nnodes*d2, :), [d2, Nnodes, Ndir])       ;...
            reshape(db_dimentional(1, (Nnodes*d2+1):end, :), [d2, Nnodes, Ndir])];
db_dim_y = [reshape(db_dimentional(2, 1:Nnodes*d2, :), [d2, Nnodes, Ndir])       ;
            reshape(db_dimentional(2, (Nnodes*d2+1):end, :), [d2, Nnodes, Ndir])];
db_dim_z = [reshape(db_dimentional(3, 1:Nnodes*d2, :), [d2, Nnodes, Ndir])       ;
            reshape(db_dimentional(3, (Nnodes*d2+1):end, :), [d2, Nnodes, Ndir])];

% -------------- Definition of all products -------------- % 

weights = repmat(formula(:,4)', 1, 1);

% Definition of weights * db coeficients
wD_x = bsxfun(@times,weights,db_dim_x); % Nnodes x d2
wD_y = bsxfun(@times,weights,db_dim_y); % Nnodes x d2
wD_z = bsxfun(@times,weights,db_dim_z); % Nnodes x d2

% Reference element coordiantes in Dirichlet faces times
% formula weights
xd=formula(:,[1 2 3])*x(T.dirichlet');  % Nnodes x Nneu
yd=formula(:,[1 2 3])*y(T.dirichlet');  % Nnodes x Nneu
zd=formula(:,[1 2 3])*z(T.dirichlet');  % Nnodes x Nneu

% Cross product coordinate definition:
cross_u_x = bsxfun(@times, Ny_D', uz(xd, yd, zd)) - ...
                        bsxfun(@times, Nz_D', uy(xd, yd, zd));
cross_u_y = bsxfun(@times, Nz_D', ux(xd, yd, zd)) - ...
                        bsxfun(@times, Nx_D', uz(xd, yd, zd));
cross_u_z = bsxfun(@times, Nx_D', uy(xd, yd, zd)) - ...
                        bsxfun(@times, Ny_D', ux(xd, yd, zd));

% Reshape per element of each coordinate evaluation
cross_u_x = reshape(cross_u_x, [Nnodes, 1, Ndir]);
cross_u_y = reshape(cross_u_y, [Nnodes, 1, Ndir]);
cross_u_z = reshape(cross_u_z, [Nnodes, 1, Ndir]);

% Product per coordinate: Boundary condition
uhd_x = reshape(pagemtimes(wD_x, cross_u_x), [2*d2, 1, Ndir]);
uhd_y = reshape(pagemtimes(wD_y, cross_u_y), [2*d2, 1, Ndir]);
uhd_z = reshape(pagemtimes(wD_z, cross_u_z), [2*d2, 1, Ndir]);

uhd_fun = uhd_x + uhd_y + uhd_z;

% Normals reshaping for compatibility of dimentions
Nx_D_dim = reshape(Nx_D, [1, 1, Ndir]);
Ny_D_dim = reshape(Ny_D, [1, 1, Ndir]);
Nz_D_dim = reshape(Nz_D, [1, 1, Ndir]);

% Left side definition:
cross_test_x = bsxfun(@times, Ny_D_dim, db_dim_z) - ...
                        bsxfun(@times, Nz_D_dim, db_dim_y);
cross_test_y = bsxfun(@times, Nz_D_dim, db_dim_x) - ...
                        bsxfun(@times, Nx_D_dim, db_dim_z);
cross_test_z = bsxfun(@times, Nx_D_dim, db_dim_y) - ...
                        bsxfun(@times, Ny_D_dim, db_dim_x);

% Product per coordinate: Test function condition
testd_x = pagemtimes(wD_x, pagetranspose(cross_test_x));
testd_y = pagemtimes(wD_y, pagetranspose(cross_test_y));
testd_z = pagemtimes(wD_z, pagetranspose(cross_test_z));

testd = testd_x + testd_y + testd_z;

% Final product construction
uhd = pagemtimes(pageinv(testd), uhd_fun);

% ------------------------ Neumann BC ------------------------ %

% Neumman elements
x12=x(T.neumann(:,2))-x(T.neumann(:,1)); %x2-x1
y12=y(T.neumann(:,2))-y(T.neumann(:,1)); %y2-y1
z12=z(T.neumann(:,2))-z(T.neumann(:,1)); %z2-z1
x13=x(T.neumann(:,3))-x(T.neumann(:,1)); %x3-x1
y13=y(T.neumann(:,3))-y(T.neumann(:,1)); %y3-y1
z13=z(T.neumann(:,3))-z(T.neumann(:,1)); %z3-z1

Neu_normals=0.5*[y12.*z13-z12.*y13,... % Neumann normals
                 z12.*x13-x12.*z13,...
                 x12.*y13-x13.*y12];   % Nneu x 3,||n^e||=|e|

% Reference element coordiantes in Neumann faces times
% formula weights
xn=formula(:,[1 2 3])*x(T.neumann');  % Nnodes x Nneu
yn=formula(:,[1 2 3])*y(T.neumann');  % Nnodes x Nneu
zn=formula(:,[1 2 3])*z(T.neumann');  % Nnodes x Nneu

% Integral quadrature definition per components.
qhn=bsxfun(@times,Neu_normals(:,1)',wD'*gx(xn,yn,zn))+...
    bsxfun(@times,Neu_normals(:,2)',wD'*gy(xn,yn,zn))+...
    bsxfun(@times,Neu_normals(:,3)',wD'*gz(xn,yn,zn));  %d2 x Nneu

end
