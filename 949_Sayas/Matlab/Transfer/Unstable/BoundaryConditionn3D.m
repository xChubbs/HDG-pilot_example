function [uhd,qhn]=BoundaryConditionn3D(ux, uy, uz, gx, gy, gz, T, k, formula)

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

% Common elements definition
d2     = nchoosek(k+2,2);    
Ndir   = size(T.dirfaces, 2);
Nneu   = size(T.neufaces, 2);
Nnodes = size(formula(:,1), 1);

Area = T.area(T.dirfaces)';

% Separation on components of coordinates
x=T.coordinates(:,1);
y=T.coordinates(:,2);
z=T.coordinates(:,3);

% Dirichlet elements
xx=formula(:,[1 2 3])*x(T.dirichlet');  % Nnodes x Ndir
yy=formula(:,[1 2 3])*y(T.dirichlet');  % Nnodes x Ndir
zz=formula(:,[1 2 3])*z(T.dirichlet');  % Nnodes x Ndir

% Dirichlet condition per coordinate
g_x = ux(xx,yy,zz); g_y = uy(xx,yy,zz); g_z = uz(xx,yy,zz); 

% Dubiner evaluation on reference element
D=dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k);  % Nnodes x d2

% - Weights times Dubiner eval on reference element
wD=bsxfun(@times,formula(:,4),D);                  % Nnodes x d2

% ------------------------ Dirichlet BC ------------------------ %

% Dirichlet elements
x12=x(T.dirichlet(:,2))-x(T.dirichlet(:,1)); %x2-x1
y12=y(T.dirichlet(:,2))-y(T.dirichlet(:,1)); %y2-y1
z12=z(T.dirichlet(:,2))-z(T.dirichlet(:,1)); %z2-z1

x13=x(T.dirichlet(:,3))-x(T.dirichlet(:,1)); %x3-x1
y13=y(T.dirichlet(:,3))-y(T.dirichlet(:,1)); %y3-y1
z13=z(T.dirichlet(:,3))-z(T.dirichlet(:,1)); %z3-z1

% Dirichlet normals: Decomposition of components
Dir_normals= [y12.*z13-z12.*y13,...
              z12.*x13-x12.*z13,...
              x12.*y13-x13.*y12];    % Ndir x 3,  0.5*||n^e|| = |e|

Dir_normals = Dir_normals' ./ vecnorm(Dir_normals)';

Nx = reshape(Dir_normals(1, :), [1, 1, Ndir]);       % Dirichlet Normals x-component
Ny = reshape(Dir_normals(2, :), [1, 1, Ndir]);       % Dirichlet Normals y-component      
Nz = reshape(Dir_normals(3, :), [1, 1, Ndir]);       % Dirichlet Normals z-component

% Definition of the Nx matrix
O = zeros(1, 1, Ndir);  % Zero vector definition

Nx = [  O,  Nz, -Ny; ...
      -Nz,   O,  Nx; ...
       Ny, -Nx,  O];

% ------- Redefinition of space: dirichlet conditions ------- %

% Let's remember the normals result after the transformation:
% N = (P2 - P1) x (P3 - P1), transform matrix asociated it's.
A_0 = reshape(T.A(1, :, T.dirfaces), [3, 1, Ndir]);
A_1 = reshape(T.A(2, :, T.dirfaces), [3, 1, Ndir]);

% ------------- Computation of < n x u_hat, \eta >_{\partial \omega} -------------- %
n_cross_A0 = pagemtimes(Nx, A_0); n_cross_A1 = pagemtimes(Nx, A_1);

% Transformation matrix definition:
A = [pagemtimes(pagetranspose(n_cross_A1), A_1), pagemtimes(pagetranspose(n_cross_A0), A_1); ...
     pagemtimes(pagetranspose(n_cross_A1), A_0), pagemtimes(pagetranspose(n_cross_A0), A_0)];

% Construction of matrix < n x u_hat, \eta >_{\partial \omega}
uhat_test = arrayfun(@(i) kron(Area(i) * A(:, :, i), wD'*D), ...
                                       1:Ndir, 'UniformOutput',false);

uhat_test = cat(3,uhat_test{:});

% ---------------- Computation of < n x g, \eta >_{\partial \omega} ---------------- %

% Definition of d2 evaluation per coordinate
wD_gx = reshape(wD' * g_x, [d2, 1, Ndir]); 
wD_gy = reshape(wD' * g_y, [d2, 1, Ndir]);
wD_gz = reshape(wD' * g_z, [d2, 1, Ndir]); 

% Definition of independent vector

% - Definition of Boundary condition product
g_test_A1 = n_cross_A1(1, :, :) .* wD_gx + ...
            n_cross_A1(2, :, :) .* wD_gy + ...
            n_cross_A1(3, :, :) .* wD_gz;

g_test_A0 = n_cross_A0(1, :, :) .* wD_gx + ...
            n_cross_A0(2, :, :) .* wD_gy + ...
            n_cross_A0(3, :, :) .* wD_gz;

g_test = reshape(Area, [1, 1, Ndir]) .* -[g_test_A1 ; g_test_A0];

% Solution of weak condition
uhd = pagemtimes(pageinv(uhat_test), g_test);       % 2*d2 x Ndir

% -------------------------- Neumman BC -------------------------- %

% Neumman elements
x12=x(T.neumann(:,2))-x(T.neumann(:,1)); %x2-x1
y12=y(T.neumann(:,2))-y(T.neumann(:,1)); %y2-y1
z12=z(T.neumann(:,2))-z(T.neumann(:,1)); %z2-z1
x13=x(T.neumann(:,3))-x(T.neumann(:,1)); %x3-x1
y13=y(T.neumann(:,3))-y(T.neumann(:,1)); %y3-y1
z13=z(T.neumann(:,3))-z(T.neumann(:,1)); %z3-z1

Neu_normals=0.5*[y12.*z13-z12.*y13,...  % Neumann normals
                 z12.*x13-x12.*z13,...
                 x12.*y13-x13.*y12];    % Nneu x 3,  ||n^e||=|e|

% Reference element coordiantes in Neumann faces times formula weights
xn=formula(:,[1 2 3])*x(T.neumann');  % Nnodes x Nneu
yn=formula(:,[1 2 3])*y(T.neumann');  % Nnodes x Nneu
zn=formula(:,[1 2 3])*z(T.neumann');  % Nnodes x Nneu

% Integral quadrature definition per components.
%qhn=bsxfun(@times,Neu_normals(:,1)',wD'*gx(xn,yn,zn))+...
%    bsxfun(@times,Neu_normals(:,2)',wD'*gy(xn,yn,zn))+...
%    bsxfun(@times,Neu_normals(:,3)',wD'*gz(xn,yn,zn));  %d2 x Nneu

qhn = zeros(d2, Nneu);

end
