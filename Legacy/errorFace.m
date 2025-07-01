function error=errorFace(T,px, py, pz, ph,k,formula)

%error = errorFaces(T,p,ph,k,formula)
%
%Input:
%          T: expanded tetrahedrization
%          p: vectorized function of three variables
%         ph: Pk function on skeleton (d2 x Nelts)
%          k: polynomial degree
%    formula: quadrature formula in 2d (N x 4 matrix)
%
%Output:
%      error: \| p - ph \|_{h,\partial T_h}
%
%Last modified: March 14, 2013

% Common elements definition
Nfaces = size(T.faces, 1);
Nnodes = size(formula,1);

d2 = nchoosek(k+2,2);  

% Coordinates and nodes evaluation
x=T.coordinates(:,1);xt=formula(:,[1 2 3])*x(T.faces(:,[1 2 3])');
y=T.coordinates(:,2);yt=formula(:,[1 2 3])*y(T.faces(:,[1 2 3])');
z=T.coordinates(:,3);zt=formula(:,[1 2 3])*z(T.faces(:,[1 2 3])');

% Evaluation of test function per coordinate
p_x = px(xt,yt,zt);         % Nnodes x Nelts
p_y = py(xt,yt,zt);         % Nnodes x Nelts
p_z = pz(xt,yt,zt);         % Nnodes x Nelts

% Definition of Base polynomials
% - Definition of base
DB = dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k);

% - Definition of matrix transformation
x12=x(T.faces(:,2))-x(T.faces(:,1)); %x2-x1
y12=y(T.faces(:,2))-y(T.faces(:,1)); %y2-y1
z12=z(T.faces(:,2))-z(T.faces(:,1)); %z2-z1

x13=x(T.faces(:,3))-x(T.faces(:,1)); %x3-x1
y13=y(T.faces(:,3))-y(T.faces(:,1)); %y3-y1
z13=z(T.faces(:,3))-z(T.faces(:,1)); %z3-z1

% Let's remember the normals result after the transformation:
% N = (P2 - P1) x (P3 - P1), transform matrix asociated it's.
A_0 = 0.5 * [x12, y12, z12]; A_1 = 0.5 * [x13, y13, z13];

% - Dimentional array of all possible transformations per element
A = zeros([3, 2, Nfaces]);
A(:, 1, :) = A_1'; A(:, 2, :) = A_0';

% - Generalization of basis
block_dof_Nodes_d2 = @(x) ((x-1)*d2*Nnodes +1):(x)*d2*Nnodes;

db_resized = zeros([2, 2 * d2 * Nnodes]);
db_resized(1, block_dof_Nodes_d2(1), :) = reshape(DB', [1, d2 * Nnodes]);
db_resized(2, block_dof_Nodes_d2(2), :) = reshape(DB', [1, d2 * Nnodes]);

% - Repetition of the definitition per element
rep_db = repmat(db_resized, Nfaces, 1);
rep_db = reshape(rep_db, [2, Nfaces, 2*d2*Nnodes]);

% - Dimentional correction to have per each element the 
%   Transposed matrix
db_dimentional_general = permute(rep_db, [1, 3, 2]);

% - Finally we transform the basis to generate or 3D db polynomial
db_dimentional = pagemtimes(A, db_dimentional_general);

% - Also for each coordiante we introduce:
db_dim_x = [reshape(db_dimentional(1, block_dof_Nodes_d2(1), :), [d2, Nnodes, Nfaces]);...
            reshape(db_dimentional(1, block_dof_Nodes_d2(2), :), [d2, Nnodes, Nfaces])];
db_dim_y = [reshape(db_dimentional(2, block_dof_Nodes_d2(1), :), [d2, Nnodes, Nfaces]);...
            reshape(db_dimentional(2, block_dof_Nodes_d2(2), :), [d2, Nnodes, Nfaces])];
db_dim_z = [reshape(db_dimentional(3, block_dof_Nodes_d2(1), :), [d2, Nnodes, Nfaces]);...
            reshape(db_dimentional(3, block_dof_Nodes_d2(2), :), [d2, Nnodes, Nfaces])];

% - Proyection of the basis to the values of Uhat
ph_dimentional = reshape(ph, [2*d2, 1, Nfaces]);

db_ph_x = reshape(pagemtimes(pagetranspose(db_dim_x), ph_dimentional), [Nnodes, Nfaces]);
db_ph_y = reshape(pagemtimes(pagetranspose(db_dim_y), ph_dimentional), [Nnodes, Nfaces]);
db_ph_z = reshape(pagemtimes(pagetranspose(db_dim_z), ph_dimentional), [Nnodes, Nfaces]);

% Error evaluation per coordinate 
error_x = db_ph_x - p_x;
error_y = db_ph_y - p_y;
error_z = db_ph_z - p_z;

error = error_x.^2 + error_y.^2 + error_z.^2;

error=sqrt(formula(:,4)'*error*(T.area).^2);
return