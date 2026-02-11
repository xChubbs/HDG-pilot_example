function error=errorFace(T, px, py, pz, ph,k,formula)

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
faces = T.faces;

Nfaces = size(faces, 1);

d2 = nchoosek(k+2,2);

% Coordinates and nodes evaluation
x=T.coordinates(:,1);xt=formula(:,[1 2 3])*x(faces(:,[1 2 3])');
y=T.coordinates(:,2);yt=formula(:,[1 2 3])*y(faces(:,[1 2 3])');
z=T.coordinates(:,3);zt=formula(:,[1 2 3])*z(faces(:,[1 2 3])');

% Evaluation of test function per coordinate
p_x = px(xt,yt,zt);         % Nnodes x Nelts
p_y = py(xt,yt,zt);         % Nnodes x Nelts
p_z = pz(xt,yt,zt);         % Nnodes x Nelts

% Definition of Base polynomials
% - Definition of base
DB = dubiner2d(2*formula(:,2)-1,2*formula(:,3)-1,k);

% - Extension of base for product:
DB = repmat(DB, 1, 2);

A_0 = reshape(T.A(2, :, :), [3, Nfaces]);
A_1 = reshape(T.A(1, :, :), [3, Nfaces]);

% Separation of results by coordinate
ph_x = [A_1(1, :) .* ph(1:d2, :); A_0(1, :) .* ph(d2+1:end, :)];
ph_y = [A_1(2, :) .* ph(1:d2, :); A_0(2, :) .* ph(d2+1:end, :)];
ph_z = [A_1(3, :) .* ph(1:d2, :); A_0(3, :) .* ph(d2+1:end, :)];

% Products per coodinate
db_ph_x = DB * ph_x;
db_ph_y = DB * ph_y;
db_ph_z = DB * ph_z;

% Error evaluation per coordinate 
error_x = db_ph_x - p_x;
error_y = db_ph_y - p_y;
error_z = db_ph_z - p_z;

error = error_x + error_y + error_z;

error= formula(:,4)'*error*T.area.^2;
return