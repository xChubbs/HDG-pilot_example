function Mass=MassMatrix(T,coeffs,k,formula)

%{M1,M2,...}=MassMatrix(T,{c1,c2,...},k,formula)
%Input:
%           T: expanded tetrahedrization
% {c1,c2,...}: cell array with vectorized functions of three variables
%           k: polynomial degree
%     formula: quadrature formula in 3d (N x 5 matrix)
%
%Output:
% {M1,M2,...}: each cell is d3 x d3 x Nelts (\int_K c{l} P_i^K P_j^K )           
%
%Last modified: March 14, 2013

% Recovery of defined elements
Nnodes=size(formula,1);
Nelts=size(T.elements,1);
d3=nchoosek(k+3,3);

% Separation of coordinates & Multiplication with weighs of coordinate
x=T.coordinates(:,1);x=formula(:,1:4)*x(T.elements');
y=T.coordinates(:,2);y=formula(:,1:4)*y(T.elements');
z=T.coordinates(:,3);z=formula(:,1:4)*z(T.elements'); % Nnd x Nelts

xhat=formula(:,2);
yhat=formula(:,3);
zhat=formula(:,4);

% Evaluation on the reference element
P=dubiner3d(2*xhat-1,2*yhat-1,2*zhat-1,k);  % Nnd x d3

% Mass Matrix Definition
nMass=size(coeffs,2);   % Number of coefficients created for Mass matrix
Mass=cell(1,nMass);     % Storage reservation for Mass matrix (Cells)
for n=1:nMass
    c=coeffs{n};        % Coeficient definition for cell
    C=bsxfun(@times,T.volume',c(x,y,z));    % Multiplication: Nnd x Nelts
    mass=zeros(d3,Nelts*d3);    % mass matrix storage reservation
    for q=1:Nnodes
        % Definition of mass matrix:
        % It's used a recursive definition of the mass matrix per nodes.
        % the product made it's T.volumes * coeficients with
        % weights times P evaluation of P times P
        % This is the definition made at (4.6):
        %               M + <C x ((w * P)' * P)>
        mass=mass+kron(C(q,:),formula(q,5)*P(q,:)'*P(q,:));
    end
    mass=reshape(mass,[d3,d3,Nelts]);   % Here it's the shape of type (a)
    Mass{n}=mass;                       % Mass matrix for element n
end
return