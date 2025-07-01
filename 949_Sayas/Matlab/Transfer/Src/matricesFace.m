function [tauPP,tauDP,nxDP,nyDP,nzDP,tauDD]=matricesFace(T,tau,k,formula)

%[tauPP,tauDP,nxDP,nyDP,nzDP,tauDD]=matricesFace(T,tau,k,formula)
%
%Input:
%             T: expanded terahedrization
%           tau: penalization parameter for HDG (Nelts x 4)
%             k: polynomial degree 
%       formula: quadrature formula in 2d (N x 4 matrix)
%Output:
%       tauPP :   d3 x d3   x Nelts, with <tau P_i,P_j>_{\partial K}
%       tauDP : 4*d2 x d3   x Nelts, with <tau D_i,P_j>_e, e\in E(K)
%        nxDP : 4*d2 x d3   x Nelts, with <nx D_i,P_j>_e, e\in E(K)
%        nyDP : 4*d2 x d3   x Nelts, with <ny D_i,P_j>_e, e\in E(K)
%        nzDP : 4*d2 x d3   x Nelts, with <nz D_i,P_j>_e, e\in E(K)
%       tauDD : 4*d2 x 4*d2 x Nelts, block diag <tau D_i, D_j>_e, e\in E(K)
%
%Last modified: March 14, 2013

% Redefinition of elements
Nelts = size(T.elements,1);
Nnodes= size(formula,1);

d2=nchoosek(k+2,2);   
d3=nchoosek(k+3,3); 

% Definition of the element area * Evaluation of Tau coeficients
TauArea= T.area(T.facebyele').*tau;     % 4 x Nelts
T.perm = T.perm';                       % 4 x Nelts

% Definition on the reference element 
s=formula(:,2);         % Reference "X-Axis"
t=formula(:,3);         % Reference "Y-Axis"
weights=formula(:,4);   % Weights of formulas

% --------------------- Computation <tau*Pi,Pj> ------------------- %

O=zeros(size(s));   % Big O definition

% Groups of quadrature points on the face K_hat
% This will need the use of the permutation matrix to conserve
% the matching orientations for these points
points3d=[s,t,O;...
          s,O,t;...
          O,s,t;...
          s,t,1-s-t];

% Evaluation of the quadrature points on face K_hat
pb=dubiner3d(2*points3d(:,1)-1, ...  % Evaluation on quadrature points
             2*points3d(:,2)-1, ...  % of the Dubiner basis.
             2*points3d(:,3)-1, k);  % (4 * Nnodes) x d3

pbweights=bsxfun(@times,[formula(:,4);...   % Weights definition
                         formula(:,4);...   % on the basis defined.
                         formula(:,4);...
                         formula(:,4)],pb);

% Indexing based evaluation: Definition of pbweights as a column vector
% Take in count the evaluation it's respectibly:
% It's operated k * Nnodes : (k + 1)(Nnodes + 1) with k = 0 : 3 
pbpb1=pbweights(1:Nnodes,:)'*pb(1:Nnodes,:);
pbpb2=pbweights(Nnodes+1:2*Nnodes,:)'*pb(Nnodes+1:2*Nnodes,:);
pbpb3=pbweights(2*Nnodes+1:3*Nnodes,:)'*pb(2*Nnodes+1:3*Nnodes,:);
pbpb4=pbweights(3*Nnodes+1:4*Nnodes,:)'*pb(3*Nnodes+1:4*Nnodes,:);

% The operation takes in count the separation before, using the sum of
% the separated integrals.
tauPP=kron(TauArea(1,:),pbpb1)+kron(TauArea(2,:),pbpb2)...
      +kron(TauArea(3,:),pbpb3)+kron(TauArea(4,:),pbpb4);

% Reshape on a type (a) matrix
tauPP=reshape(tauPP,[d3,d3,Nelts]);

% -------------- Computation <alpha*D,P>, alpha=tau,nx,ny,nz ------------- %

% This is the type (c) integral, the mixed terms evaluation.

% Definition of the matrix pb, where the columns represent the 
% Respective evaluation per coordinates of the nodes
pb=[pb(1:Nnodes,:),pb(Nnodes+1:2*Nnodes,:),...
    pb(2*Nnodes+1:3*Nnodes,:),pb(3*Nnodes+1:4*Nnodes,:)]; % Nnodes x 4 * d3

% Possible permutations of the face evaluated: 6 possible combinations
points2d=[s,t;...
          t,s;...
          1-s-t,s;...
          s,1-s-t;...
          t,1-s-t;...
          1-s-t,t];

% Definition of 2D dubiner (face) evaluation.
db=dubiner2d(2*points2d(:,1)-1,2*points2d(:,2)-1,k);      % 6 * Nnodes x d2

% Same process made with the 3D dubiner, and the matrix it's defined for
% all possible permutations
db=[db(1:Nnodes,:),db(Nnodes+1:2*Nnodes,:),...
    db(2*Nnodes+1:3*Nnodes,:),db(3*Nnodes+1:4*Nnodes,:),...
    db(4*Nnodes+1:5*Nnodes,:),db(5*Nnodes+1:6*Nnodes,:)]; % Nnodes x 6 * d2

db=bsxfun(@times,weights,db);   % weights formulas times dubiner eval
allproducts=db'*pb;                                       % 6 * d2 x 4 * d3

% Auxiliar function for indexing of the elements
block2=@(x) (1+(x-1)*d2):(x*d2);    % Face based indexing: 2D
block3=@(x) (1+(x-1)*d3):(x*d3);    % Element based indexing: 3D

% Storage reservation for elements
tauDP=zeros(4*d2,d3*Nelts); % Polynomial products

nxDP =zeros(4*d2,d3*Nelts); % X-normal components
nyDP =zeros(4*d2,d3*Nelts); % Y-normal components
nzDP =zeros(4*d2,d3*Nelts); % Z-normal components

for l=1:4

    % Normal components separation
    Nx=T.normals(:,3*(l-1)+1)';
    Ny=T.normals(:,3*(l-1)+2)';
    Nz=T.normals(:,3*(l-1)+3)';

    % The definition of corresponding permutations it's made in logical
    % expressions, it's verified that the current permutation it's
    % corresponding for each element before the kronigger product it's made
    for mu=1:6
        
        % ---------------------- Tau ------------------- % 
        % Permutation verification
        taumu=TauArea(l,:).*(T.perm(l,:)==mu);

        % The Tau DP on the block it's redifined if
        % the correct permutation it's present
        tauDP(block2(l),:)=tauDP(block2(l),:)+...
            kron(taumu,allproducts(block2(mu),block3(l)));
        
        % ------------------ n x-normal --------------- % 
        % Permutation verification
        nxmu=Nx.*(T.perm(l,:)==mu);

        % The nx on the block it's redifined if the 
        % correct permutation it's present
        nxDP(block2(l),:)=nxDP(block2(l),:)+...
            kron(nxmu,allproducts(block2(mu),block3(l)));
        
        % ------------------ n y-normal --------------- %
        % Permutation verification
        nymu=Ny.*(T.perm(l,:)==mu);

        % The ny on the block it's redifined if the 
        % correct permutation it's present
        nyDP(block2(l),:)=nyDP(block2(l),:)+...
            kron(nymu,allproducts(block2(mu),block3(l)));
        
        % ------------------ n z-normal --------------- %
        % Permutation verification
        nzmu=Nz.*(T.perm(l,:)==mu);

        % The nz on the block it's redifined if the 
        % correct permutation it's present
        nzDP(block2(l),:)=nzDP(block2(l),:)+...
            kron(nzmu,allproducts(block2(mu),block3(l)));
    end
end

% Reshape of elements
tauDP=reshape(tauDP,[4*d2,d3,Nelts]);   % Tau elements 4 * d2 x d3 x Nelts
nxDP=reshape(nxDP,[4*d2,d3,Nelts]);     % Normal X 4 * d2 x d3 x Nelts
nyDP=reshape(nyDP,[4*d2,d3,Nelts]);     % Normal Y 4 * d2 x d3 x Nelts
nzDP=reshape(nzDP,[4*d2,d3,Nelts]);     % Normal Z 4 * d2 x d3 x Nelts

% Take in count this is the type (c) matrix shape.

% ------------------------- Computation tauDD -------------------------- %

% Dubiner 2D evaluation on nodes (No permutation needed)
d=dubiner2d(2*s-1,2*t-1,k);

dweights=bsxfun(@times,d,weights);  % Dubiner times weights of formulas

dwd=dweights'*d; % This is the DD computation then it's made the d

tauDD=zeros(4*d2,4*d2,Nelts); % Storage reservation

for l=1:4
    % Tau DD takes place on the 2D blocks: and it's reshape
    tauDD(block2(l),block2(l),:)=reshape(kron(TauArea(l,:),dwd),...
                                         [d2,d2,Nelts]);
end

% This is the type (b) matrix, corresponding to the shape 4*d2 x d2 x Nelts
return