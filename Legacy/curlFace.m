function [tauPP,...
          tauDPx, tauDPy, tauDPz,...
          nxDP_y, nxDP_z, nyDP_x, nyDP_z, nzDP_x, nzDP_y,...
          tauDD_x, tauDD_y, tauDD_z] = curlFace(T, tau, k, formula)

%[tauPP,tauDPx, tauDPy, tauDPz,nxDP_y, nxDP_z, nyDP_x, nyDP_z, nzDP_x, ...
% nzDP_y,tauDD_x, tauDD_y, tauDD_z] = curlFace(T, tau, k, formula)
%
%Input:
%             T: expanded terahedrization
%           tau: penalization parameter for HDG (Nelts x 4)
%             k: polynomial degree 
%       formula: quadrature formula in 2d (N x 4 matrix)
%Output:
%       tauPP  :   d3 x d3   x Nelts, with <tau P_i,P_j>_{\partial K}
%       tauDPx : 4*d2 x d3   x Nelts, with <tau D_i,P_j>_e, e\in E(K)
%       tauDPy : 4*d2 x d3   x Nelts, with <tau D_i,P_j>_e, e\in E(K)
%       tauDPz : 4*d2 x d3   x Nelts, with <tau D_i,P_j>_e, e\in E(K)
%        nxDP  : 4*d2 x d3   x Nelts, with <nx D_i,P_j>_e, e\in E(K)
%        nyDP  : 4*d2 x d3   x Nelts, with <ny D_i,P_j>_e, e\in E(K)
%        nzDP  : 4*d2 x d3   x Nelts, with <nz D_i,P_j>_e, e\in E(K)
%
%Last modified: March 15, 2024

% Redefinition of elements
Nelts = size(T.elements,1);
Nnodes= size(formula,1);

d2=nchoosek(k+2,2);   
d3=nchoosek(k+3,3); 

% Definition of the element area * Evaluation of Tau coeficients
TauArea= T.area(T.facebyele').*tau;     % 4 x Nelts
perm = T.perm';                         % 4 x Nelts

% Definition on the reference element 
s=formula(:,2);         % Reference "X-Axis"
t=formula(:,3);         % Reference "Y-Axis"
weights=formula(:,4);   % Weights of formulas

weights_dim = reshape(repmat(weights, Nelts, 1), [1, Nnodes, Nelts]);

O=zeros(size(s));       % Big O definition

% ----------------- Definition of shared components ----------------- %

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

% Definition of the matrix pb, where the columns represent the 
% Respective evaluation per coordinates of the nodes
pb=[pb(1:Nnodes,:),pb(Nnodes+1:2*Nnodes,:),...
    pb(2*Nnodes+1:3*Nnodes,:),pb(3*Nnodes+1:4*Nnodes,:)]; % Nnodes x 4 * d3

% Dimention of the pbweights result
pb_dimentional = reshape(repmat(pb, 1, Nelts),[Nnodes, 4*d3, Nelts]);

% Possible permutations of the face evaluated: 6 possible combinations
points2d=[s    ,t    ; ...
          t    ,s    ; ...
          1-s-t,s    ; ...
          s    ,1-s-t; ...
          t    ,1-s-t; ...
          1-s-t,t   ];

% Definition of 2D dubiner (face) evaluation.
db=dubiner2d(2*points2d(:,1)-1,2*points2d(:,2)-1,k);      % 6 * Nnodes x d2

% Dimentional definition of the polynomial basis:
% - Dimentional of the basis
db_resized = zeros([2, 6 * d2 * Nnodes]);
db_resized(1, :, :) = reshape(db', [1, 6 * d2 * Nnodes]);

% - Repetition of the definitition per element
rep_db = repmat(db_resized, Nelts, 1);
rep_db = reshape(rep_db, [2, Nelts, 6*d2*Nnodes]);

% - Dimentional correction to have per each element the 
%   Transposed matrix
db_dimentional_general = permute(rep_db, [1, 3, 2]);

% Auxiliar function for indexing of the elements
block2=@(x) (1+(x-1)*d2):(x*d2);    % Face based indexing: 2D
block3=@(x) (1+(x-1)*d3):(x*d3);    % Element based indexing: 3D

% --------------------- Computation <tau*Di,n x Pj> ------------------- %

% Storage reservation for elements
tauDPx = zeros(4*d2,d3,Nelts); % Polynomial products x
tauDPy = zeros(4*d2,d3,Nelts); % Polynomial products y 
tauDPz = zeros(4*d2,d3,Nelts); % Polynomial products z

nxDP_y =zeros(4*d2,d3,Nelts); % X-normal Y-cross
nxDP_z =zeros(4*d2,d3,Nelts); % X-normal Z-cross

nyDP_x =zeros(4*d2,d3,Nelts); % Y-normal X-cross
nyDP_z =zeros(4*d2,d3,Nelts); % Y-normal Z-cross

nzDP_x =zeros(4*d2,d3,Nelts); % Z-normal X-cross
nzDP_y =zeros(4*d2,d3,Nelts); % Z-normal Y-cross

for l=1:4

    % Normal components separation
    Nx=T.normals(:,3*(l-1)+1)';
    Ny=T.normals(:,3*(l-1)+2)';
    Nz=T.normals(:,3*(l-1)+3)';

    % For each face we must define a different transformation matrix, this
    % process is done using:
    oneface=T.faces(T.facebyele(:,l),1:3);

    % Coordinates related to the face
    x=T.coordinates(oneface(:),1);
    y=T.coordinates(oneface(:),2);
    z=T.coordinates(oneface(:),3);
    
    x12=x(Nelts+1:2*Nelts)-x(1:Nelts); %x_2-x_1
    x13=x(2*Nelts+1:end)-x(1:Nelts);   %x_3-x_1
    
    y12=y(Nelts+1:2*Nelts)-y(1:Nelts); %y_2-y_1
    y13=y(2*Nelts+1:end)-y(1:Nelts);   %y_3-y_1
    
    z12=z(Nelts+1:2*Nelts)-z(1:Nelts); %z_2-z_1
    z13=z(2*Nelts+1:end)-z(1:Nelts);   %z_3-z_1
    
    % The transformation asociated it's defined as:
    A_0 = 0.5 * [x12, y12, z12]; A_1 = 0.5 * [x13, y13, z13];

    A = zeros([3, 2, Nelts]);
    A(:, 1, :) = A_1'; A(:, 2, :) = A_0';

    % Finally the polynomial is defined as:
    db_dimentional = pagemtimes(A, db_dimentional_general);

    % - Also for each coordiante we introduce:
    db_dim_x = reshape(db_dimentional(1, :, :), [d2, 6*Nnodes, Nelts]);
    db_dim_y = reshape(db_dimentional(2, :, :), [d2, 6*Nnodes, Nelts]);
    db_dim_z = reshape(db_dimentional(3, :, :), [d2, 6*Nnodes, Nelts]);

    % Now we can define for each face the all products asociated to the
    % definition of the space for the faces selected
    % - Transposition of Nnodes x 6*d2 for all products definitions
    db_dim_x = [db_dim_x(:,1:Nnodes, :)          ;db_dim_x(:,Nnodes+1:2*Nnodes,:);...
                db_dim_x(:,2*Nnodes+1:3*Nnodes,:);db_dim_x(:,3*Nnodes+1:4*Nnodes,:);...
                db_dim_x(:,4*Nnodes+1:5*Nnodes,:);db_dim_x(:,5*Nnodes+1:6*Nnodes,:)];
    db_dim_y = [db_dim_y(:,1:Nnodes, :)          ;db_dim_y(:,Nnodes+1:2*Nnodes,:);...
                db_dim_y(:,2*Nnodes+1:3*Nnodes,:);db_dim_y(:,3*Nnodes+1:4*Nnodes,:);...
                db_dim_y(:,4*Nnodes+1:5*Nnodes,:);db_dim_y(:,5*Nnodes+1:6*Nnodes,:)];
    db_dim_z = [db_dim_z(:,1:Nnodes, :)          ;db_dim_z(:,Nnodes+1:2*Nnodes,:);...
                db_dim_z(:,2*Nnodes+1:3*Nnodes,:);db_dim_z(:,3*Nnodes+1:4*Nnodes,:);...
                db_dim_z(:,4*Nnodes+1:5*Nnodes,:);db_dim_z(:,5*Nnodes+1:6*Nnodes,:)];

    % Cuadrature weights multiplication
    db_weights_dim_x = bsxfun(@times,weights_dim,db_dim_x);
    db_weights_dim_y = bsxfun(@times,weights_dim,db_dim_y);
    db_weights_dim_z = bsxfun(@times,weights_dim,db_dim_z);

    % Definition of all products:
    all_products_DP_x = pagemtimes(db_weights_dim_x, pb_dimentional);
    all_products_DP_y = pagemtimes(db_weights_dim_y, pb_dimentional);
    all_products_DP_z = pagemtimes(db_weights_dim_z, pb_dimentional);
    
    % The definition of corresponding permutations it's made in logical
    % expressions, it's verified that the current permutation it's
    % corresponding for each element before the kronigger product it's made
    for mu=1:6

        % ----------------- Permutation verification ----------------- %
        % Tau components
        taumu=TauArea(l,:).*(perm(l,:)==mu);

        % Normals 
        nxmu = Nx.*(perm(l,:)==mu);
        nymu = Ny.*(perm(l,:)==mu);
        nzmu = Nz.*(perm(l,:)==mu);

        % To correlate the dimentions involved we must reshape this
        % maps to agree the Nelts dimentional product
        taumu = reshape(taumu, [1, 1, Nelts]);

        nxmu = reshape(nxmu, [1, 1, Nelts]);
        nymu = reshape(nymu, [1, 1, Nelts]);
        nzmu = reshape(nzmu, [1, 1, Nelts]);
        
        % --------------------- Tau x-coordinate --------------------- % 
        % Definition by coordiantes of Tau, for the first face we have:
        tauDPx(block2(l),:, :) = tauDPx(block2(l),:, :) +...
        bsxfun(@times, taumu, all_products_DP_x(block2(mu),block3(l),:));

        % --------------------- Tau y-coordinate --------------------- % 
        % Definition by coordiantes of Tau, for the first face we have:
        tauDPy(block2(l),:, :) = tauDPy(block2(l),:, :) +...
        bsxfun(@times, taumu, all_products_DP_y(block2(mu),block3(l),:));

        % --------------------- Tau z-coordinate --------------------- % 
        % Definition by coordiantes of Tau, for the first face we have:
        tauDPz(block2(l),:, :) = tauDPz(block2(l),:, :) +...
        bsxfun(@times, taumu, all_products_DP_z(block2(mu),block3(l),:));
        
        % ------------------------- n x-normal ----------------------- % 
        % The nx on the block it's redifined if the 
        % correct permutation it's present
        nxDP_y(block2(l),:, :) = nxDP_y(block2(l),:, :) + ...
         bsxfun(@times, nxmu, all_products_DP_y(block2(mu),block3(l),:));

        nxDP_z(block2(l),:, :) = nxDP_z(block2(l),:, :) + ...
         bsxfun(@times, nxmu, all_products_DP_z(block2(mu),block3(l),:));
        
        % ------------------------ n y-normal ------------------------ %
        % The ny on the block it's redifined if the 
        % correct permutation it's present
        nyDP_x(block2(l),:, :) = nyDP_x(block2(l),:, :) + ...
         bsxfun(@times, nymu, all_products_DP_x(block2(mu),block3(l),:));

        nyDP_z(block2(l),:, :) = nyDP_z(block2(l),:, :) + ...
         bsxfun(@times, nymu, all_products_DP_z(block2(mu),block3(l),:));
        
        % ------------------------ n z-normal ------------------------ %
        % The nz on the block it's redifined if the 
        % correct permutation it's present
        nzDP_x(block2(l),:, :) = nzDP_x(block2(l),:, :) + ...
         bsxfun(@times, nzmu, all_products_DP_x(block2(mu),block3(l),:));

        nzDP_y(block2(l),:, :) = nzDP_y(block2(l),:, :) + ...
         bsxfun(@times, nzmu, all_products_DP_y(block2(mu),block3(l),:));

    end
end

% Following the representation of the bilinear form, we must mantain the
% form of A_2, for this is proposed the product: In relation to mantain 
% cross product relations:
%
%       |    O     DP_K^T -DP_K^T   O   |   |  n_x |
% A_2 = |  DP_K^T     O    -DP_K^T  O   | * |  n_y |
%       | -DP_K^T  DP_K^T     O     O   |   |  n_z |
%       |    O        O       O    tau  |   | -tau | 
%
% Here the result serves as the initial vector proposed, but the matrix
% relation mantains the cross product relation between every coordinate.
% It's important to ask, if the final result can be aproximated, in this
% case assumming I can compute the product, I can produce the final as:
% 
%       | (n_y * DP_K^T - n_z * DP_K^T) |
% A_2 = | (n_x * DP_K^T - n_z * DP_K^T) |
%       | (n_y * DP_K^T - n_x * DP_K^T) |
%       |        - tau * DP_K^T         | 
%
% Then the aproximated polynomial it's the same as the value proposed 
% for the product of both polynomials. We know each of this differences
% mainly because we have the products: niDP

% --------------------- Computation <tau*Pi,Pj> ------------------- %

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
pbpb1 = pbweights(1:Nnodes,:)'            * pb(1:Nnodes,:);
pbpb2 = pbweights(Nnodes+1:2*Nnodes,:)'   * pb(Nnodes+1:2*Nnodes,:);
pbpb3 = pbweights(2*Nnodes+1:3*Nnodes,:)' * pb(2*Nnodes+1:3*Nnodes,:);
pbpb4 = pbweights(3*Nnodes+1:4*Nnodes,:)' * pb(3*Nnodes+1:4*Nnodes,:);

% The operation takes in count the separation before, using the sum of
% the separated integrals.
tauPP=kron(TauArea(1,:),pbpb1)+kron(TauArea(2,:),pbpb2)...
      +kron(TauArea(3,:),pbpb3)+kron(TauArea(4,:),pbpb4);

% Reshape on a type (a) matrix
tauPP=reshape(tauPP,[d3,d3,Nelts]);



% ------------------------- Computation tauDD -------------------------- %

% Dubiner 2D evaluation on nodes (No permutation needed)
db=dubiner2d(2*s-1,2*t-1,k);

% Dimentional definition of the polynomial basis:
% - Dimentional of the basis
db_resized = zeros([2, d2 * Nnodes]);
db_resized(1, :, :) = reshape(db', [1, d2 * Nnodes]);

% - Repetition of the definitition per element
rep_db = repmat(db_resized, Nelts, 1);
rep_db = reshape(rep_db, [2, Nelts, d2*Nnodes]);

% - Dimentional correction to have per each element the 
%   Transposed matrix
db_dimentional_general = permute(rep_db, [1, 3, 2]);

tauDD_x=zeros(4*d2,4*d2,Nelts); % Storage reservation
tauDD_y=zeros(4*d2,4*d2,Nelts); % Storage reservation
tauDD_z=zeros(4*d2,4*d2,Nelts); % Storage reservation

block2=@(x) (1+(x-1)*d2):(x*d2);    % Face based indexing: 2D
TauArea= T.area(T.facebyele').*tau; % 4 x Nelts

for l=1:4

    % For each face we must define a different transformation matrix, this
    % process is done using:
    oneface=T.faces(T.facebyele(:,l),1:3);

    % Coordinates related to the face
    x=T.coordinates(oneface(:),1);
    y=T.coordinates(oneface(:),2);
    z=T.coordinates(oneface(:),3);
    
    x12=x(Nelts+1:2*Nelts)-x(1:Nelts); %x_2-x_1
    x13=x(2*Nelts+1:end)-x(1:Nelts);   %x_3-x_1
    
    y12=y(Nelts+1:2*Nelts)-y(1:Nelts); %y_2-y_1
    y13=y(2*Nelts+1:end)-y(1:Nelts);   %y_3-y_1
    
    z12=z(Nelts+1:2*Nelts)-z(1:Nelts); %z_2-z_1
    z13=z(2*Nelts+1:end)-z(1:Nelts);   %z_3-z_1
    
    % The transformation asociated it's defined as:
    A_0 = 0.5 *[x12, y12, z12]; A_1 = 0.5 *[x13, y13, z13];

    A = zeros([3, 2, Nelts]);
    A(:, 1, :) = A_1'; A(:, 2, :) = A_0';

    % Finally the polynomial is defined as:
    db_dimentional = pagemtimes(A, db_dimentional_general);

    % - Also for each coordiante we introduce:
    db_dim_x = reshape(db_dimentional(1, :, :), [d2, Nnodes, Nelts]);
    db_dim_y = reshape(db_dimentional(2, :, :), [d2, Nnodes, Nelts]);
    db_dim_z = reshape(db_dimentional(3, :, :), [d2, Nnodes, Nelts]);

    dweights_x=bsxfun(@times,weights', db_dim_x);
    dweights_y=bsxfun(@times,weights', db_dim_y); 
    dweights_z=bsxfun(@times,weights', db_dim_z);  

    dwd_x = pagemtimes(dweights_x, pagetranspose(db_dim_x));
    dwd_y = pagemtimes(dweights_y, pagetranspose(db_dim_y));
    dwd_z = pagemtimes(dweights_z, pagetranspose(db_dim_z));

    TauAreaFace = reshape(TauArea(l,:), [1, 1, Nelts]);    % Reshape of TauArea

    % Tau DD takes place on the 2D blocks: and it's reshape
    tauDD_x(block2(l),block2(l),:)= tauDD_x(block2(l),block2(l),:) + ...
                                         bsxfun(@times, TauAreaFace,dwd_x);
    tauDD_y(block2(l),block2(l),:)= tauDD_y(block2(l),block2(l),:) + ...
                                         bsxfun(@times, TauAreaFace,dwd_y);
    tauDD_z(block2(l),block2(l),:)= tauDD_z(block2(l),block2(l),:) + ...
                                         bsxfun(@times, TauAreaFace,dwd_z);
end

% This is the type (b) matrix, corresponding to the shape 4*d2 x d2 x Nelts

end

