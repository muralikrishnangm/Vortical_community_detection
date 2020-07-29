function vortical_community
% Vortical community detection in 2D vortical flows to reduce the dimension
% of the flow field into community centroids

% Flow field :      Airfoil with Gurney flap 
%                   (doi.org/10.2514/1.J056260)

% Date Created:     08/14/2018
% Date Modified:    01/21/2020; Louvain & clean-up
% Author:           Muralikrishnan Gopalakrishnan Meena, Florida State Unv.
%                   mg15h@my.fsu.edu
% (compiled using MATLAB 2018a)

% NOTE: Please use the code after verification. The code is not created
% for computational performance. MGM provides no guarantees for this code.
% Use as-is and for academic research use only, no commercial use allowed
% without permission. For citations, please use the reference below:

% Reference:    M. Gopalakrishnan Meena, A. G. Nair & K. Taira
%               "Network Community-Based Model Reduction for Vortical
%               Flows", Physical Review E, 97, 063103, 2018
%               (doi.org/10.1103/PhysRevE.97.063103)

clc
%% Load data & visualize
% xx,yy: x & y grid;    omg: vorticity; 
% xfil,yfil,omgfil: vorticity thresholded data;
% Gamma_v: circulation (omgfil*dx*dy)
load('example_airfoil.mat','xx','yy','omg','xfil','yfil','Gamma_v')  
figure(1), clf, hold all, axis off
ax1 = axes;
contour(xx,yy,omg,linspace(-20,20),'linewidth',2)
axis off, axis equal, axis([-0.5 7 -2 2]), colormap(ax1,parula)

%% Compute full adjacency matrix 
% (directed graph, edges=algebraic mean of induced velocities)
A = adjacency_mat(xfil,yfil,Gamma_v,1,'alg');

%% Identify communities in the vortical network
% (Split A into +ve & -ve vortical elements & then find communities)
[Ci,numcom] = find_vortical_communities(Gamma_v,A,1.1);

%% Find community centroids
xc = NaN(numcom,1); yc = NaN(numcom,1); Gamma_c = NaN(numcom,1);
for i = 1:numcom
    idc = Ci==i;
    xc(i) = sum(xfil(idc).*Gamma_v(idc)) / sum(Gamma_v(idc));
    yc(i) = sum(yfil(idc).*Gamma_v(idc)) / sum(Gamma_v(idc));
    Gamma_c(i) = sum(Gamma_v(idc));
end

%% Visualize dimensionally reduced flow field 
figure(1)
ax2 = axes; hold all;
scatter(xfil,yfil,10,Ci,'filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2)
plot(xc,yc,'k*','LineWidth',2)
axis off, axis equal, axis([-0.5 7 -2 2]), colormap(ax2,jet)


%% ***********************Functions**************************************
function [A] = adjacency_mat(x,y,Gamma,alpha,type)
%%
% Adjacency matrix calculator (2D fluid flow)
% Date Created  : 05/20/2016
% Date Modified : 08/24/2017, alpha, vectorized, type
% By MGM
% Reference:    A. G. Nair & K. Taira
%               "Network-Theoretic Approach to Sparsified Discrete Vortex
%               Dynamics", Journal of Fluid Mechanics, 768, 549-571, 2015
%               (doi.org/10.1017/jfm.2015.97)
%%
% [A] = adjacency_mat(x,y,gamma,alpha,type)
% Function to calculate adjacency matrix given x, y and gamma. 
% A         : adjacency matrix using algebraic mean (default)
% x         : x coordinate vector
% y         : y coordinate vector
% gamma     : circulation vector
% alpha     : network direction parameter (0.5=undirected)
% type      : algebraic or geometric mean ('alg' [default] or 'goem')

if isempty(type)
    type = 'alg';
end
pts = length(x);
x = reshape(x,pts,1); y = reshape(y,pts,1);
r(1:pts,1:2) = [x y];               % position vector matrix, r = (xi,yj)
Gamma = Gamma/(2*pi);               % scalling for Gamma in 2D flow

fprintf('\nCalculating A matrix of %i nodes...\n',length(x));
A = NaN(pts,pts);
for i = 1:pts   % can use parfor if needed
    dist_ij = r(1:pts,:) - ones(pts,1)*r(i,:);  % distance vector, r_{i->j}
    dist_mag = sqrt(sum(dist_ij.^2,2));         % ||r_{i->j}||
    dist_ij = 1./dist_mag;                      % distance ratio
    dist_ij(isinf(dist_ij)) = 0;                % update r_{i->i} = 0
    u_ij = abs(Gamma(i)) * dist_ij;             % u_{i->j}
    u_ji = abs(Gamma) .* dist_ij;               % u_{j->i}
    % A_{i->j}
    if strcmp(type,'alg') == 1
        A(i,:) = alpha*u_ij + (1-alpha)*u_ji;
    elseif strcmp(type,'geom') == 1
        A(i,:) = alpha*u_ij .* (1-alpha)*u_ji;
    else
        error('Eneter valid network definition "type"')
    end
end

function [Ci_full,numcom] = find_vortical_communities(omgfil,Afull,gamma)
%%
% Find communities in a 2D vortical flow network
% Date Created:     08/14/2018
% Date Modified:    01/22/2020; function form
% By MGM
%%
% [Ci_full,numcom] = find_vortical_communities(omgfil,Afull,gamma)
% Function to identify communities in a 2D vortical flow network.
% Communities of +ve & -ve vorticity are identified from the adjacency
% matrix of the full system
% Ci            : community index of elements
% numcom        : number of communities
% omgfil        : vorticity of the elements
% Afull         : Adjacency matrix of the elements
% gamma         : resolution parameter for community detection

% nodes of +ve & -ve vorticity
id_p = find(omgfil(:)>0);
id_n = find(omgfil(:)<0);

% +ve vorticity elements
fprintf('Finding communities of +ve rotational vortical elements...\n');
A_p = Afull(id_p,id_p);
Ci_in = 1:numel(id_p);     % IC: all nodes are in different communities
[~,Ci_p] = find_communities(gamma,Ci_in,A_p);

% -ve vorticity elements
fprintf('Finding communities of -ve rotational vortical elements...\n');
A_n = Afull(id_n,id_n);
Ci_in = 1:numel(id_n);     % IC: all nodes are in different communities
[~,Ci_n] = find_communities(gamma,Ci_in,A_n);

% combine community info of both +ve & -ve vorticity elements
Ci_full = NaN(numel(omgfil),1);
Ci_full(id_p) = Ci_p; Ci_full(id_n) = Ci_n + numel(unique(Ci_p));
numcom = numel(unique(Ci_full));


function [Q,Ci] = find_communities(gamma,Ci_in,A)
% Function to find communities in vortical (or any) networks

% Using Modularity max.
% [Ci, Q] = modularity_dir(A,gamma);

% Using Louvain algorithm
Ci  = Ci_in;                    % initial community affiliations
Q0 = -1; Q = 0;                 % initialize modularity values
while Q-Q0>1e-5                 % while modularity increases
    Q0 = Q;                     % perform community detection
    [Ci, Q] = community_louvain(A,gamma,Ci);
end


