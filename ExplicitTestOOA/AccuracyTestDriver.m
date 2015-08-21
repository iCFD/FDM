%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Solving 1-D wave equation with Finite Difference Methods
%
%                 dq/dt + df/dx = 0,  for x \in [a,b]
%                   where f = u*q :: linear flux
%
%              coded by Manuel Diaz, NTU, 2012.12.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: A flux reconstruction approach to high-order schemes including
% Discontinuous Galerkin methods. H.T. Huynh, AIAA 2007.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% Fixed Parameters
tEnd = 2; % One cycle for every test
IC = 2; % sinusoidal function

% Parameters
mth = [1,2,3,4,5]; % methods: {1}Upwind,{2}TVD,{3}WENO3,{4}WENO5,{5}WENO7.
cfl = [0.95,0.95,0.8,0.8,0.1]; % one CFL for every method!
nc  = [20,40,80,160,320]; % number of cells to use in every test.

% Number of parameters
p1 = length(mth);
p2 = length(nc);

% Allocate space for results
table = zeros(p2,2,p1,3);
Norm = zeros(size(table));
OOA = zeros(size(table));

%% Compute L1 and L\infty norms

for l = 1:p1
    for n = 1:p2
        tic
        [Norm(n,1,l),Norm(n,2,l)] = ...
            TestFDMfun(cfl(l),tEnd,IC,nc(n),mth(l));
        toc
    end
end

%% Compute the Order of Accuracy (OOA)

for l = 1:p1
    for n = 2:p2
        OOA(n,1,l) = log(Norm(n-1,1,l)/Norm(n,1,l))/log(2);
        OOA(n,2,l) = log(Norm(n-1,2,l)/Norm(n,2,l))/log(2);
    end
end

%% Plot figure with results
loglog(nc,Norm(:,:,1),'-s',...
    nc,Norm(:,:,2),'-o',...
    nc,Norm(:,:,4),'-h',...
    nc,Norm(:,:,4),'-<',...
    nc,Norm(:,:,5),'->')

%% Display Result
for l = 1:p1
    fprintf('***************************************************************\n')
    fprintf(' Method %d\n',mth(l));
    fprintf('***************************************************************\n')
    fprintf(' nE \t L1-Norm \t Degree \t Linf-Norm \t Degree\n');
    for n = 1:p2
        fprintf('%3.0f \t %1.2e \t %2.2f \t\t %1.2e \t %2.2f \n',...
        nc(n),Norm(n,1,l),OOA(n,1,l),Norm(n,2,l),OOA(n,2,l));
    end
end
fprintf('\n');
% By observing the degree of accuaracy with respect to the CFL number, it
% is suggested that we should tune the CFL condition for each case in order
% to get the fastest, stable and most accuarate solution.

% Manuel Diaz, NTU, 2013
% End of Test