function [L1_nodes,Linf_nodes,CPUtime,it] = ...
    TestFDMfun(CFL,tEnd,Ic,nc,method,eps,rk)
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
% clear all; close all; clc;

%% Parameters
     a = +1.0;  % scalar velocity in x direction
%    CFL =  0.95;  % CFL condition
%   tEnd =  2.0;  % Final time
% method =    5;	% {1}Upwind,{2}TVD,{3}WENO3,{4}WENO5,{5}WENO7.
%     nc =   20;	% cells
%     rk = 'IMEX-SSP3';

%% Preprocess
% Domain discretization
FDM = FDMethods('periodic',nc);
[x,dx] = FDM.mesh1d(-1,1,nc+1);
xc = (x(2:end)+x(1:end-1))/2; 

% Load RK coefs
RK = RKMethods(rk); stages = RK.s;

% set IC
%q0=TestingIC(x);
q0=IC(xc,Ic); qe=-exp(eps*xc-(-tEnd+xc)*eps).*sin(pi*(xc-tEnd));

% Time discretization
dt = CFL*dx/abs(a); tsteps = ceil(tEnd/dt); 
dt = tEnd/tsteps;   time=dt:dt:tEnd;

% set plot range
%plotrange=[-1,1,min(q0)-0.1,max(q0)+0.1];
    
%% Solver Loop 

% Create storage arrays
dF = zeros(1,nc,stages);
qs = zeros(1,nc,stages);

% load initial conditions
q=q0; it=0;

tic;
for t = time
    % last step value
    qo=q;
    
    % Update iteration counter
    it=it+1;
       
    % WENO with IMEX/ARK implementations
    switch method
            
        case 3 % WENO3 

            % 1st stage
            qs(:,:,1) = qo/(1-dt*eps*RK.A(1,1));
            
            % From 2nd to 6th stage
            for s=2:stages
                dF(:,:,s-1) = FDM.WENO3residual1d(qs(:,:,s-1),a,dx);
                temp = qo;
                for j=1:s-1
                    temp = temp-dt*(RK.Ahat(s,j)*dF(:,:,j)-eps*RK.A(s,j)*qs(:,:,j));
                end
                qs(:,:,s) = temp/(1-dt*eps*RK.A(s,s));
            end
            
            % Next time step
            dF(:,:,stages) = FDM.WENO3residual1d(qs(:,:,stages),a,dx);
            q = qo;
            for j=1:stages
                q = q - dt*(RK.bhat(j)*dF(:,:,j)-eps*RK.b(j)*qs(:,:,j));
            end
            
        case 4 % WENO5
            
            % 1st stage
            qs(:,:,1) = qo/(1-dt*eps*RK.A(1,1));
            
            % From 2nd to 6th stage
            for s=2:stages
                dF(:,:,s-1) = FDM.WENO5residual1d(qs(:,:,s-1),a,dx);
                temp = qo;
                for j=1:s-1
                    temp = temp-dt*(RK.Ahat(s,j)*dF(:,:,j)-eps*RK.A(s,j)*qs(:,:,j));
                end
                qs(:,:,s) = temp/(1-dt*eps*RK.A(s,s));
            end
            
            % Next time step
            dF(:,:,stages) = FDM.WENO5residual1d(qs(:,:,stages),a,dx);
            q = qo;
            for j=1:stages
                q = q - dt*(RK.bhat(j)*dF(:,:,j)-eps*RK.b(j)*qs(:,:,j));
            end
            
        case 5 % WENO7
            
            % 1st stage
            qs(:,:,1) = qo/(1-dt*eps*RK.A(1,1));
            
            % From 2nd to 6th stage
            for s=2:stages
                dF(:,:,s-1) = FDM.WENO7residual1d(qs(:,:,s-1),a,dx);
                temp = qo;
                for j=1:s-1
                    temp = temp-dt*(RK.Ahat(s,j)*dF(:,:,j)-eps*RK.A(s,j)*qs(:,:,j));
                end
                qs(:,:,s) = temp/(1-dt*eps*RK.A(s,s));
            end
            
            % Next time step
            dF(:,:,stages) = FDM.WENO7residual1d(qs(:,:,stages),a,dx);
            q = qo;
            for j=1:stages
                q = q - dt*(RK.bhat(j)*dF(:,:,j)-eps*RK.b(j)*qs(:,:,j));
            end
            
        otherwise
            error('method not available');
    end
    
% 	% plot
%     plot(xc,q0,xc,q,'.'); colormap Copper; axis(plotrange);
%     title(['Upwind, dx = ',num2str(dx),' time: ',num2str(t)])
%     xlabel('x points'); ylabel('q(x)');
%     
%     % update plot
%     if(rem(it,1)==0)
%         drawnow
%     end
end
CPUtime=toc;

%% Compute Norms

% L1 Error
L1_nodes = sum(abs(q-qe))/nc;

% L\infty Error
Linf_nodes = max(abs(q-qe));