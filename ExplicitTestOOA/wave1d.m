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
clear; close all; clc;

%% Parameters
     a = +1.0;  % scalar velocity in x direction
   CFL =  0.95;  % CFL condition
  tEnd =  2.0;  % Final time
method =    2;	% {1}Upwind,{2}TVD,{3}WENO3,{4}WENO5,{5}WENO7.
    nc =   20;	% cells

%% Preprocess
% Domain discretization
FDM = FDMethods('periodic',nc);
[x,dx] = FDM.mesh1d(-1,1,nc+1);
xc = (x(2:end)+x(1:end-1))/2; 

% set IC
%q0=TestingIC(x);
q0=IC(xc,2); qe=q0;

% Time discretization
dt = CFL*dx/abs(a); tsteps = ceil(tEnd/dt); 
dt = tEnd/tsteps;   time=dt:dt:tEnd;

% set plot range
plotrange=[-1,1,min(q0)-0.1,max(q0)+0.1];
    
%% Solver Loop 

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Using a 4rd Order 5-stage SSPRK time integration
res_q = zeros(1,nc); % Runge-Kutta residual storage
     
% load initial conditions
q=q0; it=0;

for t=time
    qo=q;
    
    % Update iteration counter
    it=it+1;
       
    switch method
        case 1  % Upwind
            
            q = qo-dt*FDM.Upwindresidual1d(q,a,dx); % cfl ~ 0.5
            
        case 2 % TVD
            
            % 1st stage  % Carefull RK2 is not working correctly!
            %dF = FDM.TVDresidual1d(q,a,dt,dx,'Koren');
            %q = qo-dt*dF/2;

            % 2nd Stage
            %dF = FDM.TVDresidual1d(q,a,dt,dx,'Koren');
            %q = qo-dt*dF;

            q = qo-dt*FDM.TVDresidual1d(q,a,dt,dx,'Koren');
            
        case 3 % WENO3 % SSP-RK3-3-stages 

            % 1st stage
            dF = FDM.WENO3residual1d(q,a,dx);
            q = qo-dt*dF;

            % 2nd Stage
            dF = FDM.WENO3residual1d(q,a,dx);
            q = 0.75*qo+0.25*(q-dt*dF);

            % 3rd stage
            dF = FDM.WENO3residual1d(q,a,dx);
            q = (qo+2*(q-dt*dF))/3;
            
        case 4 % WENO5 % LE-RK4-5-stages
            
            for RKs = 1:5
                t_local = t + rk4c(RKs)*dt;
                dF = FDM.WENO5residual1d(q,a,dx);
                res_q = rk4a(RKs)*res_q + dt*dF;
                q = q - rk4b(RKs)*res_q;
            end
            
        case 5 % WENO7 % LE-RK4-5-stages
            
            for RKs = 1:5
                t_local = t + rk4c(RKs)*dt;
                dF = FDM.WENO7residual1d(q,a,dx);
                res_q = rk4a(RKs)*res_q + dt*dF;
                q = q - rk4b(RKs)*res_q;
            end
            
        otherwise
            error('method not available');
    end
    
	% plot
    plot(xc,qe,xc,q,'.'); colormap Copper; axis(plotrange);
    title(['Upwind, dx = ',num2str(dx),' time: ',num2str(t)])
    xlabel('x points'); ylabel('q(x)');
    
    if(rem(it,1)==0)
        drawnow
    end
end

%% Compute Norms

% L1 Error
L1_nodes = sum(abs(q-qe))/nc;

% L\infty Error
Linf_nodes = max(abs(q-qe));