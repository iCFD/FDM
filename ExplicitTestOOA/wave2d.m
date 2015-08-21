%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         Solving 2-D wave equation with 1st order Upwind Method
%
%            dq/dt + df/dx + dg/dy = 0,  for x,y \in [a,b;c,d]
%                     where f = u*q  and  g = v*q
%
%              coded by Manuel Diaz, NTU, 2012.12.18
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all; close all; clc;

%% Parameters
     u = +1.0;  % scalar velocity in x direction
     v = +0.5;  % scalar velocity in y direction
   CFL =  0.5;  % CFL condition
  tEnd = 0.10;  % Final time
method = 'RK33';%'fEuler'
    nx = 400;   % cells in x
    ny = 400;   % cells in y

%% Preprocess
% Domain discretization
FDM = FDMethods('non-periodic',[nx,ny]);
%FDM = FDMethods2;
a=0; b=1; c=0; d=1; [x,dx,y,dy] = FDM.mesh2d(a,b,nx,c,d,ny);

% set IC
q0=IC2d(x,y,1); %{1} 4 Quadrants, {2} Square Jump

% Time discretization
dt=CFL*min(dy,dx)/max(abs(v),abs(u)); t=0:dt:tEnd; 

% set plot range
plotrange=[a,b/dx,c,d/dy,min(min(q0)),max(max(q0))];
    
%% Solver Loop 
% load initial conditions
q=q0; it=0;

tic
for tstep=t
    qo=q;
    
    % Update iteration counter
    it=it+1;
    
    % plot
%     subplot(1,2,1); mesh(q); colormap Copper; axis(plotrange);
%     %colorbar('location','EastOutside'); 
%     title(['Upwind, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(tstep)])
%     xlabel('x points'); ylabel('y points'); zlabel('q(x,y)');
%     subplot(1,2,2); contourf(q); colormap Copper; 
%     title(['Upwind, dx = ',num2str(dx),', dy = ',num2str(dy),', time: ',num2str(tstep)])
%     xlabel('x points'); ylabel('y points');
    
    switch method
        case 'fEuler'   % Forward Euler 
            %dF = FDM.WENO5residual2d(q,u,v,dx,dy);
            %dF = FDM.TVDresidual2d(q,u,v,dt,dx,dy,'Koren');
            %dF = FDM.Upwindresidual2d(q,u,v,dx,dy);
            q = qo-dt*dF;
        case 'RK33'     % SSP-RK3-3-stages (not good for this 1st order Up)
            % 1st stage
            %dF = FDM.WENO7residual2d(q,u,v,dx,dy);
            dF = FDM.WENO5residual2d(q,u,v,dx,dy);
            %dF = FDM.TVDresidual2d(q,u,v,dt,dx,dy,'Koren');
            %dF = FDM.Upwindresidual2d(q,u,v,dx,dy);
            q = qo-dt*dF;

            % 2nd Stage
            %dF = FDM.WENO7residual2d(q,u,v,dx,dy);
            dF = FDM.WENO5residual2d(q,u,v,dx,dy);
            %dF = FDM.TVDresidual2d(q,u,v,dt,dx,dy,'Koren');
            %dF = FDM.Upwindresidual2d(q,u,v,dx,dy);
            q = 0.75*qo+0.25*(q-dt*dF);

            % 3rd stage
            %dF = FDM.WENO7residual2d(q,u,v,dx,dy);
            dF = FDM.WENO5residual2d(q,u,v,dx,dy);
            %dF = FDM.TVDresidual2d(q,u,v,dt,dx,dy,'Koren');
            %dF = FDM.Upwindresidual2d(q,u,v,dx,dy);
            q = (qo+2*(q-dt*dF))/3;
        otherwise
            error('method not available');
    end
%     if(rem(it,10)==0)
%         drawnow
%     end
end
toc