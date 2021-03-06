classdef FDMethods2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Finite Difference Methods for Linear conservation equations
    %
    %   This class cotain fixed subroutines for building cartesian grids
    %   and  computing the residuals for the following methods:
    %               1st order Double sidedd Upwind,
    %               2nd order TVD,
    %               5th order WENO reconstruction.
    %               7th order WENO reconstruction.
    %
    %             coded by Manuel Diaz, NTU, 2013.12.22
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE: 
    % This is the original implementation using Matlab's circshift function.
    % Use this formulation exclusivelly when Periodic BCs are required!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
    end
        
    methods (Static)
        function [x,dx,idx] = mesh1d(a,b,nx)
            % Build a 1D grid for DFM
            % coded by Manuel Diaz, NTU, 2012.12.22
            %
            % INPUTS
            % [a,b]: x range, nx: the desided number of nodes.
            %
            % OUTPUTS
            % uniformly distributed mesh grid for computation.
            % x: vector x's coordinates, dx: cell size.
            
            % Compute cell size
            dx = (b-a)/(nx-1);
            
            % Compute vector arrays
            x = a:dx:b; idx=1:nx;
        end
        
        function [x,dx,y,dy,idxy,idx,idy] = mesh2d(a,b,nx,c,d,ny)
            % Build a 2D grid for DFM
            % coded by Manuel Diaz, NTU, 2012.12.22
            %
            % INPUTS
            % [a,b]: x range, nx: the desided number of nodes.
            % [c,d]: y range, ny: the desided number of nodes.
            %
            % OUTPUTS
            % x: vector x's coordinates, dx: cell size.
            % y: vector y's coordinates, dy: cell size.
            
            % Compute cell size
            dx = (b-a)/(nx-1); dy = (d-c)/(ny-1);
            
            % Compute vector arrays
            x = a:dx:b; y = c:dy:d; [x,y] = meshgrid(x,y);
            idx = 1:nx; idy = 1:ny; idxy  = reshape(1:nx*ny,ny,nx);
        end
        
        function [x,dx,y,dy,z,dz,idxyz,idx,idy,idz] = mesh3d(a,b,nx,c,d,ny,e,f,nz)
            % Build a 3D grid for DFM
            % coded by Manuel Diaz, NTU, 2012.12.22
            %
            % INPUTS
            % [a,b]: x range, nx: the desided number of nodes.
            % [c,d]: y range, ny: the desided number of nodes.
            % [e,f]: z range, nz: the desided number of nodes.
            %
            % OUTPUTS
            % x: vector x's coordinates, dx: cell size.
            % y: vector y's coordinates, dy: cell size.
            % z: vector z's coordinates, dz: cell size.
            
            % Compute cell size
            dx = (b-a)/(nx-1); dy = (d-c)/(ny-1); dz = (f-e)/(nz-1);
            
            % Compute vector arrays
            x = a:dx:b; y = c:dy:d; z = e:dz:f; [x,y,z] = meshgrid(x,y,z);
            idx = 1:nx; idy = 1:ny; idz = 1:nx; idxyz = reshape(1:nx*ny*nz,ny,nx,nz);
        end
        
        function dF = Upwindresidual1d(q,ax,dx)
            % Compute the Residual for 2d wave equation using Upwind Method
            %
            %           Residual: RHS of dq/dt
            % 
            %   coded by Manuel Diaz, NTU, 2012.12.18
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.

            % Along x
            dq = FDMethods2.Upwindflux(ax,q,[0,1])/dx;
            
            % Impose BC
            dq(:,1)=dq(:,2); dq(:,end)=dq(:,end-1);
            
            % The Residual
            dF = dq;
        end
        
        function dF = Upwindresidual2d(q,ax,ay,dx,dy)
                % Compute the Residual for 2d wave equation using Upwind Method
                %
                %           Residual: RHS of dq/dt
                %
                %   coded by Manuel Diaz, NTU, 2012.12.18
                %
                % ax,zy: scalar advection velocities in x and y directions respectively.
                % ensure array is 2-D!
                qo(:,:) = q;
                
                % Along x
                dq=   FDMethods2.Upwindflux(ax,qo,[0,1])/dx;
                
                % Along y
                dq=dq+FDMethods2.Upwindflux(ay,qo,[1,0])/dy;
                
                % Impose BC
                dq(:,1)=dq(:,2); dq(:,end)=dq(:,end-1);
                dq(1,:)=dq(2,:); dq(end,:)=dq(end-1,:);
                
                % The Residual
                dF = dq;
            end
               
        function [dF] = Upwindflux(a,u,m,p)
            % Vectorized version of Two sided upwinding flux
            %
            % coded by Manuel Diaz, NTU, 2012.12.20
            %
            % Domain cells (I{i}) reference:
            %
            %                |           |   u(i)    |           |
            %                |  u(i-1)   |___________|           |
            %                |___________|           |   u(i+1)  |
            %                |           |           |___________|
            %             ...|-----0-----|-----0-----|-----0-----|...
            %                |    i-1    |     i     |    i+1    |
            %                |-         +|-         +|-         +|
            %              i-3/2       i-1/2       i+1/2       i+3/2
            %
            % Vector variables, (periodic BC behavior included)
            %um =circshift(uo,ind);  %u(i-1)
            %up =circshift(uo,-ind); %u(i+1)

            % for flux splitting
            a_p = max(a,0); % a{+}
            a_m = min(a,0); % a{-}

            % Double sided upwind fluxes
            fl = a_p.*u(m) + a_m.*u; % left
            fr = a_p.*u + a_m.*u(p); % right

            % df values
            dF = fr-fl;
        end
        
        function dF = TVDresidual1d(q,ax,dt,dx,limiter)
            % Compute the Residual for 2d wave equation using TVD method
            %
            %           Residual: RHS of dq/dt
            %
            %   coded by Manuel Diaz, NTU, 2012.12.18
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            
            % Along x

            % Compute smoothness factor r
            [r] = FDMethods2.theta(q,ax,[0,1]);

            % Compute the Flux Limiter.
            [phix] = FDMethods2.fluxlimiter(r,limiter);
            
            % flux in x.
            dq = FDMethods2.TVDflux(ax,q,dt/dx,phix,[0,1]);
            
            % Impose BC
            dq(:,1)=dq(:,2); dq(:,end)=dq(:,end-1);
            
            % The Residual
            dF = dq/dx;
        end
        
        function dF = TVDresidual2d(q,ax,ay,dt,dx,dy,limiter)
            % Compute the Residual for 2d wave equation using TVD method
            %
            %           Residual: RHS of dq/dt
            %
            %   coded by Manuel Diaz, NTU, 2012.12.18
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            % ensure array is 2-D!
            qo(:,:) = q;
            
            % Along x
            
            % Compute smooth factor along x
            [r_x] = FDMethods2.theta(qo,ax,[0,1]);
            
            % Compute the Flux Limiter.
            [phix] = FDMethods2.fluxlimiter(r_x,limiter);
            
            % flux in x.
            dq = FDMethods2.TVDflux(ax,qo,dt/dx,phix,[0,1])/dx;
            
            % Along x
            
            % Compute smooth factor along x
            [r_y] = FDMethods2.theta(qo,ay,[1,0]);
            
            % Compute the Flux Limiter.
            [phiy] = FDMethods2.fluxlimiter(r_y,limiter);
            
            % flux in x.
            dq=dq + FDMethods2.TVDflux(ay,qo,dt/dy,phiy,[1,0])/dy;
            
            % Impose BC
            dq(:,1)=dq(:,2); dq(:,end)=dq(:,end-1);
            dq(1,:)=dq(2,:); dq(end,:)=dq(end-1,:);
            
            % The Residual
            dF = dq;
        end
        
        function [r_x] = theta(uo,a,ind)
            % Vectorized version of TVD theta: smooth measurement function
            %
            % coded by Manuel Diaz, NTU, 2012.12.20
            
            r_x=zeros(size(uo));
            
            up =circshift(uo,-ind);   %u(i+1)
            upp=circshift(uo,-2*ind); %u(i+2)
            um =circshift(uo,+ind);   %u(i-1)
            
            ido=find(uo==up); % smooth cells
            r_x(ido)=1;
            
            idx=find(uo~=up); % non-smooth cells
            if a>0
                r_x(idx)=(uo(idx)-um(idx))./(up(idx)-uo(idx));
            else % a<0
                r_x(idx)=(upp(idx)-up(idx))./(up(idx)-uo(idx));
            end
        end
        
        function phi_x = fluxlimiter(r_x,limiter)
            % Compute by limiter case:
            switch limiter
                case 'Van Leer'
                    %fprintf('\n using Van Leer Limiter \n\n');
                    phi_x = (r_x + abs(r_x))./(1 + abs(r_x));
                    
                case 'Superbee'
                    %fprintf('\n using Superbee Limiter \n\n');
                    %phi = max(0,min(2.*r,0),min(r,2));
                    phi_x_star = max(min(2.*r_x,0),min(r_x,2));
                    phi_x = max(0,phi_x_star);
                    
                case 'minmod'
                    %fprintf('\n using Minmod Limiter \n\n');
                    phi_x = max(0,min(1,r_x));
                    
                case 'Koren'
                    %fprintf('\n using Koren Limiter \n\n');
                    %phi = max(0,min(2*r,(2/3)*r+1/3,2));
                    phi_x_star = min(2*r_x,(2/3)*r_x+1/3);
                    phi_x_hat  = min(phi_x_star,2);
                    phi_x = max(0,phi_x_hat);
                    
                otherwise
                    error('case not supported');
            end
        end
        
        function [dF] = TVDflux(a,uo,dtdx,phix,ind)
            % Vectorized version of scalar TVD flux
            %
            % coded by Manuel Diaz, NTU, 2012.12.20
            %
            % Domain cells (I{i}) reference:
            %
            %                |           |   u(i)    |           |
            %                |  u(i-1)   |___________|           |
            %                |___________|           |   u(i+1)  |
            %                |           |           |___________|
            %             ...|-----0-----|-----0-----|-----0-----|...
            %                |    i-1    |     i     |    i+1    |
            %                |-         +|-         +|-         +|
            %              i-3/2       i-1/2       i+1/2       i+3/2
            %
            % vector variables, (periodic BC behavior included)
            um = circshift(uo,ind);  %u(i-1)
            up = circshift(uo,-ind); %u(i+1)
            phixm = circshift(phix,ind); %phix(i-1)
            
            % for flux splitting
            a_p = max(a,0); % a{+}
            a_m = min(a,0); % a{-}
            
            % Compute TVD fluxes
            % Notation: l = low flux  & h = high flux
            %           s = Left flux & r = Right flux
            
            % Flux to the right(r)
            frl = a_p*uo + a_m.*up;
            frh = a.*(uo + up)/2 - dtdx*(a.^2).*(up - uo)/2;
            fr  = frl + phix.*(frh - frl);
            
            % FLux to the left (s)
            fsl = a_p*um + a_m*uo;
            fsh = a.*(um + uo)/2 - dtdx*(a.^2)*(uo - um)/2;
            fs  = fsl + phixm.*(fsh - fsl);
            
            % df values
            dF = fr-fs;
        end
        
        function dF = WENO5residual1d(q,ax,dx)
            % Compute the Residual for 1d wave equation using WENO5
            %
            %               Residual = dq/dx
            %
            %     coded by Manuel Diaz, NTU, 2012.08.20
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            
            % Along x
            [axp,axm]=FDMethods2.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods2.WENO5_reconstruction(axm,q,[0 -1]);
            dq=dq+FDMethods2.WENO5_reconstruction(axp,q,[0 +1]);
            
            % Impose BC
            dq(:,3)=dq(:,4); dq(:,end-2)=dq(:,end-3);
            dq(:,2)=dq(:,3); dq(:,end-1)=dq(:,end-2);
            dq(:,1)=dq(:,2); dq(:, end )=dq(:,end-1);
                        
            % The Residual
            dF = -dq/dx;
        end
        
        function dF = WENO5residual2d(q,ax,ay,dx,dy)
            % Compute the Residual for 2d wave equation using Upwind
            %
            %             Residual: RHS of dq/dt
            %
            % 	coded by Manuel Diaz, NTU, 2012.12.18
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            % ensure array is 2-D!
            qo(:,:) = q;
            
            % Along x
            [axp,axm]=FDMethods2.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods2.WENO5_reconstruction(axm,qo,[0 -1])/dx;
            dq=dq+FDMethods2.WENO5_reconstruction(axp,qo,[0 +1])/dx;
            
            % Along y
            [ayp,aym]=FDMethods2.fluxsplitting_scalar(ay,'LF');
            aym=circshift(aym,[1 0]);
            
            dq=dq+FDMethods2.WENO5_reconstruction(aym,qo,[-1 0])/dy;
            dq=dq+FDMethods2.WENO5_reconstruction(ayp,qo,[+1 0])/dy;
            
            % Impose BC in x
            dq(:,3)=dq(:,4); dq(:,end-2)=dq(:,end-3);
            dq(:,2)=dq(:,3); dq(:,end-1)=dq(:,end-2);
            dq(:,1)=dq(:,2); dq(:, end )=dq(:,end-1);
            
            % Impose BC in y
            dq(3,:)=dq(4,:); dq(end-2,:)=dq(end-3,:);
            dq(2,:)=dq(3,:); dq(end-1,:)=dq(end-2,:);
            dq(1,:)=dq(2,:); dq( end ,:)=dq(end-1,:);
            
            % The Residual
            dF = -dq;
        end
        
        function outputf = WENO5_reconstruction(a,qo,ind)
            %
            % coded by Manuel Diaz, NTU, 2012.12.20
            %
            % Domain cells (I{i}) reference:
            %
            %                |           |   u(i)    |           |
            %                |  u(i-1)   |___________|           |
            %                |___________|           |   u(i+1)  |
            %                |           |           |___________|
            %             ...|-----0-----|-----0-----|-----0-----|...
            %                |    i-1    |     i     |    i+1    |
            %                |-         +|-         +|-         +|
            %              i-3/2       i-1/2       i+1/2       i+3/2
            %
            % ENO stencils (S{r}) reference:
            %
            %
            %                               |___________S2__________|
            %                               |                       |
            %                       |___________S1__________|       |
            %                       |                       |       |
            %               |___________S0__________|       |       |
            %             ..|---o---|---o---|---o---|---o---|---o---|...
            %               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
            %                                      -|
            %                                     i+1/2
            %
            %
            %               |___________S0__________|
            %               |                       |
            %               |       |___________S1__________|
            %               |       |                       |
            %               |       |       |___________S2__________|
            %             ..|---o---|---o---|---o---|---o---|---o---|...
            %               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
            %                               |+
            %                             i-1/2
            %
            % WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Note: by using circshift over our domain, we are implicitly creating
            % favorable code that includes periodical boundary conditions. MD. 12.2012.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Set Variables
            qmm=circshift(qo,ind*2);
            qm =circshift(qo,ind);
            qp =circshift(qo,-ind);
            qpp=circshift(qo,-ind*2);
            
            % Reconstruction Polynomials
            up0 = (2*qmm - 7*qm + 11*qo)/6.;
            up1 = ( -qm  + 5*qo + 2*qp )/6.;
            up2 = (2*qo  + 5*qp - qpp  )/6.;
            
            % Smooth parameters
            b0 = 13/12*(qmm-2*qm+qo ).^2 + 1/4*(qmm-4*qm+3*qo).^2;
            b1 = 13/12*(qm -2*qo+qp ).^2 + 1/4*(qm-qp).^2;
            b2 = 13/12*(qo -2*qp+qpp).^2 + 1/4*(3*qo-4*qp+qpp).^2;
            
            % Constants
            g0=1/10; g1=6/10; g2=3/10; epsilon=1e-6;
            
            % weigths
            wt0 = g0./(epsilon+b0).^2;
            wt1 = g1./(epsilon+b1).^2;
            wt2 = g2./(epsilon+b2).^2;
            sum_wt=wt0+wt1+wt2;
            
            % Non-linear weigths
            w0 = wt0./sum_wt;
            w1 = wt1./sum_wt;
            w2 = wt2./sum_wt;
            
            % WENO5 polynomial
            qp = w0.*up0 + w1.*up1 + w2.*up2;
            
            % Compute flux
            fp = qp.*a; outputf = -(fp-circshift(fp,ind));
            % End of reconstruction.
        end
        
        function dF = WENO7residual1d(q,ax,dx)
            % Compute the Residual for 1d wave equation using WENO5
            %
            %               Residual = dq/dx
            %
            %     coded by Manuel Diaz, NTU, 2012.08.20
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
                       
            % Along x
            [axp,axm]=FDMethods2.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods2.WENO7_reconstruction(axm,q,[0 -1]);
            dq=dq+FDMethods2.WENO7_reconstruction(axp,q,[0 +1]);
            
            % Impose BC
            dq(:,4)=dq(:,5); dq(:,end-3)=dq(:,end-4);
            dq(:,3)=dq(:,4); dq(:,end-2)=dq(:,end-3);
            dq(:,2)=dq(:,3); dq(:,end-1)=dq(:,end-2);
            dq(:,1)=dq(:,2); dq(:, end )=dq(:,end-1);
                        
            % The Residual
            dF = -dq/dx;
        end
        
        function dF = WENO7residual2d(q,ax,ay,dx,dy)
            % Compute the Residual for 2d wave equation using Upwind
            %
            %             Residual: RHS of dq/dt
            %
            % 	coded by Manuel Diaz, NTU, 2012.12.18
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            % ensure array is 2-D!
            qo(:,:) = q;
            
            % Along x
            [axp,axm]=FDMethods2.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods2.WENO7_reconstruction(axm,qo,[0 -1])/dx;
            dq=dq+FDMethods2.WENO7_reconstruction(axp,qo,[0 +1])/dx;
            
            % Along y
            [ayp,aym]=FDMethods2.fluxsplitting_scalar(ay,'LF');
            aym=circshift(aym,[1 0]);
            
            dq=dq+FDMethods2.WENO7_reconstruction(aym,qo,[-1 0])/dy;
            dq=dq+FDMethods2.WENO7_reconstruction(ayp,qo,[+1 0])/dy;
            
            % Impose BC in x
            dq(:,4)=dq(:,5); dq(:,end-3)=dq(:,end-4);
            dq(:,3)=dq(:,4); dq(:,end-2)=dq(:,end-3);
            dq(:,2)=dq(:,3); dq(:,end-1)=dq(:,end-2);
            dq(:,1)=dq(:,2); dq(:, end )=dq(:,end-1);
            
            % Impose BC in y
            dq(4,:)=dq(5,:); dq(end-3,:)=dq(end-4,:);
            dq(3,:)=dq(4,:); dq(end-2,:)=dq(end-3,:);
            dq(2,:)=dq(3,:); dq(end-1,:)=dq(end-2,:);
            dq(1,:)=dq(2,:); dq( end ,:)=dq(end-1,:);
            
            % The Residual
            dF = -dq;
        end
        
        function outputf = WENO7_reconstruction(a,qo,ind)
            %
            % coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
            %
            % Domain cells (I{i}) reference:
            %
            %                |           |   u(i)    |           |
            %                |  u(i-1)   |___________|           |
            %                |___________|           |   u(i+1)  |
            %                |           |           |___________|
            %             ...|-----0-----|-----0-----|-----0-----|...
            %                |    i-1    |     i     |    i+1    |
            %                |-         +|-         +|-         +|
            %              i-3/2       i-1/2       i+1/2       i+3/2
            %
            % ENO stencils (S{r}) reference:
            %
            %                               |_______________S3______________|
            %                               |                               |
            %                       |______________S2_______________|       |
            %                       |                               |       |
            %               |______________S1_______________|       |       |
            %               |                               |       |       |
            %       |_______________S0______________|       |       |       |
            %     ..|---o---|---o---|---o---|---o---|---o---|---o---|---o---|...
            %       | I{i-3}| I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}| I{i+3}|
            %                                      -|
            %                                     i+1/2
            %
            %       |______________S0_______________|
            %       |                               |
            %       |       |______________S1_______________|
            %       |       |                               |
            %       |       |       |______________S2_______________|
            %       |       |       |                               |
            %       |       |       |       |_______________S3______________|
            %     ..|---o---|---o---|---o---|---o---|---o---|---o---|---o---|...
            %       | I{i-3}| I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|| I{i+3}
            %                               |+
            %                             i-1/2
            %
            % WENO stencil: S{i} = [ I{i-3},...,I{i+3} ]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Note: by using circshift over our domain, we are implicitly creating
            % faqorable code that includes peridical boundary conditions. MD. 12.2012.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set Variables
            qmmm=circshift(qo,ind*3);
            qmm =circshift(qo,ind*2);
            qm  =circshift(qo,ind);
            qp  =circshift(qo,-ind);
            qpp =circshift(qo,-ind*2);
            qppp=circshift(qo,-ind*3);
            
            % Reconstruction Polynomials
            p0 = (-3*qmmm + 13*qmm - 23*qm  + 25*qo  )/12;
            p1 = ( 1*qmm  - 5*qm   + 13*qo  +  3*qp  )/12;
            p2 = (-1*qm   + 7*qo   +  7*qp  -  1*qpp )/12;
            p3 = ( 3*qo   + 13*qp  -  5*qpp +  1*qppp)/12;
            
            % Smooth Indicators
            B0 = qm.*(134241*qm-114894*qo) ...
                +qmmm.*(56694*qm-47214*qmm+6649*qmmm-22778*qo)...
                +25729*qo.^2  +qmm.*(-210282*qm+85641*qmm+86214*qo);
            B1 = qo.*(41001*qo-30414*qp) ...
                +qmm.*(-19374*qm+3169*qmm+19014*qo-5978*qp)...
                +6649*qp.^2   +qm.*(33441*qm-70602*qo+23094*qp);
            B2 = qp.*(33441*qp-19374*qpp) ...
                +qm.*(6649*qm-30414*qo+23094*qp-5978*qpp)...
                +3169*qpp.^2  +qo.*(41001*qo-70602*qp+19014*qpp);
            B3 = qpp.*(85641*qpp-47214*qppp) ...
                +qo.*(25729*qo-114894*qp+86214*qpp-22778*qppp)...
                +6649*qppp.^2 +qp.*(134241*qp-210282*qpp+56694*qppp);
            
            % Constants
            epsilon = 1e-6; g0 = 1/35; g1 = 12/35; g2 = 18/35; g3 = 4/35;
            
            % Alpha weights
            alpha0 = g0./(epsilon + B0/2880).^2;
            alpha1 = g1./(epsilon + B1/2880).^2;
            alpha2 = g2./(epsilon + B2/2880).^2;
            alpha3 = g3./(epsilon + B3/2880).^2;
            alphasum = alpha0 + alpha1 + alpha2 + alpha3;
            
            % Non-linear weigths
            w0 = alpha0./alphasum;
            w1 = alpha1./alphasum;
            w2 = alpha2./alphasum;
            w3 = alpha3./alphasum;
            
            % Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
            qp = w0.*p0 + w1.*p1 + w2.*p2 + w3.*p3;
            
            % Compute flux
            fp = qp.*a; outputf = -(fp-circshift(fp,ind));
            % End of reconstruction.
        end
        
        function [ap,am]=fluxsplitting_scalar(a,strategy)
            % fort splitting the flux into right-going and left-going
            % OUTPUTS:
            % ap: positive dflux, v^{+}, which corresponds to f_{i+1/2}^{-}
            % an: negative dflux  v^{-}, which corresponds to f_{i+1/2}^{+}
            
            switch strategy
                case 'Upwind'
                    ap= max(a,0); %a^{+}
                    am=-min(a,0); %a^{-}
                case 'LF'
                    au=max(abs(a(:)));
                    ap =  0.5*(a+au); %a^{+}
                    am = -0.5*(a-au); %a^{-}
                otherwise
                    error('case not available')
            end
        end
    end % static methods
end % class