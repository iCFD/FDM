classdef FDMethods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Finite Difference Methods for Linear conservation equations
    %
    %   This class cotains fixed subroutines for computing the cartesian
    %   grids and the residuals for the following methods:
    %               1st order Double sided Upwind,
    %               2nd order TVD,
    %               5th order WENO reconstruction.
    %               7th order WENO reconstruction.
    %
    %             coded by Manuel Diaz, NTU, 2013.12.22
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convention:
    %   circshift(uo,[0 +3]) = u(mmm{:});  %u(i-3)
    %   circshift(uo,[0 +2]) = u(mm{:});   %u(i-2)
    %   circshift(uo,[0 +1]) = u(m{:});    %u(i-1)
    %   circshift(uo,[0 -1]) = u(p{:});    %u(i+1)
    %   circshift(uo,[0 -2]) = u(pp{:});   %u(i+2)
    %   circshift(uo,[0 -3]) = u(ppp{:});  %u(i+3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        mmmx
        mmx
        mx
        px
        ppx
        pppx
        mmmy
        mmy
        my
        py
        ppy
        pppy
        mmmz
        mmz
        mz
        pz
        ppz
        pppz
    end
    
    methods
        function obj = FDMethods(BC,NX) % Constructor
            % Check Matlab Tips and Tricks for more reference.
            dim = size(NX,2); % max(varargin)=3 is expected
            switch BC
                case 'periodic'
                    switch dim
                        case 1 %'1d'
                            dimension = 1; nx = NX(1);
                            idx = repmat({':'},dimension,1);
                            obj.mmmx = idx;
                            obj.mmx = idx;
                            obj.mx = idx;
                            obj.px = idx;
                            obj.ppx = idx;
                            obj.pppx = idx;
                            % x
                            obj.mmmx{2} = [ nx-2:nx 1:nx-3];
                            obj.mmx{2}  = [ nx-1:nx 1:nx-2];
                            obj.mx{2}   = [ nx 1:nx-1];
                            obj.px{2}   = [ 2:nx 1];
                            obj.ppx{2}  = [ 3:nx 1:2];
                            obj.pppx{2}	= [ 4:nx 1:3];
                        case 2 %'2d'
                            dimension = 2; nx = NX(1); ny = NX(2);
                            idx = repmat({':'},dimension,1);
                            obj.mmmx = idx;
                            obj.mmx = idx;
                            obj.mx = idx;
                            obj.px = idx;
                            obj.ppx = idx;
                            obj.pppx = idx;
                            obj.mmmy = idx;
                            obj.mmy = idx;
                            obj.my = idx;
                            obj.py = idx;
                            obj.ppy = idx;
                            obj.pppy = idx;
                            % x
                            obj.mmmx{2} = [ nx-2:nx 1:nx-3];
                            obj.mmx{2}  = [ nx-1:nx 1:nx-2];
                            obj.mx{2}   = [ nx 1:nx-1];
                            obj.px{2}   = [ 2:nx 1];
                            obj.ppx{2}  = [ 3:nx 1:2];
                            obj.pppx{2}	= [ 4:nx 1:3];
                            % y
                            obj.mmmy{1} = [ ny-2:ny 1:ny-3];
                            obj.mmy{1}  = [ ny-1:ny 1:ny-2];
                            obj.my{1}   = [ ny 1:ny-1];
                            obj.py{1}   = [ 2:ny 1];
                            obj.ppy{1}  = [ 3:ny 1:2];
                            obj.pppy{1}	= [ 4:ny 1:3];
                        case 3 %'3d'
                            dimension = 3; nx = NX(1); ny = NX(2); nz = NX(3);
                            idx = repmat({':'},dimension,1);
                            obj.mmmx = idx;
                            obj.mmx = idx;
                            obj.mx = idx;
                            obj.px = idx;
                            obj.ppx = idx;
                            obj.pppx = idx;
                            obj.mmmy = idx;
                            obj.mmy = idx;
                            obj.my = idx;
                            obj.py = idx;
                            obj.ppy = idx;
                            obj.pppy = idx;
                            obj.mmmz = idx;
                            obj.mmz = idx;
                            obj.mz = idx;
                            obj.pz = idx;
                            obj.ppz = idx;
                            obj.pppz = idx;
                            % x
                            obj.mmmx{2} = [ nx-2:nx 1:nx-3];
                            obj.mmx{2}  = [ nx-1:nx 1:nx-2];
                            obj.mx{2}   = [ nx 1:nx-1];
                            obj.px{2}   = [ 2:nx 1];
                            obj.ppx{2}  = [ 3:nx 1:2];
                            obj.pppx{2}	= [ 4:nx 1:3];
                            % y
                            obj.mmmy{1} = [ ny-2:ny 1:ny-3];
                            obj.mmy{1}  = [ ny-1:ny 1:ny-2];
                            obj.my{1}   = [ ny 1:ny-1];
                            obj.py{1}   = [ 2:ny 1];
                            obj.ppy{1}  = [ 3:ny 1:2];
                            obj.pppy{1}	= [ 4:ny 1:3];
                            % z
                            obj.mmmz{3} = [ nz-2:nz 1:nz-3];
                            obj.mmz{3}  = [ nz-1:nz 1:nz-2];
                            obj.mz{3}   = [ nz 1:nz-1];
                            obj.pz{3}   = [ 2:nz 1];
                            obj.ppz{3}  = [ 3:nz 1:2];
                            obj.pppz{3}	= [ 4:nz 1:3];
                        otherwise
                            error('incorrect arguments');
                    end
                case 'non-periodic'
                    switch dim
                        case 1 %'1d'
                            dimension = 1; nx = NX(1); 
                            idx = repmat({':'},dimension,1);
                            obj.mmmx = idx;
                            obj.mmx = idx;
                            obj.mx = idx;
                            obj.px = idx;
                            obj.ppx = idx;
                            obj.pppx = idx;
                            % x
                            obj.mmmx{2} = [ 1 1 1 1:nx-3];
                            obj.mmx{2}  = [ 1 1 1:nx-2];
                            obj.mx{2}   = [ 1 1:nx-1];
                            obj.px{2}   = [ 2:nx nx];
                            obj.ppx{2}  = [ 3:nx nx nx];
                            obj.pppx{2}	= [ 4:nx nx nx nx];
                        case 2 %'2d'
                            dimension = 2; nx = NX(1); ny = NX(2); 
                            idx = repmat({':'},dimension,1);
                            obj.mmmx = idx;
                            obj.mmx = idx;
                            obj.mx = idx;
                            obj.px = idx;
                            obj.ppx = idx;
                            obj.pppx = idx;
                            obj.mmmy = idx;
                            obj.mmy = idx;
                            obj.my = idx;
                            obj.py = idx;
                            obj.ppy = idx;
                            obj.pppy = idx;
                            % x
                            obj.mmmx{2} = [ 1 1 1 1:nx-3];
                            obj.mmx{2}  = [ 1 1 1:nx-2];
                            obj.mx{2}   = [ 1 1:nx-1];
                            obj.px{2}   = [ 2:nx nx];
                            obj.ppx{2}  = [ 3:nx nx nx];
                            obj.pppx{2}	= [ 4:nx nx nx nx];
                            % y
                            obj.mmmy{1} = [ 1 1 1 1:ny-3];
                            obj.mmy{1}  = [ 1 1 1:ny-2];
                            obj.my{1}   = [ 1 1:ny-1];
                            obj.py{1}   = [ 2:ny ny];
                            obj.ppy{1}  = [ 3:ny ny ny];
                            obj.pppy{1}	= [ 4:ny ny ny ny];
                        case 3 %'3d'
                            dimension = 3; nx = NX(1); ny = NX(2); nz = NX(3);
                            idx = repmat({':'},dimension,1);
                            obj.mmmx = idx;
                            obj.mmx = idx;
                            obj.mx = idx;
                            obj.px = idx;
                            obj.ppx = idx;
                            obj.pppx = idx;
                            obj.mmmy = idx;
                            obj.mmy = idx;
                            obj.my = idx;
                            obj.py = idx;
                            obj.ppy = idx;
                            obj.pppy = idx;
                            obj.mmmz = idx;
                            obj.mmz = idx;
                            obj.mz = idx;
                            obj.pz = idx;
                            obj.ppz = idx;
                            obj.pppz = idx;
                            % x
                            obj.mmmx{2} = [ 1 1 1 1:nx-3];
                            obj.mmx{2}  = [ 1 1 1:nx-2];
                            obj.mx{2}   = [ 1 1:nx-1];
                            obj.px{2}   = [ 2:nx nx];
                            obj.ppx{2}  = [ 3:nx nx nx];
                            obj.pppx{2}	= [ 4:nx nx nx nx];
                            % y
                            obj.mmmy{1} = [ 1 1 1 1:ny-3];
                            obj.mmy{1}  = [ 1 1 1:ny-2];
                            obj.my{1}   = [ 1 1:ny-1];
                            obj.py{1}   = [ 2:ny ny];
                            obj.ppy{1}  = [ 3:ny ny ny];
                            obj.pppy{1}	= [ 4:ny ny ny ny];
                            % z
                            obj.mmmz{3} = [ 1 1 1 1:nz-3];
                            obj.mmz{3}  = [ 1 1 1:nz-2];
                            obj.mz{3}   = [ 1 1:nz-1];
                            obj.pz{3}   = [ 2:nz nz];
                            obj.ppz{3}  = [ 3:nz nz nz];
                            obj.pppz{3}	= [ 4:nz nz nz nz];
                        otherwise
                            error('incorrect arguments');
                    end
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Residuals Collection: %
    %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function dF = Upwindresidual1d(obj,q,ax,dx)
            % Compute the Residual for 2d wave equation using Upwind Method
            %
            %           Residual: RHS of dq/dt
            % 
            %   coded by Manuel Diaz, NTU, 2012.12.18
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.

            % Along x
            dq = FDMethods.Upwindflux(ax,q,obj.mx,obj.px)/dx;
            
            % The Residual
            dF = dq;
        end
        
        function dF = Upwindresidual2d(obj,q,ax,ay,dx,dy)
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
            %dq=   FDMethods.Upwindflux(ax,qo,[0,1])/dx;
            dq=   FDMethods.Upwindflux(ax,qo,obj.mx,obj.px)/dx;

            % Along y
            %dq=dq+FDMethods.Upwindflux(ay,qo,[1,0])/dy;
            dq=dq+FDMethods.Upwindflux(ay,qo,obj.my,obj.py)/dy;

            % The Residual
            dF = dq;
        end
            
        function dF = TVDresidual1d(obj,q,ax,dt,dx,limiter)
            % Compute the Residual for 2d wave equation using TVD method
            %
            %           Residual: RHS of dq/dt
            %
            %   coded by Manuel Diaz, NTU, 2012.12.18
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            
            % Along x

            % Compute smoothness factor r
            [r] = FDMethods.theta(q,ax,obj.mx,obj.px,obj.ppx);

            % Compute the Flux Limiter.
            [phix] = FDMethods.fluxlimiter(r,limiter);
            
            % flux in x.
            dq = FDMethods.TVDflux(ax,q,dt/dx,phix,obj.mx,obj.px);
            
            % The Residual
            dF = dq/dx;
        end
        
        function dF = TVDresidual2d(obj,q,ax,ay,dt,dx,dy,limiter)
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
            [r_x] = FDMethods.theta(qo,ax,obj.mx,obj.px,obj.ppx);
            
            % Compute the Flux Limiter.
            [phix] = FDMethods.fluxlimiter(r_x,limiter);
            
            % flux in x.
            dq = FDMethods.TVDflux(ax,qo,dt/dx,phix,obj.mx,obj.px)/dx;
            
            % Along x
            
            % Compute smooth factor along x
            [r_y] = FDMethods.theta(qo,ay,obj.my,obj.py,obj.ppy);
            
            % Compute the Flux Limiter.
            [phiy] = FDMethods.fluxlimiter(r_y,limiter);
            
            % flux in x.
            dq=dq + FDMethods.TVDflux(ay,qo,dt/dy,phiy,obj.my,obj.py)/dy;
            
            % The Residual
            dF = dq;
        end
        
        function dF = WENO3residual1d(obj,q,ax,dx)
            % Compute the Residual for 1d wave equation using WENO5
            %
            %               Residual = dq/dx
            %
            %     coded by Manuel Diaz, NTU, 2012.08.20
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            
            % Along x
            [axp,axm]=FDMethods.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods.WENO3_reconstruction(axm,q,obj.px,obj.mx);
            dq=dq+FDMethods.WENO3_reconstruction(axp,q,obj.mx,obj.px);
                        
            % The Residual
            dF = -dq/dx;
        end
        
        function dF = WENO3residual2d(obj,q,ax,ay,dx,dy)
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
            [axp,axm]=FDMethods.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods.WENO3_reconstruction(axm,qo,obj.px,obj.mx)/dx;
            dq=dq+FDMethods.WENO3_reconstruction(axp,qo,obj.mx,obj.px)/dx;
            
            % Along y
            [ayp,aym]=FDMethods.fluxsplitting_scalar(ay,'LF');
            aym=circshift(aym,[1 0]);
                        
            dq=dq+FDMethods.WENO3_reconstruction(aym,qo,obj.py,obj.my)/dy;
            dq=dq+FDMethods.WENO3_reconstruction(ayp,qo,obj.my,obj.py)/dy;
                        
            % The Residual
            dF = -dq;
        end
        
        function dF = WENO5residual1d(obj,q,ax,dx)
            % Compute the Residual for 1d wave equation using WENO5
            %
            %               Residual = dq/dx
            %
            %     coded by Manuel Diaz, NTU, 2012.08.20
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
            
            % Along x
            [axp,axm]=FDMethods.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods.WENO5_reconstruction(axm,q,obj.ppx,obj.px,obj.mx,obj.mmx);
            dq=dq+FDMethods.WENO5_reconstruction(axp,q,obj.mmx,obj.mx,obj.px,obj.ppx);
                        
            % The Residual
            dF = -dq/dx;
        end
        
        function dF = WENO5residual2d(obj,q,ax,ay,dx,dy)
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
            [axp,axm]=FDMethods.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods.WENO5_reconstruction(axm,qo,obj.ppx,obj.px,obj.mx,obj.mmx)/dx;
            dq=dq+FDMethods.WENO5_reconstruction(axp,qo,obj.mmx,obj.mx,obj.px,obj.ppx)/dx;
            
            % Along y
            [ayp,aym]=FDMethods.fluxsplitting_scalar(ay,'LF');
            aym=circshift(aym,[1 0]);
                        
            dq=dq+FDMethods.WENO5_reconstruction(aym,qo,obj.ppy,obj.py,obj.my,obj.mmy)/dy;
            dq=dq+FDMethods.WENO5_reconstruction(ayp,qo,obj.mmy,obj.my,obj.py,obj.ppy)/dy;
                        
            % The Residual
            dF = -dq;
        end
        
        function dF = WENO7residual1d(obj,q,ax,dx)
            % Compute the Residual for 1d wave equation using WENO5
            %
            %               Residual = dq/dx
            %
            %     coded by Manuel Diaz, NTU, 2012.08.20
            %
            % ax,zy: scalar advection velocities in x and y directions respectively.
                       
            % Along x
            [axp,axm]=FDMethods.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
                       
            dq=   FDMethods.WENO7_reconstruction(axm,q,obj.pppx,obj.ppx,obj.px,obj.mx,obj.mmx,obj.mmmx);
            dq=dq+FDMethods.WENO7_reconstruction(axp,q,obj.mmmx,obj.mmx,obj.mx,obj.px,obj.ppx,obj.pppx);
                                    
            % The Residual
            dF = -dq/dx;
        end
        
        function dF = WENO7residual2d(obj,q,ax,ay,dx,dy)
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
            [axp,axm]=FDMethods.fluxsplitting_scalar(ax,'LF');
            axm=circshift(axm,[0 1]);
            
            dq=   FDMethods.WENO7_reconstruction(axm,qo,obj.pppx,obj.ppx,obj.px,obj.mx,obj.mmx,obj.mmmx)/dx;
            dq=dq+FDMethods.WENO7_reconstruction(axp,qo,obj.mmmx,obj.mmx,obj.mx,obj.px,obj.ppx,obj.pppx)/dx;
            
            % Along y
            [ayp,aym]=FDMethods.fluxsplitting_scalar(ay,'LF');
            aym=circshift(aym,[1 0]);
                        
            dq=dq+FDMethods.WENO7_reconstruction(aym,qo,obj.pppy,obj.ppy,obj.py,obj.my,obj.mmy,obj.mmmy)/dy;
            dq=dq+FDMethods.WENO7_reconstruction(ayp,qo,obj.mmmy,obj.mmy,obj.my,obj.py,obj.ppy,obj.pppy)/dy;
                        
            % The Residual
            dF = -dq;
        end
        
    end % methods

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Build Cartesian Grids %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    end % static methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Advection subroutines %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static)
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

            % for flux splitting
            a_p = max(a,0); % a{+}
            a_m = min(a,0); % a{-}

            % Double sided upwind fluxes
            fl = a_p.*u(m{:}) + a_m.*u; % left
            fr = a_p.*u + a_m.*u(p{:}); % right

            % df values
            dF = fr-fl;
        end
                
        function [r_x] = theta(uo,a,m,p,pp)
            % Vectorized version of TVD theta: smooth measurement function
            %
            % coded by Manuel Diaz, NTU, 2012.12.20
            
            r_x=zeros(size(uo));
            
            up =uo(p{:});  %u(i+1)
            upp=uo(pp{:}); %u(i+2)
            um =uo(m{:});  %u(i-1)
            
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
        
        function [dF] = TVDflux(a,u,dtdx,phix,m,p)
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
            %um = circshift(uo,ind);  %u(i-1)
            %up = circshift(uo,-ind); %u(i+1)
            %phixm = circshift(phix,ind); %phix(i-1)
            
            % for flux splitting
            a_p = max(a,0); % a{+}
            a_m = min(a,0); % a{-}
            
            % Compute TVD fluxes
            % Notation: l = low flux  & h = high flux
            %           s = Left flux & r = Right flux
            
            % Flux to the right(r)
            frl = a_p*u + a_m.*u(p{:});
            frh = a.*(u + u(p{:}))/2 - dtdx*(a.^2).*(u(p{:}) - u)/2;
            fr  = frl + phix.*(frh - frl);
            
            % FLux to the left (s)
            fsl = a_p*u(m{:}) + a_m*u;
            fsh = a.*(u(m{:}) + u)/2 - dtdx*(a.^2)*(u - u(m{:}))/2;
            fs  = fsl + phix(m{:}).*(fsh - fsl);
            
            % df values
            dF = fr-fs;
        end
        
        function outputf = WENO3_reconstruction(a,q,m,p)
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
            %                               |______S1_______|
            %                               |               |
            %                       |______S0_______|       |
            %             ..|---o---|---o---|---o---|---o---|---o---|...
            %               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
            %                                      -|
            %                                     i+1/2
            %
            %
            %                       |______S0_______|
            %                       |               |
            %                       |       |______S1_______|
            %             ..|---o---|---o---|---o---|---o---|---o---|...
            %               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
            %                               |+
            %                             i-1/2
            %
            % WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            % Reconstruction Polynomials
            up0 = (-q(m{:}) + 3*q     )/2.;
            up1 = ( q       + q(p{:}) )/2.;
            
            % Smooth indicators parameters
            b0 = (q(m{:}) - q).^2;
            b1 = (q - q(p{:})).^2;
            
            % Constants
            g0=1/3; g1=2/3; epsilon=1e-6;
            
            % weigths
            wt0 = g0./(epsilon+b0).^2;
            wt1 = g1./(epsilon+b1).^2;
            sum_wt=wt0+wt1;
            
            % Non-linear weigths
            w0 = wt0./sum_wt;
            w1 = wt1./sum_wt;
            
            % WENO5 polynomial
            qp = w0.*up0 + w1.*up1;
            
            % Compute flux
            %fp = qp.*a; outputf = -(fp-circshift(fp,ind));
            fp = qp.*a; outputf = -(fp-fp(m{:}));
            % End of reconstruction.
        end
        
        function outputf = WENO5_reconstruction(a,q,mm,m,p,pp)
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
                        
            % Reconstruction Polynomials
            up0 = (2*q(mm{:}) - 7*q(m{:}) + 11*q     )/6.;
            up1 = ( -q(m{:})  + 5*q       + 2*q(p{:}))/6.;
            up2 = (2*q        + 5*q(p{:}) -  q(pp{:}))/6.;
            
            % Smooth parameters
            b0 = 13/12*(q(mm{:})-2*q(m{:})+q ).^2       + 1/4*(q(mm{:})-4*q(m{:})+3*q).^2;
            b1 = 13/12*(q(m{:}) -2*q+q(p{:}) ).^2       + 1/4*(q(m{:})-q(p{:})).^2;
            b2 = 13/12*(q       -2*q(p{:})+q(pp{:})).^2 + 1/4*(3*q-4*q(p{:})+q(pp{:})).^2;
            
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
            %fp = qp.*a; outputf = -(fp-circshift(fp,ind));
            fp = qp.*a; outputf = -(fp-fp(m{:}));
            % End of reconstruction.
        end
        
        function outputf = WENO7_reconstruction(a,q,mmm,mm,m,p,pp,ppp)
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
            
            % Reconstruction Polynomials
            p0 = (-3*q(mmm{:}) + 13*q(mm{:}) - 23*q(m{:})  + 25*q        )/12;
            p1 = ( 1*q(mm{:})  - 5*q(m{:})   + 13*q        +  3*q(p{:})  )/12;
            p2 = (-1*q(m{:})   + 7*q         +  7*q(p{:})  -  1*q(pp{:}) )/12;
            p3 = ( 3*q        + 13*q(p{:})   -  5*q(pp{:}) +  1*q(ppp{:}))/12;
            
            % Smooth Indicators
            B0 = q(m{:}).*(134241*q(m{:})-114894*q) ...
                +q(mmm{:}).*(56694*q(m{:})-47214*q(mm{:})+6649*q(mmm{:})-22778*q)...
                +25729*q.^2  +q(mm{:}).*(-210282*q(m{:})+85641*q(mm{:})+86214*q);
            B1 = q.*(41001*q-30414*q(p{:})) ...
                +q(mm{:}).*(-19374*q(m{:})+3169*q(mm{:})+19014*q-5978*q(p{:}))...
                +6649*q(p{:}).^2   +q(m{:}).*(33441*q(m{:})-70602*q+23094*q(p{:}));
            B2 = q(p{:}).*(33441*q(p{:})-19374*q(pp{:})) ...
                +q(m{:}).*(6649*q(m{:})-30414*q+23094*q(p{:})-5978*q(pp{:}))...
                +3169*q(pp{:}).^2  +q.*(41001*q-70602*q(p{:})+19014*q(pp{:}));
            B3 = q(pp{:}).*(85641*q(pp{:})-47214*q(ppp{:})) ...
                +q.*(25729*q-114894*q(p{:})+86214*q(pp{:})-22778*q(ppp{:}))...
                +6649*q(ppp{:}).^2 +q(p{:}).*(134241*q(p{:})-210282*q(pp{:})+56694*q(ppp{:}));
            
            % Constants
            g0 = 1/35; g1 = 12/35; g2 = 18/35; g3 = 4/35; epsilon = 1e-6;
            
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
            %fp = qp.*a; outputf = -(fp-circshift(fp,ind));
            fp = qp.*a; outputf = -(fp-fp(m{:}));
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
        
    end % Advection/Reconstruction Subroutines
end % class