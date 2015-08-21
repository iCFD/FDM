function u0 = IC(x,ICcase)
% Create vector u0 with an initial condition. 
% 3 cases are available.
%**************************************************************************
%
% 3 cases available: 
% {1} Gaussian wave
% {2} Sinusoidal wave
% {3} Hyperbolic Tangent wave
%
% Coded by Manuel Diaz 2014.04.12
%**************************************************************************

% Create the selected IC
switch ICcase
    case 1 % Gaussian wave for Convection-Diffusion solution
        xmid = (x(end)+x(1))/2;
        u0 = exp(-4*(x-xmid).^2);
               
    case 2 % Sinusoidal wave for Burges solution
        u0 = -sin(pi*x);
           
    case 3 % Tanh wave for Convection-Diffusion solution
        a = x(1); b = x(end);
        xi = (4-(-4))/(b-a)*(x - a) - 4;
        u0 = 1/2*(tanh(-4*xi)+1);
    
    otherwise
        error('case not in the list')
end