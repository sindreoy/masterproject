function [ Vl, rhoc, rhod, sigma, Vmax, P, nu ] = setParams( Rmax, flag )
% Set all parameters needed for the program to function

Vl    = 725e-6;      % Volume of liquid in tank           [m^3]
rhoc  = 1e3;         % Density continuous phase (water)   [kg/m^3]
%rhod  = 0.837e3;    % Density of dispersed phase (oil)   [kg/m^3]
%nu    = 16.88e-3;   % Kinematic viscosity                [m^2/s]
%sigma = 22e-3;      % Surface tension                    [N/m]
%Rmax  = 120e-6;      % Maximum allowed radius of bubbles  [m]
Vmax  = 4/3*Rmax^3;  % Maximum volume                     [m^3]
%P     = 0.366;      % Power usage                        [W]

switch flag
    case 1
        % We chose crude oil B
        rhod  = 0.837e3;
        nu    = 16.88e-3;
        sigma = 22e-3;
        P     = 0.366;
    case 2
        % We chose crude oil C
        rhod  = 0.911e3;
        nu    = 81.67e-3;
        sigma = 19e-3;
        P     = 0.152;
    otherwise
        ME = MException('MATLAB:IllegalFlag', ...
            ['You have passed in an illegal argument. '...
            'Flag must be either 1 or 2. Received %i.\n'], flag);
        throw(ME);
% Crude B:
%   rhod  = 0.837e3
%   nu    = 16.88e-3
%   sigma = 22e-3
%   P     = 0.366
% Crude C:
%   rhod  = 0.911e3
%   nu    = 81.67e-3
%   sigma = 19e-3
%   P     = 0.152
% 
% eps = P / (rhod * Vl)
end

