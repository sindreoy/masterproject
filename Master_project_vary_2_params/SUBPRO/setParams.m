function [ Vl, rhoc, rhod, nu, sigma, Vmax ] = setParams( )
% Set all parameters needed for the program to function

global tf Rmax
Vl    = 1;           % Volume of liquid in tank           [m^3]
rhoc  = 1e3;         % Density continuous phase (water)   [kg/m^3]
rhod  = 0.830e3;     % Density of dispersed phase (oil)   [kg/m^3]
nu    = 1;           % Kinematic viscosity                [m^2/s]
sigma = 20e-3;       % Surface tension                    [N/m]
Rmax  = 700e-6;     % Maximum allowed radius of bubbles  [m]
tf    = 5*60*60;          % Residence time of 15 seconds       [s]
Vmax  = 4/3*Rmax^3;  % Maximum volume                     [m^3]

end

