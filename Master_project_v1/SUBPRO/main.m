% Title: Solution of the transient nondimensionalized PBE
% Author: Sindre Bakke Oyen
% Date (started): 13.06.2017
% Description: Main script for solving the transient nondimensionalized 
%              PBE. It rescales initial distribution, f*, to fv and sets
%              the collocation points at the roots of Jacobi polynomials.
%              The points are orthogonally collocated in xi, xi' and xi''.
%
% Notation:
%   BB  :: Birth brekage
%   DB  :: Death brekage
%   BC  :: Birth coalescence
%   DC  :: Death coalescence
%   vp  :: Generic variable v prime (v')
%   vpp :: Generic variable v double prime (v'')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clc

% Global parameters are used for things that will never change.
% I have not used global for the constants to make the program more robust
% in case dependencies are added, making them functions of temperature etc
global Rmax w xi N xipBB xipBC xippBC xipDC tf


%% Initialize some variables needed for program
cont = '';
iterations = 1;
legendArray = cell(1,1); % To be used in plotting section
greekmu = char(956);    % Greek letter mu
greeksigma = char(963); % Greek letter sigma
greekepsilon = char(949); % Greek letter epsilon
greekphi = char(966);   % Greek letter phi

fv_init_mat = [];
fv_end = [];
psi_end = [];
psi0_mat = [];
mean = [];              % First moment of distribution

breakflag_outer = 0;
breakflag_eps = 0;
breakflag_phi = 0;

%% Menu based program
while (true) % Whole program loops
if iterations > 4
    breakflag_outer = 1;
    break;
end %if
if iterations == 1
    while ~(breakflag_eps)
        fprintf('1) Use %s=1\n', greekepsilon)
        fprintf('2) Use %s=2\n', greekepsilon)
        fprintf('3) Use %s=4\n', greekepsilon)
        P_choice = input('Choose an option (press 0 to exit): ');
        switch P_choice
            case 0
                breakflag_outer = 1;
                break;
            case 1
                P = 1e3;
                breakflag_eps = 1;
            case 2
                P = 2e3;
                breakflag_eps = 1;
            case 3
                P = 4e3;
                breakflag_eps = 1;
            otherwise
                fprintf('Key press %i is invalid\n\n', P_choice)
        end % switch eps
    end %while eps
    if (breakflag_outer)
        break;
    end %if
    while ~(breakflag_phi)
        fprintf('\n1) Use %s=0.05\n', greekphi)
        fprintf('2) Use %s=0.15\n', greekphi)
        fprintf('3) Use %s=0.30\n', greekphi)
        phi_choice = input('Choose an option (press 0 to exit): ');
        switch phi_choice
            case 0
                breakflag_outer = 1;
                break;
            case 1
                phi = 0.05;
                breakflag_phi = 1;
            case 2
                phi = 0.15;
                breakflag_phi = 1;
            case 3
                phi = 0.3;
                breakflag_phi = 1;
            otherwise
                fprintf('Key press %i is invalid\n\n', phi_choice)
        end %switch phi
    end %while phi
    if breakflag_outer
        break;
    end %if
end %if
switch iterations
    case 1
        % mean = 350-6, s.d. = 60-6
        legendArray{iterations} = sprintf...
            ('m=350e-6, s.d.=60e-6');
        logNorm3 = 'raw_logNormal3.txt';
        [fn, fv, r] = rescaleInitial(logNorm3, phi, 0);
    case 2
        % m = 400e-6, s.d. = 70e-6
        legendArray{iterations} = sprintf...
            ('m=400e-6, s.d.=70e-6');
        logNorm4 = 'raw_logNormal4.txt';
        [fn, fv, r] = rescaleInitial(logNorm4, phi, 0);
    case 3
        % m = 450e-6, s.d. = 50e-6
        legendArray{iterations} = sprintf...
            ('m=450e-6, s.d.=50e-6');
        logNorm5 = 'raw_logNormal5.txt';
        [fn, fv, r] = rescaleInitial(logNorm5, phi, 0);
    case 4
        % m = 200e-6, s.d. = 45e-6
        legendArray{iterations} = sprintf...
            ('m=200e-6, s.d.=45e-6');
        logNorm6 = 'raw_logNormal6.txt';
        [fn, fv, r] = rescaleInitial(logNorm6, phi, 0);
    otherwise
        break;
end %switch

%% Setting parameters
Vl    = 1;           % Volume of liquid in tank           [m^3]
rhoc  = 1e3;         % Density continuous phase (water)   [kg/m^3]
rhod  = 0.830e3;     % Density of dispersed phase (oil)   [kg/m^3]
eps   = P/(rhoc*Vl); % Energy dissipation rate            [m^2/s^3]
nu    = 1;           % Kinematic viscosity                [m^2/s]
sigma = 20e-3;       % Surface tension                    [N/m]
Rmax  = 1500e-6;     % Maximum allowed radius of bubbles  [m]
tf    = 500;          % Residence time of 15 seconds       [s]
Vmax  = 4/3*Rmax^3;  % Maximum volume                     [m^3]

%% Setting constants
kb1 = 1;          % Model fitted parameter 1, g        [-]
                     % Magnitude of breakage source
kb2 = 5e-1;          % Model fitted parameter 2, g        [-]
                     % Increase to increase width of breakage source
kc1 = 1;          % Model fitted parameter, probability[-]
                     % Magnitude of coalescence source
kc2 = 1e2;           % Model fitted parameter, efficiency [-]
                     % Increase to shorten width of coalescence source

% Final constants
k1 = tf*kb1*eps^(1/3)/(2^(2/3)*Rmax^(2/3))*sqrt(rhod/rhoc);
k2 = kb2*sigma/(rhod*2^(5/3)*eps^(2/3)*Rmax^(5/3));
k3 = tf/Vmax * Rmax^(7/3)*4*2^(1/3)*kc1*eps^(1/3);
k4 = kc2*Rmax^(5/6)*rhoc^(1/2)*eps^(1/3)/(2*sigma^(1/2));

%% Initializing and discretizing
[xi, A, B, w] = Collocation(148,1,1);
xi(1) = 1e-10;
N = length(xi);
alpha = xi;
gamma = xi;

% BB
xipBB = (1-xi)*gamma'+xi*ones(1, N);

% DC
xipDC  = xi;

%BC
xipBC  = 2^(-1/3)*xi*alpha';
xippBC = xi*(1-alpha.^3/2).^(1/3)';


% The matrices are created as such, storing all values of xi' for the
% current xi in rows
%         xi1' xi2' ... xiN'
%         --             --
%    xi1 |                 |
%    xi2 |                 |
%      . |                 |
%      . |                 |
%      . |                 |
%    xiN |                 |
%         --             --

%% Evaluate kernels at all xi
% The kernels are stored in a struct that has the following structure
% BB        :: Breakage birth (stores g and b for this phenomenon)
% DB        :: Breakage death (stores g for this phenomenon)
% BC        :: Coalescence birth (stores x for this phenomenon)
% DC        :: Coalescence death (stores x for this phenomenon)
%
% 1st term  :: Kernel for BB
% 2nd term  :: Kernel for DB
% 3rd term  :: Kernel for BC
% 4th term  :: Kernel for DC

% Start timer to find best algorithms
tic
[kern.BB.k, kern.DB.k, kern.BC.k, kern.DC.k] = evalKernels(k2, k4, 2);
t_kern = toc

%% Solving
% Preparing constants to be passed into ode15s
const.k1 = k1;
const.k3 = k3;
const.phi = phi;

% Create the initial psi over the correct domain
psi0 = pchip(r/Rmax, fv*Rmax, xi);

% Setting time interval to integrate over and ODE options
tau = 0:0.05:1;
options = odeset();
[tau, psi] = ode15s(@evalSource, tau, psi0, options, kern, const, 2);
t_ode = toc - t_kern

%% Dimensionalizing results
fv_init = fv;           % Keep initial distribution
fv_init_mat = [fv_init_mat fv]; % Store all initial distributions
fv = psi/Rmax;          % Volume density

psi0_mat = [psi0_mat psi0]; % Keep initial distribution (dimensionless)
fvmax = max(max(fv));   % Maximum value of fv
psimax = max(max(psi)); % Maximum value of psi

t = tau*tf / 3600;      % Time in hours
T = length(tau);        % Number of time discretizations

fontProps.FontName = 'Calibri'; % Set the plotting style
fontProps.FontSize = 14;
fontProps.FontWeight = 'bold';

t_struct.t = t;         % Create a struct of time for plots
t_struct.T = T;
t_struct.tau = tau;

%% Store final distribution for this iteration
mean = [mean getCurrentMean(fv, phi)];

fv_end = [fv_end; fv(end, :)];
psi_end = [psi_end; psi(end, :)];

%% Checking mass conservations
for i = 1:T
    conservation = phi - psi(i, :)*w;
    conservation_percent = abs(conservation/phi) * 100;
    if conservation_percent > 5
        fprintf('The mass conservation has NOT been satisfied in ')
        fprintf('the odesolver.\nThe phase difference between the ')
        fprintf('original phase fraction and the current is ')
        fprintf('%2.5f.\n\n', conservation_percent)
    else
        fprintf('The mass conservation is satisfied ')
        fprintf('to 5 percent relative error.\n')
    end %if
end %for
% Create a 3D surface plot
% hFig = plot3D(fv, fvmax, t, fontProps);

% Plot for each time step
% plot2D(fv, fvmax, t_struct, psi, fontProps);

iterations = iterations + 1;
end %while outer loop
if iterations == 1
    return
end %if

%% Plot the steady state means of different initial distributions
hfig3 = plotMean(mean, fontProps, tau);
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');

%% Plot the desired distributions
fvmax = max( max(max(fv_end)), max(max(fv_init_mat)) );
fv_struct.fv_end = fv_end;
fv_struct.init.fv_init = fv_init_mat;
fv_struct.fvmax = fvmax;
fv_struct.init.r = r;

psi_struct.psi_end = psi_end;
psi_struct.psi0_mat = psi0_mat;
psi_struct.psimax = psimax;

hFig4 = plotFinalDistribution(fv_struct, psi_struct, fontProps);
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');

% save2 = lower(input('Save steady state plot? ("y" or "n"): ', 's'));
% if strcmp(save2, 'y')
%     plot_name_handle = input('Enter filename: ', 's');
%     saveas(hFig, plot_name_handle, 'epsc')
% end %if
% END PROGRAM
