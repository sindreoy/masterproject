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
% clc

% Global parameters are used for things that will never change.
% I have not used global for the constants to make the program more robust
% in case dependencies are added, making them functions of temperature etc
global Rmax w xi N xipBB xipBC xippBC xipDC tf

%% Setting parameters and constants
[ Vl, rhoc, rhod, nu, sigma, Vmax ] = setParams();
eps = 1;
phi = 0.05;

tau = 0:1e-3:1;

kb1 = 1e-4;          % Model fitted parameter 1, g        [-]
                     % Magnitude of breakage source
kb2 = 1e-1;             % Model fitted parameter 2, g        [-]
                     % Increase to increase width of breakage source
kc1 = 1e-4;          % Model fitted parameter, probability[-]
                     % Magnitude of coalescence source
kc2 = 5e1;           % Model fitted parameter, efficiency [-]
                     % Increase to shorten width of coalescence source

%% Initializing and discretizing
% Setting ODE options
options = odeset();

[xi, A, B, w] = Collocation(188,1,1);
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

% Fetch initial distribution and interpolate onto xi domain
f0 = 'raw_logNormal3.txt';
[fn, fv, r] = rescaleInitial(f0, phi, 0);
psi0 = pchip(r/Rmax, fv*Rmax, xi);

%% Initialize some variables needed for program
fontProps.FontName = 'Calibri';
fontProps.FontSize = 14;
fontProps.FontWeight = 'bold';

% Set plot properties
colorMatrix{1} = sprintf('b');
colorMatrix{2} = sprintf('r');
colorMatrix{3} = sprintf('g');

% Store all evaluations
means    = [];
psi_cell = cell(1, 1);

legendArray = cell(1,1); % To be used in plotting section
colorMatrix = cell(1, 1);
k = 1;

% kb1, kc1 varying
kb1s = [kb1, kb1*2, kb1*3];
kc1s = [kc1, kc1*2, kc1*3];

% kb2, kc2 varying
kb2s = [kb2, kb2*2, kb2*3];
kc2s = [kc2, kc2*2, kc2*3];

fprintf('1) Vary dynamics (kb1, kc1)\n')
fprintf('2) Vary steady state distribution (kb2, kc2)\n')
choice = input(...
    'Choose whether you want to vary dynamics or steady state distribution: ');
for i = 1:3
    switch choice
        case 1 % Vary dyanamics
        kb1 = kb1s(i);
        kc1 = kc1s(i);
        legendArray{i} = sprintf('k_{b1}=%1.6f, k_{c1}=%1.6f', kb1, kc1);
        case 2 % Vary steady state distribution
            kb2 = kb2s(i);
            kc2 = kc2s(i);
            legendArray{i} = sprintf('k_{b2}=%1.2f, k_{c2}=%1.2f', kb2, kc2);
    end %switch choice
    % Final constants
    k1 = tf*kb1*eps^(1/3)/(2^(2/3)*Rmax^(2/3))*sqrt(rhod/rhoc);
    k2 = kb2*sigma/(rhod*2^(5/3)*eps^(2/3)*Rmax^(5/3));
    k3 = tf/Vmax * Rmax^(7/3)*4*2^(1/3)*kc1*eps^(1/3);
    k4 = kc2*Rmax^(5/6)*rhoc^(1/2)*eps^(1/3)/(2*sigma^(1/2));

    tic
    [kern.BB.k, kern.DB.k, kern.BC.k, kern.DC.k] =...
        evalKernels(k2, k4, 2);
    t_kern = toc;

    % Used in solving later
    const.k1 = k1;
    const.k3 = k3;
    const.phi = phi;

    % Set plot properties
    colorMatrix{1} = sprintf('b');
    colorMatrix{2} = sprintf('r');
    colorMatrix{3} = sprintf('g');
    colorMatrix{4} = sprintf('m');

    tic
    [tau, psi] = ...
        ode15s(@evalSource, tau, psi0, options, kern, const, 2);
    t_ode = toc
    % Store time evolution of psi for current distribution
    psi_cell{k}  = psi;
    means        = [means, getCurrentMean(psi, phi)];
    % Next initial distribution
    k = k + 1;
end %for
ax = plotMean(means, tau, fontProps, colorMatrix);
legend(ax, legendArray)
% Save variables for later plotting
now  = strsplit(char(datetime())); % Cell of date and time
time = strsplit(now{2}, ':');
now  = strcat('Results/Vary2params/', ...
    now{1}, '-', time{1}, '_', time{2}, '_', time{3});
save(now)

% END PROGRAM
