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

clear all
% close all
clc

rng default;
% Global parameters are used for things that will never change.
% I have not used global for the constants to make the program more robust
% in case dependencies are added, making them functions of temperature etc

%% Get initial distribution and discretize
phi = 0.7e-1;
f0 = 'Experimental/august/crudeB.csv';
[fn, fv, r, t] = rescaleInitial(f0, phi, 0);
Rmax = 120e-6;

[xi, A, B, w] = Collocation(98,1,1);
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
%% Setting parameters and constants
kb1 = 1.5e-6;        % Model fitted parameter 1, g        [-]
                     % Dynamics, breakage
kb2 = 1e-2;          % Model fitted parameter 2, g        [-]
                     % Steady state settlement, breakage
kc1 = 1.5e-5;        % Model fitted parameter, probability[-]
                     % Dynamics, coalescence
kc2 = 5e2;           % Model fitted parameter, efficiency [-]
                     % Steady state settlement, coalescence
ratio = kb1/kc1;

kb1 = 1e3;
kb2 = 1e3;
kc1 = 1e3;
kc2 = 1e3;

consts.Rmax       = Rmax;
consts.xis.xi     = xi;
consts.xis.xipBB  = xipBB;
consts.xis.xipBC  = xipBC;
consts.xis.xippBC = xippBC;
consts.xis.xipDC  = xipDC;
consts.fv         = fv;
consts.t          = t;
consts.tf         = t(end);
consts.r          = r;
consts.Rmax       = Rmax;
consts.w          = w;

%% What search area to chart?
Niter = 60;     % Number of iterations (Niter x Niter SSE matrix produced)
stepSize = 1.4; % The parameters charted will be multiplied by this
fprintf('1 : kb1, kb2\n')
fprintf('2 : kb1, kc1\n')
fprintf('3 : kb1, kc2\n')
fprintf('4 : kb2, kc1\n')
fprintf('5 : kb2, kc2\n')
fprintf('6 : kc1, kc2\n')
flag = input('Which parameters would you explore? ');
steps = zeros(Niter, 1);
steps(1) = 1;
for i = 1:Niter-1
    steps(i+1) = steps(i) / stepSize;
end %for
switch flag
    case 1
        k1_vec = kb1 * steps;
        k2_vec = kb2 * steps;
        p1   = kc1;
        p2   = kc2;
    case 2
        k1_vec = kb1 * steps;
        p1   = kb2;
        k2_vec = kc1 * steps;
        p2   = kc2;
    case 3
        k1_vec = kb1 * steps;
        p1   = kb2;
        p2   = kc1;
        k2_vec = kc2 * steps;
    case 4
        p1   = kb1;
        k1_vec = kb2 * steps;
        k2_vec = kc1 * steps;
        p2   = kc2;
    case 5
        p1   = kb1;
        k1_vec = kb2 * steps;
        p2   = kc1;
        k2_vec = kc2 * steps;
    case 6
        p1   = kb1;
        p2   = kb2;
        k1_vec = kc1 * steps;
        k2_vec = kc2 * steps;
    otherwise
        error('Flag does not match any of the given\n');
end %switch

%% Solve and chart sensitivity
eSquared = zeros(Niter, Niter);
parfor i = 1:Niter
    k1 = k1_vec(i);
    tmpSquared = zeros(Niter, 1);
    for j = 1:Niter
        k2 = k2_vec(j);
        tmpSquared(j) = getSSE(k1, k2, p1, p2, consts, flag);
    end %for
    eSquared(i, :) = tmpSquared;
end %parfor

createFigure();
ax = axes();
hold(ax, 'on');
set(ax, 'Xscale', 'log');
set(ax, 'Yscale', 'log');
set(ax, 'Zscale', 'log');
surf(k1_vec, k2_vec, eSquared)
xlabel('k_{b,1}')
ylabel('k_{b,2}')
zlabel('SSE')
view(3);
%% Set plot properties
greekeps = char(949); % Greek letter epsilon
greekphi = char(966);     % Greek letter phi

colorMatrix{1} = sprintf('b');
colorMatrix{2} = sprintf('r');
colorMatrix{3} = sprintf('g');
colorMatrix{4} = sprintf('m');

fontProps.FontName = 'Calibri';
fontProps.FontSize = 14;
fontProps.FontWeight = 'bold';

%% Save variables for later plotting
now  = strsplit(char(datetime())); % Cell of date and time
time = strsplit(now{2}, ':');
now  = strcat('Results/Experimental/crudeB/', ...
     now{1}, '-', time{1}, '_', time{2}, '_', time{3});
save(now)

% END PROGRAM
