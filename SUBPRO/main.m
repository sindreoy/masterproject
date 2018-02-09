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

% clear all
% close all
clc

rng default;
% Global parameters are used for things that will never change.
% I have not used global for the constants to make the program more robust
% in case dependencies are added, making them functions of temperature etc
global Rmax w xi N xipBB xipBC xippBC xipDC tf fv

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
%% Get initial distribution
phi = 0.7e-1;
f0 = 'Experimental/august/crudeB.csv';
[fn, fv, r] = rescaleInitial(f0, phi, 0);

%% Discretizing
[xi, A, B, w] = Collocation(198,1,1);
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

%% Ordinary least squares (OLS) for finding optimal parameters
fv_final = fv(end, :);

% kb1 = [3.08410133782123e-06];
% kb2 = [0.000154672061824731];
% ratio = [0.391897433107201];
% kc2 = [545.570881854125];
% kb1 = [1.72919052347659e-06];
% kb2 = [0.000254404605142969];
% ratio = [2.22689077028660];
% kc2 = [751.403465356064];
beta0 = [kb1, kb2, ratio, kc2];
lb    = [0,   0,   0.1, 0];
ub    = [Inf, Inf, 100, Inf];
options = optimoptions('lsqcurvefit',...
   'Algorithm', 'trust-region-reflective',...
   'StepTolerance', 1e-6,...
   'MaxIterations', 5e3, ...
   'MaxFunctionEvaluations', 5e3, ...
   'OptimalityTolerance', 1e-12,...
   'FunctionTolerance', 1e-12);
beta = beta0;
while options.StepTolerance >= 1e-10
[beta,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = ...
    lsqcurvefit(@modelfun, beta, r, fv_final, lb, ub, options);
    options.StepTolerance = options.StepTolerance * 1e-2;
end % while

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
% save(now)

% END PROGRAM
