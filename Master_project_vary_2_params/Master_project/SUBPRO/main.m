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

%% Setting parameters and constants
[ Vl, rhoc, rhod, nu, sigma, Vmax ] = setParams();
epsmat = [1, 2, 4];
phimat = [0.05, 0.15, 0.3];
tau = 0:1e-3:1;

kb1 = 1.5e-2;          % Model fitted parameter 1, g        [-]
                     % Magnitude of breakage source
kb2 = 9e-2;          % Model fitted parameter 2, g        [-]
                     % Increase to increase width of breakage source
kc1 = 1.5e-2;          % Model fitted parameter, probability[-]
                     % Magnitude of coalescence source
kc2 = 2.5e1;           % Model fitted parameter, efficiency [-]
                     % Increase to shorten width of coalescence source

%% Initializing and discretizing
[xi, A, B, w] = Collocation(128,1,1);
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

%% Initialize some variables needed for program
greekeps = char(949); % Greek letter epsilon
greekphi = char(966);     % Greek letter phi

fontProps.FontName = 'Calibri';
fontProps.FontSize = 14;
fontProps.FontWeight = 'bold';

fprintf('1) Equilibriate 4 initial distributions\n')
fprintf('2) Keep epsilon or phi constant and check final distribution\n')
cont = input('Which case do you wish to initiate (press 0 to exit)? ');
%% Menu based program
while (cont ~= 0) % Whole program loops
    k = 1;
    psi_end = [];
    psi0_mat = [];
    mean = [];              % First moment of distribution
    psi_cell = cell(1, 1);
    
    legendArray = cell(1,1); % To be used in plotting section
    colorMatrix = cell(1, 1);
switch cont
    case 1
        %% We have chosen to run 4 initial distributions
        phi = phimat(1);
        eps = epsmat(1);
        
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
        
        while k < 5
            switch k
                case 1
                    % mean = 350-6, s.d. = 60-6
                    legendArray{k} = sprintf...
                        ('m=350e-6, s.d.=60e-6');
                    logNorm3 = 'raw_logNormal3.txt';
                    [fn, fv, r] = rescaleInitial(logNorm3, phi, 0);
                case 2
                    % m = 400e-6, s.d. = 70e-6
                    legendArray{k} = sprintf...
                        ('m=400e-6, s.d.=70e-6');
                    logNorm4 = 'raw_logNormal4.txt';
                    [fn, fv, r] = rescaleInitial(logNorm4, phi, 0);
                case 3
                    % m = 450e-6, s.d. = 50e-6
                    legendArray{k} = sprintf...
                        ('m=450e-6, s.d.=50e-6');
                    logNorm5 = 'raw_logNormal5.txt';
                    [fn, fv, r] = rescaleInitial(logNorm5, phi, 0);
                case 4
                    % m = 200e-6, s.d. = 45e-6
                    legendArray{k} = sprintf...
                        ('m=200e-6, s.d.=45e-6');
                    logNorm6 = 'raw_logNormal6.txt';
                    [fn, fv, r] = rescaleInitial(logNorm6, phi, 0);
                otherwise
                    break;
            end %switch
            
            [psi0_mat, psi_end, mean, psi] = getDistributionAndMean...
                (r, fv, phi, psi0_mat, psi_end, mean, kern, const, tau);
            
            % Store time evolution of psi for current distribution
            psi_cell{k} = psi;
            
            % Next initial distribution
            k = k + 1;
        end %while initial distributions
        % Plotting
            psimax   = max( max(max(psi0_mat)), max(max(psi_end)) );
           
            psi_struct.psi0_mat = psi0_mat;
            psi_struct.psi_end = psi_end;
            psi_struct.psimax = psimax;
            
            legendArray{end+1} = sprintf('Equilibrated distributions');
            
            % Save variables for later plotting
            now = strsplit(char(datetime())); % Cell of date and time
            now = strcat('Results/equilibrate/', now{1}, '-', now{2});
            save(now)

    case 2
        %% Single distribution, holding eps or phi constant
        % Choose a distribution
        logNorm3 = 'raw_logNormal3.txt';
        
        fprintf('\n1) Hold %s constant\n', greekeps)
        fprintf('2) Hold %s constant\n', greekphi)
        epsOrPhiConst = ...
            input('Which parameter do you wish to hold constant? ');
        switch epsOrPhiConst
            case 1
                % Hold epsilon constant
                fprintf('\n1) %s=1\n', greekeps)
                fprintf('2) %s=2\n', greekeps)
                fprintf('3) %s=4\n', greekeps)
                epsChoice = input('Choose the constant value of epsilon: ');
                switch epsChoice
                    case 1
                        eps = 1;
                    case 2 
                        eps = 2;
                    case 3
                        eps = 4;
                    otherwise
                        fprintf('Illegal value %i. Only 1, 2 or 3 is accepted\n', epsChoice)
                end %switch choose eps
            case 2
                % Hold phi constant
                fprintf('\n1) %s=0.05\n', greekphi)
                fprintf('2) %s=0.15\n', greekphi)
                fprintf('3) %s=0.30\n', greekphi)
                phiChoice = input('Choose the constant value of phi: ');
                switch phiChoice
                    case 1
                        phi = 0.05;
                    case 2 
                        phi = 0.15;
                    case 3
                        phi = 0.30;
                    otherwise
                        fprintf('Illegal value %i. Only 1, 2 or 3 is accepted\n', phiChoice)
                end %switch choose eps
            otherwise
                fprintf('Illegal value %i. Only 1 or 2 is accepted\n', epsOrPhiConst)
        end %switch eps or phi constant
        
        % We have covered what to hold constant at this point and are ready
        % to loop through all 3 cases.
        while k < 4
            switch epsOrPhiConst
                case 1
                    % eps is constant
                    phi = phimat(k);
                    legendArray{k} = sprintf('%s=const=%1.2f, %s=%1.2f',...
                        greekeps, eps, greekphi, phi);
                case 2
                    % phi is constant
                    eps = epsmat(k);
                    legendArray{k} = sprintf('%s=const=%1.2f, %s=%1.2f',...
                        greekphi, phi, greekeps, eps);
            end %switch epsorphi
            % Start calculations
            [fn, fv, r] = rescaleInitial(logNorm3, phi, 0);
            
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
            
            [psi0_mat, psi_end, mean, psi] = getDistributionAndMean...
                (r, fv, phi, psi0_mat, psi_end, mean, kern, const, tau);
            
            % Store time evolution of psi for current distribution
            psi_cell{k} = psi;
            
            k = k + 1;
        end %while holding something constant
        
        % Set plot properties
        colorMatrix{1} = sprintf('b');
        colorMatrix{2} = sprintf('r');
        colorMatrix{3} = sprintf('g');
        
        psi0_mat = psi0_mat(:, 1);
        psimax   = max( max(max(psi0_mat)), max(max(psi_end)) );
           
        psi_struct.psi0_mat = psi0_mat;
        psi_struct.psi_end = psi_end;
        psi_struct.psimax = psimax;
        
        legendArray{end+1} = sprintf('Initial distribution');
        
        % Save variables for later plotting
        now = strsplit(char(datetime())); % Cell of date and time
        now = strcat('Results/epsOrPhiConst/', now{1}, '-', now{2});
        save(now)
        
    otherwise
        fprintf('Illegal value %i. Only 1 or 2 is accepted\n', cont)
end %switch continue
fprintf('\n1) Equilibriate 4 initial distributions\n')
fprintf('2) Keep epsilon or phi constant and check final distribution\n')
cont = input('Which case do you wish to initiate (press 0 to exit)? ');
end %continue loop

% END PROGRAM
