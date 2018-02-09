function [p, runtime]=CMAES(y_experimental, initP, odesys, r, fitnessfct, varargin)
% Large parts of the function is taken from Hansen, N. "The CMA Evolution
% Strategy: A Tutorial".

% The function takes the following inputs:
% y_experimental: Experimental values to compare to the ones calculated by
%   the model.
% time: Time interval over which the model will be integrated
% y0: Vector containing initial states of the system. Used by the ODE solver
% initP: Vector containing initial guesses for parameters
% odesys: System of ODEs that make the model. Needs to take parameters,
%   initial state values and timespan (the "time" input) as arguments, and
%   return the calculated states.
% fitnessfct: Fitness function, i.e. sum of squared error or similar

% Additional options:
% 'adaptiveStep': True or false, to toggle the adaptive step method, which
%   uses additional negative weights to shift the mean.
% 'negativeWeights': True or false, to toggle the use of additional
%   negative weights to adapt the covariance matrix.


% -------------------- Initialization --------------------------------
inputs = inputParser; % Instance inputParser class for varargin
inputs.addParameter('adaptiveStep',0,@(x) x == 0 || x == 1); % Toggle adaptive step
inputs.addParameter('negativeWeights',0,@(x) x == 0 || x == 1); % Toggle negative weights
inputs.addParameter('pmin',NaN,@(x) isnumeric(x) && length(x) == length(initP));
inputs.addParameter('pmax',NaN,@(x) isnumeric(x) && length(x) == length(initP));
inputs.parse(varargin{:}); 

adaptiveStep = inputs.Results.adaptiveStep;
negativeWeights = inputs.Results.negativeWeights;
pmin = inputs.Results.pmin;
pmax = inputs.Results.pmax;

% INPUT PARAMETERS

N = length(initP); % Number of dimensions/parameters to estimate
if isrow(initP)   
    xmean = initP'; % objective variables initial point
else
    xmean = initP;
end

sigma = 0.5; % coordinate wise standard deviation (step-size)
stopfitness = 1e-10; % stop if fitness < stopfitness (minimization)
stopeval = 1e3*N^2; % stop after stopeval number of function evaluations
alpha_cov = 2;
errorflag = 0;

% Strategy parameter setting: Selection
lambda = 4+floor(3*log(N)); % population size, offspring number

mu = floor(lambda/2); % lambda=12; mu=3; weights = ones(mu,1); would be (3_I,12)-ES

wTmp = log((lambda+1)/2) - log(1:lambda); %Initial (un-normalized) weights


mueff = sum(wTmp(1:mu))^2 / sum(wTmp(1:mu).^2);%Effective selection mass
muNegeff = sum(wTmp(mu+1:lambda))^2 / sum(wTmp(mu+1:lambda).^2);

c1 = alpha_cov / ((N + 1.3)^2 + mueff);   %Learning rate for rank-one 
                                            %update of C
cmu = min([1 - c1, alpha_cov * (mueff - 2 + 1/mueff) / ...
    ((N + 2)^2 + alpha_cov*mueff/2)]);     %Learning rate for rank-mu
                                            %update of C
                                            
alpha_mu = 1 + c1/cmu;                        %The alphas are multipliers
alpha_mu_eff = 1 + 2*muNegeff/(mueff + 2);    %for scaling the negative
alpha_posDef = (1 - c1 - cmu)/(N*cmu);       %weights

posWeights = (1/sum(wTmp(wTmp>0))) * wTmp(wTmp>0)';
negWeights = (min([alpha_mu, alpha_mu_eff, alpha_posDef]) / (abs(sum(wTmp(wTmp<0)))))...
            * wTmp(wTmp<=0)';


switch adaptiveStep % Whether to use negative weights for adjusting mean or not
    case 0
        nMeanWeights = mu;
        meanWeights = posWeights;
    case 1
        nMeanWeights = lambda;
        meanWeights = [posWeights; negWeights./100];
end

switch negativeWeights % Whether to use negative weights for adjusting C or not
    case 0
        nCovWeights = mu;
        covWeights = posWeights;
    case 1
        nCovWeights = lambda;
        covWeights = [posWeights; negWeights];
end

% Strategy parameter setting: Adaptation
cc = (4+mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
cs = (mueff+2)/(N+mueff+5); % t-const for cumulation for sigma control
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma

% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1); % evolution paths for C and sigma

if (pmin < pmax)
    B = eye(N);
    D = diag(sqrt(0.5*(pmax - pmin)));
else
    B = eye(N); % B defines the coordinate system
    D = eye(N); % diagonal matrix D defines the scaling
end
C = B*D*(B*D)'; % covariance matrix


eigeneval = 0; % B and D updated at counteval == 0
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2)); % expectation of ||N(0,I)|| == norm(randn(N,1))

% -------------------- Generation Loop --------------------------------
counteval = 0;
% try
while counteval < stopeval

% Generate and evaluate lambda offspring
for k=1:lambda
arz(:,k) = randn(N,1); % standard normally distributed vector
arx(:,k) = xmean + sigma * (B*D * arz(:,k)); % add mutation % Eq. 40
try % To solve ODE system
    calc_values = feval(odesys, arx(:, k), r);
    arfitness(k) = feval(fitnessfct, y_experimental, calc_values); % objective function call
catch % If system is unsolveable at the chosen parameters
    arfitness(k) = inf;
end
counteval = counteval+1;
end

% Sort by fitness and compute weighted mean into xmean
[arfitness, arindex] = sort(arfitness); % minimization
xmean = arx(:,arindex(1:nMeanWeights))*meanWeights; % recombination % Eq. 42
zmean = arz(:,arindex(1:nMeanWeights))*meanWeights; % == D^-1*B'*(xmean-xold)/sigma

% Cumulation: Update evolution paths
ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean); % Eq. 43
hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4+2/(N+1);
pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (B*D*zmean); % Eq. 45

if negativeWeights == 1 % Rescaling vector lengths associated with negative weights
    covWeights = [covWeights(covWeights>=0)', ...
        covWeights(covWeights<0)'*N/norm(B*inv(D)*B'* ...
        (B*D * arz(:,arindex(1:lambda))))^2]';
end

% Adapt covariance matrix C
C = (1-c1-cmu) * C ... % regard old matrix % Eq. 47
    + c1 * (pc*pc' ... % plus rank one update
    + (1-hsig) * cc*(2-cc) * C) ... % minor correction
    + cmu ... % plus rank mu update
    * (B*D*arz(:,arindex(1:nCovWeights))) ...
    * diag(covWeights) * (B*D*arz(:,arindex(1:nCovWeights)))';

% Adapt step-size sigma
sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1)); % Eq. 44

% Update B and D from C
if counteval - eigeneval > lambda/(c1+cmu)/N/10 % to achieve O(N^2)
eigeneval = counteval;
C=triu(C)+triu(C,1)'; % enforce symmetry
[B,D] = eig(C); % eigen decomposition, B==normalized eigenvectors
D = diag(sqrt(diag(D))); % D contains standard deviations now
end

% Break, if fitness is good enough
if arfitness(1) <= stopfitness
break;
end

% Escape flat fitness, or better terminate?
if arfitness(1) == arfitness(ceil(0.7*lambda))
sigma = sigma * exp(0.2+cs/damps);
disp('warning: flat fitness, consider reformulating the objective');
end

disp([num2str(counteval) ': ' num2str(arfitness(1))]);

if adaptiveStep == 1
    meanWeights(mu+1:lambda) = meanWeights(mu+1:lambda)*0.95;
end

end % while, end generation loop
% catch
%     errorflag = 1;
% end
% -------------------- Final Message ---------------------------------
disp([num2str(counteval) ': ' num2str(arfitness(1))]);
if errorflag == 0
    p = arx(:, arindex(1)); % Return best point of last generation.
else 
    p = NaN;
end

runtime = counteval;
end




