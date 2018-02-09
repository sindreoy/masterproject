function [fv_modeled] = modelfun(beta, r)

global tf Rmax xi fv w
%% Fetch betas and experimental values
kb1   = beta(1);
kb2   = beta(2);
ratio = beta(3);
kc2   = beta(4);

kc1 = kb1 / ratio;
fv0 = fv(1, :);
% Check for weird parameter values
if any(beta) < 0
    fv_modeled = ones(1, length(fv0)) * Inf;
    fprintf('Some parameter is less than 0.\n')
    return
end %if

%% Set the parameters and constants needed
[ Vl, rhoc, rhod, ~, sigma, Vmax, P ] = setParams( 1 ); % Crude B
eps = P / (rhod * Vl);
phi = 0.7e-1;

% Final constants
k1 = tf*kb1*eps^(1/3)/(2^(2/3)*Rmax^(2/3))*sqrt(rhod/rhoc);
k2 = kb2*sigma/(rhod*2^(5/3)*eps^(2/3)*Rmax^(5/3));
k3 = tf/Vmax * Rmax^(7/3)*4*2^(1/3)*kc1*eps^(1/3);
k4 = kc2*Rmax^(5/6)*rhoc^(1/2)*eps^(1/3)/(2*sigma^(1/2));

%% Solve program
% Find kernels
[kern.BB.k, kern.DB.k, kern.BC.k, kern.DC.k] =...
    evalKernels(k2, k4, 2);

% Store some constants needed for the source evaluation
const.k1 = k1;
const.k3 = k3;
const.phi = phi;

tau = 0:1e-3:1;
psi0 = pchip(r/Rmax, fv0*Rmax, xi);

% Setting ODE options
options = odeset();
tic
[~, psi] = ...
    ode15s(@evalSource, tau, psi0, options, kern, const, 2);
if size(psi, 1) ~= length(tau)
    fv_modeled = ones(1, length(fv0)) * Inf;
    fprintf('kb1=%2.10f, kb2=%2.10f\nkc1=%2.10f, kc2=%2.10f\n',kb1, kb2, kc1, kc2)
    return
end %if
t_ode = toc
fv_modeled = psi(end, :) / Rmax;
fv_modeled = pchip(xi*Rmax, fv_modeled, r);

%% Check mass conservation
% T = length(tau);
% for i = 1:T
%     conservation = phi - psi(i, :)*w;
%     conservation_percent = abs(conservation/phi) * 100;
%     if conservation_percent > 5
%         fprintf('The mass conservation has NOT been satisfied in ')
%         fprintf('the odesolver.\nThe phase difference between the ')
%         fprintf('original phase fraction and the current is ')
%         fprintf('%2.5f.\n\n', conservation_percent)
%     end %if
% end %for

end %function