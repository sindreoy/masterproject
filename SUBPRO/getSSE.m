function [eSquared] = getSSE(k_1, k_2, p1, p2, consts, flag)
% Input flag :: decides what pair of parameters received
%               1 : kb1, kb2
%               2 : kb1, kc1
%               3 : kb1, kc2
%               4 : kb2, kc1
%               5 : kb2, kc2
%               6 : kc1, kc2
Rmax = consts.Rmax;
r    = consts.r;
xi   = consts.xis.xi;
w    = consts.w;
fv   = consts.fv;
t    = consts.t;
tf   = consts.tf;
phi  = consts.phi;
%% Fetch betas and experimental values
switch flag
    case 1 % chosen kb1, kb2
        kb1 = k_1;
        kb2 = k_2;
        kc1 = p1;
        kc2 = p2;
    case 2 % chosen kb1, kc1
        kb1 = k_1;
        kb2 = p1;
        kc1 = k_2;
        kc2 = p2;
    case 3 % chosen kb1, kc2
        kb1 = k_1;
        kb2 = p1;
        kc1 = p2;
        kc2 = k_2;
    case 4 % chosen kb2, kc1
        kb1 = p1;
        kb2 = k_1;
        kc1 = k_2;
        kc2 = p2;
    case 5 % chosen kb2, kc2
        kb1 = p1;
        kb2 = k_1;
        kc1 = p2;
        kc2 = k_2;
    case 6 % chosen kc1, kc2
        kb1 = p1;
        kb2 = p2;
        kc1 = k_1;
        kc2 = k_2;
    otherwise
        error('Flag does not match any of the given\n');
end %switch
fv0 = fv(1, :);

%% Set the parameters and constants needed
[ Vl, rhoc, rhod, sigma, Vmax, P, ~ ] = setParams( Rmax, 1 ); % Crude B
eps = P / (rhod * Vl);

% Final constants
k1 = tf*kb1*eps^(1/3)/(2^(2/3)*Rmax^(2/3))*sqrt(rhod/rhoc);
k2 = kb2*sigma/(rhod*2^(5/3)*eps^(2/3)*Rmax^(5/3));
k3 = tf/Vmax * Rmax^(7/3)*4*2^(1/3)*kc1*eps^(1/3);
k4 = kc2*Rmax^(5/6)*rhoc^(1/2)*eps^(1/3)/(2*sigma^(1/2));

%% Solve program
% Find kernels
[kern.BB.k, kern.DB.k, kern.BC.k, kern.DC.k] =...
    evalKernels(k2, k4, consts.xis, 2);

% Store some constants needed for the source evaluation
const.k1 = k1;
const.k3 = k3;
const.w = w;

tau  = t / tf;
psi0 = pchip(r/Rmax, fv0*Rmax, xi);

% Setting ODE options
options = odeset();
tic
[~, psi] = ...
    ode15s(@evalSource, tau, psi0, options, kern, const, consts.xis, 2);
t_ode = toc
deviation = abs((psi(end, :)*w - phi) / phi * 100);
if deviation > 5
    fprintf('The mass is not conserved. Phase fraction deviation %4.3f\n', deviation)
end %if

if size(psi) ~= size(fv)
    eSquared = NaN;
else
    fv_modeled = psi / Rmax;
    fv_modeled = pchip(xi*Rmax, fv_modeled, r);
    eSquared   = sum(sum((fv_modeled - fv).^2));
end
end %function