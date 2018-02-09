function [ mean ] = getCurrentMean( psi, phi )
% Title: Get Current Mean
% Author: Sindre Bakke Oyen
% Date (started): 08.09.2017
% Description: Gets the first moment (mean) for the current distribution at time t (that
%              is, it gets mean(fv(t, r))

% Output args:
%       mean (scalar) :: the first moment of the distribution
% Input args:
%       fv (array)    :: the final distribution
%       phi (scalar)  :: phase fraction

% The mean is defined by the first moment:
% mu = integral(x*f(x)dx)
% trapezoids => sum((f(k-1)+f(k))/2 * dx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global xi Rmax

[row, col] = size(psi);
mean = zeros(row, 1);

for i = 1:row
    sum = 0;
%     for j = 2:col
%         dxi = xi(j) - xi(j-1);
%         xidpsi = xi(j)*(psi(i, j-1) + psi(i, j));
%         sum = sum + xidpsi/2 * dxi;
%     end %for j
    sum = trapz(xi, xi'.*psi(i, :));
    mean(i) = 1/phi * sum;
end %for i

mean = mean*Rmax;

end %function