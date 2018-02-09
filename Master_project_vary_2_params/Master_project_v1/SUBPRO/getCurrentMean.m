function [ mean ] = getCurrentMean( fv, phi )
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

[row, col] = size(fv);
mean = zeros(row, 1);
r = xi*Rmax;

for i = 1:row
%     integrand = fv(i, :).*xi'*Rmax;
%     mean(i) = 1/phi * integrand*w;
    sum = 0;
    for j = 2:col
        dr = r(j) - r(j-1);
        rdf = r(j)*(fv(i, j-1) + fv(i, j));
        sum = sum + rdf/2 * dr;
    end %for j
    mean(i) = 1/phi * sum;
end %for i

end %function