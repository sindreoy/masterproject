
function [ fn, fv, r ] = rescaleInitial( f0, phi, flag )
% Title: Rescale Initial
% Author: Sindre Bakke Oyen
% Date (started): 07.06.2017
% Description: This function should rescale an initial distribution
%              and return both the number distribution function
%              and its corresponding volumetric number distribution
%              function.
%
% Output args:
%       fn (array)   :: number distribution function [1/m^3*m]
%       fv (array)   :: volumetric number distribution function [1/m]
%       r  (array)   :: radial discretization
% Input args:
%       f0 (csv)     :: initial distribution with radii
%       phi (scalar) :: volume fraction of dispersed phase
%       flag (bool)  :: whether f0 is a volumetric
%                       or a normal number distribution
%                       flag == 0 means f0 is a fv and 
%                       flag == 1 means f0 is a fn
%
% This function should rescale the integral since
% I = integral (fv*dr) from rmin to rmax
% I = integral (v(r)*fn*dr) from rmin to rmax
% Rescale: fi = phi/I * f0, where i is either n or v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We shall utilize the trapezoidal rule:
% integral (f(x)dx) from a to b = h/2*sum(f(x_(k+1))+f(x_k))
% = (b-a)/2N * (f(x_1)+2f(x_2)+2f(x_3)+...+2f(x_N)+f(x_(N+1)))

%% Fetching out data from table
A = importdata(f0);
% Find index of density distribution
r = A.data(1, :) * 1e-6; % r is in micrometers
f = A.data(2:end, :);

r = r(1:80);        % After this f-values are irrelevant
f = f(:, 1:80);     % After this f-values are irrelevant

r = r / 2; % originally in diameters, now in radii

%% Setting variables and preparing for integrating
rmin = r(1);
rmax = r(end);
N = length(r) - 1;

[rows, cols] = size(f);
fv = zeros(rows, cols);
fn = zeros(rows, cols);
switch flag
    case 0
    %% The initial distribution was a fv
    for j=1:rows
        % The variable I will represent the approximated integral
%         I = f(j,1) + f(j,N+1);
%         for i=2:N
%             I = I + 2*f(j,i);
%         end %for
%         I = I * (rmax - rmin) / (2 * N);
        I = 0;
        for i = 1:N
            I = I + ( r(i+1) - r(i) ) * ( f(j, i+1) + f(j, i));
        end %for i
        I = I / 2;
%         trapz(r, f(j, :))
        fv(j, :) = phi / I * f(j,:);
        for i=1:cols
            if r(i) == 0
                fn(j, i) = 0; % There are no particles of radius = 0
            else
            fn(j, i) = fv(j, i)/(4/3*pi*r(i)^3);
            end %if
        end %for i
    end %for j
    case 1
        %% Initial distribution was fn
        % The integrand is v(r)*fn
        v = 4/3*pi*r.^3;
        integrand = v.*f;
        I = integrand(1) + integrand(end);
        for i=2:N
            I = I + 2*integrand(i);
        end %for
        I = I * (rmax - rmin) / (2 * N);
        
        fn = phi / I * f;
        fv = fn*4/3*pi.*r.^3;
        
    otherwise
        %% The received f0 was neither fv nor fn
        ME = MException('MATLAB:IllegalFlag',...
            ['You have passed in an illegal argument. '...
            'Flag must be either 0 or 1. Received %i.\n'], flag);
        throw(ME);

end % switch
end % function