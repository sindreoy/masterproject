
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
for k=1:length(A.colheaders)
    if strcmp(A.colheaders{k}, 'r')
        r = A.data(:, k);
    end %if
    if strcmp(A.colheaders{k}, 'f_0')
        f = A.data(:,k);
        break;
    end %if
end %if

%% Setting variables and preparing for integrating
rmin = r(1);
rmax = r(end);
N = size(r,1) - 1;

switch flag
    case 0
    %% The initial distribution was a fv
        % The variable I will represent the approximated integral
        I = f(1) + f(end);
        for i=2:N
            I = I + 2*f(i);
        end %for
        I = I * (rmax - rmin) / (2 * N);
        fprintf('The integral was evaluated to %2.3f.\n', I)
        fv = phi / I * f;
        fn = zeros(length(fv), 1);
        for i=1:length(fn)
            if r(i) == 0
                fn(i) = 0; % There are no particles of radius = 0
            else
            fn(i) = fv(i)/(4/3*pi*r(i)^3);
            end %if
        end %for
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
        fprintf('The integral was evaluated to %2.3f.\n', I)
        
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