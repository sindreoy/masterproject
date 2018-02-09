function [ kBB, kDB, kBC, kDC ] = evalKernels(k2, k4, flag)
% Title: evalKernels
% Author: Sindre Bakke Oyen
% Date (started): 14.06.2017
% Description: This function evaluates all breakage and coalescence kernels
%              associated with them. The kernels for birth and death are
%              evaluated individually as they are over different domains.
%
% Output args:
%       kBB  (2D array):: Birth breakage rate of breakage
%       kDB  (1D array):: Death breakage rate of breakage
%       kBC  (2D array):: Birth coalescence rate of aggregation
%       kDC  (2D array):: Death coalescence rate of aggregation
%
% Input args:
%       kg2  (scalar)  :: Parameter for brekage frequency
%       k4   (scalar)  :: Parameter for coalescence exponent
%       flag (scalar)  :: contains information about which algorithm to run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global xi N xipBB xipBC xippBC xipDC

switch flag
    case 0 % Double for loop
        % Preallocate space for kernels
        kBB = zeros(N, N);
        kDB = zeros(N, 1);
        kBC = zeros(N, N);
        kDC = zeros(N, N);
        for row=2:N
            % DB
            kDB(row) = 1/xi(row)^(2/3) * exp(-k2/xi(row)^(5/3));
            for col=1:N
                % BB
                kBB(row, col) = 2 * 1/xipBB(row, col)^(2/3) ...
                    *exp(-k2/xipBB(row, col)^(5/3)) ...
                    *2.4/xipBB(row, col)^3 ...
                    *exp(-4.5 * (2*xi(row)^3-xipBB(row, col)^3)^2 ...
                    /xipBB(row, col)^6)...
                    *3*xi(row)^2;

                % BC
                kBC(row, col) = (xipBC(row,col)+xippBC(row,col))^2 ...
                    *(xipBC(row,col)^(2/3)+xippBC(row,col)^(2/3))^(1/2) ...
                    *exp(-k4*(1/xipBC(row,col) + 1/xippBC(row,col)) ...
                    ^(-5/6));

                % DC
                kDC(row,col) = (xipDC(col)+xi(row))^2 ...
                    *(xipDC(col)^(2/3)+xi(row)^(2/3))^(1/2) ...
                    *exp(-k4*(1/xipDC(col) + 1/xi(row))^(-5/6));
            end %col
        end %row
        
    case 1 % Single for loop
        % Preallocate space for kernels
        kBB = zeros(N, N);
        kBC = zeros(N, N);
        kDC = zeros(N, N);
        for row = 2:N
        kBB(row, :) = 2 * 1./xipBB(row,:).^(2/3) ...
                .*exp(-k2./xipBB(row,:).^(5/3)) ...
                .*2.4./xipBB(row,:).^3 ...
                .*exp(-4.5*(2.*xi(row)^3-xipBB(row,:).^3).^2 ...
                ./xipBB(row,:).^6) ...
                .*3.*xi(row)^2;

        kBC(row,:) = (xipBC(row, :)+xippBC(row, :)).^2 ...
                .*(xipBC(row, :).^(2/3)+xippBC(row, :).^(2/3)).^(1/2) ...
                .*exp(-k4*(1./xipBC(row, :) + 1./xippBC(row, :)).^(-5/6));

        kDC(row, :) = (xipDC+xi(row)).^2 ...
                .*(xipDC.^(2/3)+xi(row).^(2/3)).^(1/2) ...
                .*exp(-k4*(1./xipDC + 1/xi(row)).^(-5/6));
        end %row
        kDB = 1./xi.^(2/3) .* exp(-k2./xi.^(5/3));
        
    case 2 % No for loops
        xir = repmat(xi, 1, N);
        xiprDC = repmat(xipDC, 1, N)';
        
        kBB = 2 * 1./xipBB.^(2/3) ...
                .*exp(-k2./xipBB.^(5/3)) ...
                .*2.4./xipBB.^3 ...
                .*exp(-4.5*(2.*xir.^3-xipBB.^3).^2./xipBB.^6) ...
                .*3.*xir.^2;

        kBC = (xipBC+xippBC).^2 ...
                .*(xipBC.^(2/3)+xippBC.^(2/3)).^(1/2) ...
                .*exp(-k4*(1./xipBC + 1./xippBC).^(-5/6));

        kDC = (xiprDC+xir).^2 ...
                .*(xiprDC.^(2/3)+xir.^(2/3)).^(1/2) ...
                .*exp(-k4*(1./xiprDC + 1./xir).^(-5/6));
        kDB = 1./xi.^(2/3) .* exp(-k2./xi.^(5/3));
        
    otherwise
        ME = MException('MATLAB:IllegalFlag', ...
            ['The flag received is illegal. '...
            'Supported: 0, 1 or 2. Received: %i'], flag);
        throw(ME);
        
end %function