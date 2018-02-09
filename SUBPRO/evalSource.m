function RHS = evalSource( tau, psi, kern, const, xis, flag )
% Title: Evaluate Source
% Author: Sindre Bakke Oyen
% Date (started): 20.06.2017
% Description: This function should evaluate the right hand side of the
%              nondimensionalized PBE. It will evaluate it for each radial
%              discretization, meaning it will return an array of
%              length = number of discretization points.
%
% Output args:
%       RHS   (array)  :: source of bubbles of radius xi
% Input args:
%       tau   (array)  :: dimensionless time
%       psi   (array)  :: dimensionless volumetric density distribution
%       kern  (struct) :: contains all birth and death kernels
%       const (struct) :: contains all constants, k1, k2, k3, k5 and phi
%       flag  (scalar) :: contains information about which algorithm to run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xi = xis.xi;
xipBB = xis.xipBB;
xipBC = xis.xipBC;
xippBC = xis.xippBC;
xipDC = xis.xipDC;
N = length(xi);

% Fetch all kernels
kBB = kern.BB.k;
kDB = kern.DB.k;
kBC = kern.BC.k;
kDC = kern.DC.k;

% Fetch all constants
k1 = const.k1;
k3 = const.k3;
w  = const.w;

% Interpolate onto domains in terms of volume
psipBB   = pchip(xi, psi, xipBB);
psipBC  = pchip(xi, psi, xipBC);
psippBC = pchip(xi, psi, xippBC);
psipDC  = pchip(xi, psi, xipDC);

switch flag
    case 0 % Double for loop
        % Preallocate space for integrands
        IBB = zeros(N, N);
        IBC = zeros(N, N);
        IDC = zeros(N, N);

        % Preallocate space for the source terms
        BB = zeros(N, 1);
        DB = zeros(N, 1);
        BC = zeros(N, 1);
        DC = zeros(N, 1);
        for row=2:N % loop on xi
            for col=1:N % loop on xi' and xi''
                % BB
                IBB(row, col) = kBB(row,col) ...
                    *psipBB(row, col)/xipBB(row,col)^3 ...
                    *(1-xi(row));
                % BC
                IBC(row, col) = kBC(row, col) ...
                    *psipBC(row, col)/xipBC(row,col)^3 ...
                    *psippBC(row,col)/xippBC(row,col)^3 ...
                    *xi(row)^2/xippBC(row, col)^2 ...
                    *xi(row)*2^(-1/3);
                % DC
                IDC(row, col) = kDC(row, col) * psipDC(col)/xipDC(col)^3;
            end %col
            BB(row) = k1 * xi(row)^3 * IBB(row, :)*w;
            DB(row) = k1*kDB(row)*psi(row);
            BC(row) = k3*xi(row)^3 * IBC(row, :)*w;
            DC(row) = k3*psi(row) * IDC(row, :)*w;
        end %row

        B = BB - DB; % Net breakage
        C = BC - DC; % Net coalescence

        RHS = B + C;
        
    case 1 % Single for loop
        % Preallocate space for integrands
        IBB = zeros(N, N);
        IBC = zeros(N, N);
        IDC = zeros(N, N);
        for row=2:N
            % BB
            IBB(row, :) = kBB(row,:) ...
                .*psipBB(row, :)./xipBB(row,:).^3 ...
                 *(1-xi(row));
            % BC
            IBC(row, :) = kBC(row, :) ...
                .*psipBC(row, :)./xipBC(row, :).^3 ...
                .*psippBC(row, :)./xippBC(row, :).^3 ...
                 *xi(row)^2./xippBC(row, :).^2 ...
                 *xi(row)*2^(-1/3);
            % DC
            IDC(row, :) = kDC(row, :) .* (psipDC./xipDC.^3)';
        end %row
        BB = k1 * xi.^3. .* (IBB*w);
        DB = k1*kDB.*psi;
        BC = k3*xi.^3 .* (IBC*w);
        DC = k3*psi .* (IDC*w);
        
        B = BB - DB; % Net breakage
        C = BC - DC; % Net coalescence

        RHS = B + C;
        
    case 2 % No for loops
        xir = repmat(xi, 1, N); %xi without loops must have same dimensions
        
        IBB = kBB.*psipBB./xipBB.^3.*(1-xir);
        IBC = kBC.* ...
            psipBC./xipBC.^3 ...
            .*psippBC./xippBC.^3 ...
            .*xir.^2./xippBC.^2 ...
            .*xir*2^(-1/3);
        IDC = kDC .* repmat((psipDC./xipDC.^3)', N, 1);
        IDC(1, :) = 0;
        
        BB = k1 *  xi.^3 .* (IBB*w);
        DB = k1*kDB.*psi;
        BC = k3*xi.^3 .* (IBC*w);
        DC = k3*psi .* (IDC*w);
        
        B = BB - DB; % Net breakage
        C = BC - DC; % Net coalescence

        RHS = B + C;
        RHS;
    otherwise
        ME = MException('MATLAB:IllegalFlag', ...
            ['The flag received is illegal. '...
            'Supported: 0, 1 or 2. Received: %i'], flag);
        throw(ME);
        
end %function