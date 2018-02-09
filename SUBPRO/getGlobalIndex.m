function [idx] = getGlobalIndex(N, n, i)
% Title: getGlobalIndex
% Author: Sindre Bakke Oyen
% Date (started): 07.02.2018
% Description: This function is meant to be used for a collocation grid
%              with the method of elements. That means the collocation grid
%              is split into subgrids for locally higher precision.
%
% Output args:
%       idx (scalar) :: The global index of the grid
% Input args:
%       N (array)    :: vector that states how many points are located in
%                       each subgrid
%       n (scalar)   :: states which subgrid we are looking at
%       i (scalar)   :: states which point in subgrid n we are looking at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you pass in something with length <= 0 or index, i, is greater than
% length of subgrid, N(n): throw error
if (isempty(N) && isempty(n) && isempty(i)) || i > N(n)
    error('Received illegal input arguments.')
end %if
% Previous subgrids
k = n - 1;

% Sum all points in previous subgrids
s = sum(N(1:k));

% Add index from current grid and remove intersections between subgrids
% k previous subgrids have k 
idx = s + i - k;

end

