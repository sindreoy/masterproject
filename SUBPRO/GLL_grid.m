function x = GLL_grid(N)
% Takes a number of grid points N and returns the roots of Legendre
% polynomials, x

% Three term recurrence for alpha=beta=0 (Legendre polynomials):
% aj = (2j+1)/(j+1)
% bj = 0
% cj = j/(j+1)
j = 1:N;
a = (2*j+1)./(j+1);
b = 0;
c = j./(j+1);

% Construct T matrix
superDiag = diag(1./a(1:end-1), 1);
mainDiag  = diag(-b./a, 0);
subDiag   = diag(c(2:end)./a(2:end), -1);
T = superDiag + mainDiag + subDiag;

% Find the abscissas as eigenvalues of T
x = [-1; sort(eig(T)); 1];

end %function