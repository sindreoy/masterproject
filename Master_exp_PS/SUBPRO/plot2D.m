function plot2D( psi, tau, fontProps )
% Title: Plot 2D
% Author: Sindre Bakke Oyen
% Description: This function should take in distributions at different time
%              steps and plot the distribution in the r-psi plane at each
%              20th time step
% 
% Output args:
%
% Input args:
%       psi (cell)           :: Cell array with each initial distribution
%                               and its evolution with time
%       tau (array)          :: different dimensionless time steps
%       fontProps (struct)   :: contains all properties for fonts in plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Rmax xi
warning('OFF', 'MATLAB:legend:IgnoringExtraEntries')

greektau = char(964);
greekmu = char(956);

psimax = max(max(psi));

rmicrons = xi*Rmax / 1e-6;

T = length(tau);
twentieth = floor(T / 20);

k = 1;
legendArray = cell(round(T/twentieth), 1);

createFigure();
hAxes = axes('Xscale', 'log');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on');

title('Volume Density Distribution as Time Increases')
xlabel(sprintf('Radius [%sm]', greekmu))
ylabel('Volume Density Distribution [-]')

xlim([1e1 rmicrons(end)]);
ylim([0 psimax + psimax/10]);
for i = 1:twentieth:T % loop over time
    plot(rmicrons, psi(i,:), 'LineWidth', 2)
    legendArray{k} = sprintf('%s=%1.2f', greektau, tau(i));
    hold on
    k = k + 1;
end %for time
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');
%saveas(gcf, sprintf('Plots/2Dplots/%s', ));
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');

end %function
