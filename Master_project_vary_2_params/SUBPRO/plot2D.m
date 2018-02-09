function ax = plot2D( fv, tau, radius, fontProps )
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

global Rmax xi tf
warning('OFF', 'MATLAB:legend:IgnoringExtraEntries')

greektau = char(964);
greekmu = char(956);

psimax = 0;
for i=1:length(fv)
    tempMax = max(max(fv{i}));
    if tempMax > psimax
        psimax = tempMax;
    end %if
end %if

rmicrons = radius / 1e-6;

T = length(tau);
twentieth = floor(T / 20);

k = 1;
legendArray = cell(round(T/twentieth), 1);


for r=1:length(fv)
    psi_plot = fv{r};
    createFigure();
    ax = axes('Xscale', 'log');
    set(ax, fontProps);
    box(ax, 'on');
    hold(ax, 'on');

    title('Volume Density Distribution as Time Increases')
    xlabel(sprintf('Radius [%sm]', greekmu))
    ylabel('Volume Density Distribution [1/m]')

    xlim([1e0 rmicrons(end)]);
    ylim([0 psimax + psimax/10]);
    for i = 1:twentieth:T % loop over time
        plot(rmicrons, psi_plot(i,:), 'LineWidth', 2)
        legendArray{k} = ...
            sprintf('t=%1.2f h', tau(i));
        hold on
        k = k + 1;
    end %for time
    hLegend = legend(legendArray);
    set(hLegend, fontProps, 'Location', 'NorthEastOutside');
end %for #distributions
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');

end %function
