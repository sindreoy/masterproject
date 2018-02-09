function [ hFig ] = plotMean( mean, tau, fontProps, colorMatrix )
% Title: Plot mean
% Author: Sindre Bakke Oyen
% Date (started): 08.09.2017
% Description: Plot the development of the mean for several distributions.
%
% Output args:
%       hFig (object)      :: created figure with plots in it
% Input args:
%       mean (array)       :: Array of means as they develop with time for
%                             several distributions
%       fontProps (object) :: Properties for fonts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

greekmu = char(956);

hFig = createFigure();
hAxes = axes('Xscale', 'log');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on');

title('Development of mean with time')
xlabel('Time [-]')
ylabel(sprintf('Mean [%sm]', greekmu))

% xlim([5e-2 1])

cols = size(mean, 2);

for i=1:cols
    plot(tau, mean(:, i)/1e-6, 'LineWidth', 2, 'Color', colorMatrix{i})
end %for

end %function

