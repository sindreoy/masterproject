function ax = plotMean( mean, tau, fontProps, colorMatrix )
% Title: Plot mean
% Author: Sindre Bakke Oyen
% Date (started): 08.09.2017
% Description: Plot the development of the mean for several distributions.
%
% Output args:
%       axes (object)      :: created axes with plots in it
% Input args:
%       mean (array)       :: Array of means as they develop with time for
%                             several distributions
%       fontProps (object) :: Properties for fonts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global tf

greekmu = char(956);

createFigure();
ax = axes('Xscale', 'log');
set(ax, fontProps);
box(ax, 'on');
hold(ax, 'on');

title('Development of mean with time')
xlabel('Time [h]')
ylabel(sprintf('Mean [%sm]', greekmu))

xlim([tau(2)-tau(1) 1]*tf/3600)

cols = size(mean, 2);

for i=1:cols
    plot(tau*tf/3600, mean(:, i)/1e-6, 'LineWidth', 2, 'Color', colorMatrix{i})
end %for

end %function

