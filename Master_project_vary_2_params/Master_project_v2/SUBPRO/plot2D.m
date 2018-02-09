function [ hFig ] = plot2D( psi, tau, fontProps )

global Rmax xi

greektau = char(964);
psimax = max(max(psi));

hFig = createFigure();
hAxes = axes('Xscale', 'log');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on');

title('Volume Density Distribution as Time Increases')
xlabel('Radius [m]')
ylabel('Volume Density Distribution [-]')

rmicrons = xi*Rmax / 1e-6;

xlim([1e1 rmicrons(end)]);
ylim([0 psimax + psimax/10]);

T = length(tau);
twentieth = floor(T / 20);

k = 1;
legendArray = cell(round(T/twentieth), 1);
for i = 1:twentieth:T % loop over time
    plot(rmicrons, psi(i,:), 'LineWidth', 2)
    legendArray{k} = sprintf('%s=%1.2f', greektau, tau(i));
    hold on
    k = k + 1;
end %for
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');
% save = lower(input('Save 2D plot ("y" or "n"): ', 's'));
% if strcmp(save, 'y')
%     plot_name_handle = input('Enter filename: ', 's');
%     saveas(hFig, plot_name_handle, 'png')
% end %if
end %function
