function [ hFig ] = plot2D( fv, fvmax, t_struct, psi, fontProps )

global Rmax xi

greektau = char(964);

t = t_struct.t;
T = t_struct.T;
tau = t_struct.tau;

hFig = createFigure();
hAxes = axes('Xscale', 'linear');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on');

title('Volume Density Distribution as Time Increases')
xlabel('Radius [m]')
ylabel('Volume Density Distribution [1/m]')

xlim([0 Rmax]);
ylim([0 fvmax + fvmax/10]);

legendArray = cell(T, 1);
for i = 1:T % loop over time
    plot(xi*Rmax, fv(i,:), 'LineWidth', 2)
    legendArray{i} = sprintf('%s=%1.2f', greektau, tau(i));
    hold on
end %for
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');
save2 = lower(input('Save 2D plot ("y" or "n"): ', 's'));
if strcmp(save2, 'y')
    plot_name_handle = input('Enter filename: ', 's');
    saveas(hFig, plot_name_handle, 'epsc')
end %if
end %function
