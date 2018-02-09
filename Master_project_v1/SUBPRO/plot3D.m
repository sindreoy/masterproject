function [ hFig ] = plot3D( fv, fvmax, t, fontProps )

global Rmax xi tf

hFig = createFigure(1);
hAxes = axes('Xscale', 'linear');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on');

title('Volume Density Distribution')
xlabel('Radius [m]')
ylabel('Time [h]')
zlabel('Volume Density Distribution [1/m]')

xlim([0 Rmax]);
ylim([0 tf/3600]);
zlim([0 fvmax + fvmax/10]);

surf(hAxes, xi*Rmax, t, fv)
view(3);
save1 = lower(input('Save surface plot ("y" or "n"): ', 's'));
if strcmp(save1, 'y')
    surf_name_handle = input('Enter filename: ', 's');
    saveas(hFig, surf_name_handle, 'epsc')
end %if
    
end %function
