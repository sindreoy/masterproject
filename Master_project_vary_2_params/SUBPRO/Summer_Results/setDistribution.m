% Title: Set Distribution
% Author: Sindre Bakke Oyen
% Date (started): 07.06.2017
% Description: Sets a normal distribution with mean mu
%              and standard deviation sigma. The program plots the
%              distribution and writes it to a tab separated textfile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Start timer to find runtime
tic

%% Set initials and create normal distribution
mu = 200e-6;
sigma = 50e-6;

x = 0:10e-6:400e-6;
f = 1/sqrt(2*pi*sigma^2) * exp(-(x-mu).^2/(2*sigma^2));

%% Plot and check conservations
hFig = createFigure(1);
fontProps.FontName = 'Calibri';
fontProps.FontSize = 14;
fontProps.FontWeight = 'bold';
hAxes = axes('Xscale', 'linear');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on')

title('Normal Distribution')
xlabel('radius, R [m]')
ylabel('Number Density, f [-]')

plot(x, f, 'Color', 'r', 'LineWidth', 2, 'Marker',...
    'o', 'MarkerEdgeColor','r', 'MarkerFaceColor', 'none',...
    'DisplayName', 'f(r)')

legend('\sigma = 50, \mu = 200')
hLegend = legend(hAxes, 'show');
set(hLegend, fontProps, 'Location', 'NorthEastOutside');


saveas(hFig, 'intial_distro', 'epsc') 

%% Write to normal distribution with radii to file
result_matrix = [x; f];
fid = fopen('raw_distribution.txt', 'w');
fprintf(fid, '%1s %10s \n', 'r', 'f_0');
fprintf(fid, '%1.6f %10.6f \n', result_matrix);
fclose(fid);

t_elapsed = toc
% END OF PROGRAM