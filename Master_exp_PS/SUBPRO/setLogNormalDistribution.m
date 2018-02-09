% Title: Set Distribution
% Author: Sindre Bakke Oyen
% Date (started): 07.06.2017
% Description: Sets a log normal distribution with mean mu
%              and standard deviation sigma. The program plots the
%              distribution and writes it to a tab separated textfile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Start timer to find runtime
tic

%% Set initials and create normal distribution
iterations = 0;
legendArray = cell(4,1);
x_mat = [];
f_mat = [];
while iterations <= 4
    iterations = iterations + 1;
switch iterations
    case 1
        mean  = 350e-6;
        sd    = 60e-6;
        fname = 'raw_logNormal3.txt';
        legendArray{iterations} = sprintf('m=350e-6, s.d.=60e-6');
    case 2
        mean  = 400e-6;
        sd    = 70e-6;
        fname = 'raw_logNormal4.txt';
        legendArray{iterations} = sprintf('m=400e-6, s.d.=70e-6');
    case 3
        mean  = 450e-6;
        sd    = 50e-6;
        fname = 'raw_logNormal5.txt';
        legendArray{iterations} = sprintf('m=450e-6, s.d.=50e-6');
    case 4
        mean  = 200e-6;
        sd    = 45e-6;
        fname = 'raw_logNormal6.txt';
        legendArray{iterations} = sprintf('m=200e-6, s.d.=45e-6');
end
var  = sd^2;

mu = log(mean/sqrt(1+var/mean^2));
sigma = sqrt(log(1+var/mean^2));

x = 0:5e-6:1500e-6;
f = 1./(x*sigma*sqrt(2*pi)).*exp(-(log(x)-mu).^2/(2*sigma^2));
f(1) = 0;

%% Write to normal distribution with radii to file
result_matrix = [x; f];
fid = fopen(fname, 'w');
fprintf(fid, '%1s %10s \n', 'r', 'f_0');
fprintf(fid, '%1.6f %10.6f \n', result_matrix);
fclose(fid);

% Make sure results from last iterations are not left
clear result_matrix
x_mat = [x_mat; x];
f_mat = [f_mat; f];
end

%% Plot and check conservations
hFig = createFigure(1);
fontProps.FontName = 'Calibri';
fontProps.FontSize = 14;
fontProps.FontWeight = 'bold';
hAxes = axes('Xscale', 'linear');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on')

title('Log Normal Distribution')
xlabel('radius, R [m]')
ylabel('Number Density, f [-]')

for i = 1:size(f_mat, 1)
    plot(x_mat(i, :), f_mat(i, :), 'LineWidth', 2, 'Marker', 'o')
end
hLegend = legend(legendArray);
set(hLegend, fontProps, 'Location', 'NorthEastOutside');

%saveas(hFig, 'intial_logNormal2', 'epsc') 

t_elapsed = toc
% END OF PROGRAM