function plotFinalDistribution( psi_struct, fontProps, colorMatrix, flag )
global xi Rmax

greekmu = char(956);

psi_end = psi_struct.psi_end;
psi0_mat = psi_struct.psi0_mat;
psimax = psi_struct.psimax;

hFig = createFigure();
hAxes = axes('Xscale', 'log');
set(hAxes, fontProps);
box(hAxes, 'on');
hold(hAxes, 'on');

title('Final Distributions')
xlabel(sprintf('Radius [%sm]', greekmu))
ylabel('Volume Density Distribution [-]')

rmicrons = xi*Rmax/1e-6;

xlim([1e1 rmicrons(end)]);
ylim([0 psimax + psimax/10]);
if flag == 0
    % Chosen to plot several initial distributions
    for i=1:size(psi0_mat, 2)
        plot(rmicrons, psi0_mat(:, i), 'LineWidth', 2, ...
            'LineStyle', '--', 'Color', colorMatrix{i})
%         plot(rmicrons, psi_end(i, :), 'LineWidth', 2, ...
%         'Color', colorMatrix{i})
    end %for
    plot(rmicrons, psi_end(1, :), 'LineWidth', 2, ...
        'Color', 'k')
elseif flag == 1
    % Chosen to plot one distribution, several finals
    for i=1:size(psi_end, 1)
        plot(rmicrons, psi_end(i, :),...
            'LineWidth', 2, 'Color', colorMatrix{i})
    end %for
    plot(rmicrons, psi0_mat(:, 1), 'LineWidth', 2, ...
        'LineStyle', '--', 'Color', 'k')
end %if
end %function

