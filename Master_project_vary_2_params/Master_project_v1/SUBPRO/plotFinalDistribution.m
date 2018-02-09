function [ hFig ] = plotFinalDistribution( fv_struct, psi_struct, fontProps )
global xi Rmax w

greekmu = char(956);

fv = fv_struct.fv_end;
fvmax = fv_struct.fvmax;
fv_init = fv_struct.init.fv_init;
r_init = fv_struct.init.r;

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

colorMatrix = cell(size(psi0_mat, 2), 1);
colorMatrix{1} = sprintf('b');
colorMatrix{2} = sprintf('r');
colorMatrix{3} = sprintf('y');
colorMatrix{4} = sprintf('m');

for i=1:size(psi0_mat, 2)
    plot(rmicrons, psi0_mat(:, i), 'LineStyle', '--',...
        'Color', colorMatrix{i})
end %for

for i=1:size(fv, 1)
plot(rmicrons, psi_end(i, :)', 'LineWidth', 2,...
    'Color', colorMatrix{i})
end %for

end %function

