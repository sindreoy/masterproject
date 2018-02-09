function [ psi ] = getDistributionAndMean...
    ( r, fv, phi, kern, const, tau )

global Rmax xi w

psi0 = pchip(r/Rmax, fv*Rmax, xi);

% Setting ODE options
options = odeset();
tic
[tau, psi] = ...
    ode15s(@evalSource, tau, psi0, options, kern, const, 2);
t_ode = toc

% fontProps.FontName = 'Calibri';
% fontProps.FontSize = 14;
% fontProps.FontWeight = 'bold';
% 
% plot2D(psi, tau, fontProps);

% Check that mass is conserved
T = length(tau);
for i = 1:T
    conservation = phi - psi(i, :)*w;
    conservation_percent = abs(conservation/phi) * 100;
    if conservation_percent > 5
        fprintf('The mass conservation has NOT been satisfied in ')
        fprintf('the odesolver.\nThe phase difference between the ')
        fprintf('original phase fraction and the current is ')
        fprintf('%2.5f.\n\n', conservation_percent)
%     else
%         fprintf('The mass conservation is satisfied ')
%         fprintf('to 5 percent relative error.\n')
    end %if
end %for

end

