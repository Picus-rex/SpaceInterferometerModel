function [eta_mod, IWA] = ...
    compute_IWA(PSF, theta_range, eta_opt, areas, autoplot)
%COMPUTE_IWA Compute the inner working angle given the PSF and the range of
%angluar apertures.
%
% INPUTS:
%   PSF[MxM]            2D PSF as output from compute_psf.
%   theta_range[Mx1]    Range of angular offsets from the star. [rad]
%   eta_opt[Nx1]        Optical efficiency of every aperture.
%   areas[Nx1]          Areas of every aperture. [m^2]
%   autoplot[bool]      If true, create plot.
%
% OUTPUTS:
%   eta_mod[1]          Asymptotic rms modulation efficiency of array;
%   IWA[1]              Internal working angle. [rad]
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   - Check plot for correctness of fmincon function.
%
% VERSION HISTORY:
%   2025-02-18 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using formula on page 8.
den = 0;
for i = 1:length(eta_opt)
    den = den + eta_opt(i) * areas(i);
end
eta_mod = max(PSF, [], "all") / den;

% Extract the position of the main peak and the main depth.
[i, j] = find(PSF == max(PSF, [], "all"));
[~, jmin] = find(PSF == min(PSF, [], "all"));

i = i(1);
j = j(1);
jmin = jmin(1);

% 2D horizontal representation
signal = PSF(i, :);

% Solver options
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...     
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'Display', 'off');      % Change to "iteration" for debugging

% Apply a scale to avoid problems.
fun = @(x) find_IWA(x, theta_range, signal, eta_mod);
x0 = theta_range(j) / 1e-7;
IWA = fmincon(fun, x0/2, [], [], [], [], theta_range(jmin)/1e-7, ...
    theta_range(j)/1e-7, [], options);
IWA = IWA * 1e-7;

if autoplot

    figure; hold on;
    plot(theta_range, signal, "LineWidth", 1.5, "DisplayName", "PSF along maximum peak")
    yline(eta_mod, '--', "LineWidth", 1.5, "DisplayName", "Asymptotic modulation efficiency");
    xline(IWA, '--', "LineWidth", 1.5, "DisplayName", "Internal working angle");

    xline(theta_range(j), "r--", "DisplayName", "Max \theta")
    xline(theta_range(jmin), "b--", "DisplayName", "Min \theta")

    grid minor;
    legend(Location="bestoutside")

end

end


function delta = find_IWA(theta_scaled, theta_range, signal, eta_mod)
theta = theta_scaled * 1e-7;
signal_interp = interp1(theta_range, signal, theta, 'linear');
delta = 1e15 * abs(signal_interp - eta_mod);
end