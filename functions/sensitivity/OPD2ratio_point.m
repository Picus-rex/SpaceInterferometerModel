function ratios = OPD2ratio_point(N, amplitudes, phases, positions, ...
    theta_star, lambdas, point, autoplot)
%OPD2RATIO_POINT By varying the wavelength with defined shifts on the two
%branches of a double Bracewell telescopes, find the resulting nulling
%ratios for the given interval of wavelengths.
%
% INPUTS:
%   N[1]                Number of apertures. [-]
%   amplitudes[Nx1]     Amplitudes associated to each aperture.
%   phases[Nx1]         Vector of phase shifts for each aperture.
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambdas[Mx1]        Wavelengths of observation. By default, 1000 points
%                       between 0.1 and 10 micro m are chosen. [m]
%   point[2x1]          Shift to impose on the 1-3 and 2-4 branch
%                       respectively. By default, 2, 1 micro m. [m]
%   autoplot[bool]      If given and true, create plots. Default: true
%
% OUTPUTS:
%   ratios[M]           Nulling ratios for each wavelength in the specified
%                       point.
%
% NOTES:
%   - OPDs are periodic, therefore each multiple of n * lambda, with n
%     integer number will produce the same nulling ratio. Therefore, for an
%     optimal result, the point should be below the wavelength.
%
% VERSION HISTORY:
%   2025-03-17 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    lambdas = linspace(1e-6, 1e-4, 1000);
    point = [2; 1] * 1e-7;
    autoplot = true;
elseif nargin < 7
    point = [2; 1] * 1e-7;
    autoplot = true;
elseif nargin < 8
    autoplot = true;
end

points = length(lambdas);
ratios = zeros(points, 1);
opds = [point(1), point(2), 0, 0];

if max(opds) > min(lambdas)
    warning("The plot will be inconsistent because the maximum specified OPD point is larger than the period of shifting, returning results at a later repetition.")
end

for k = 1:points
    lambda = lambdas(k);
    
    % The phases are bound to 2pi.
    shift_phases = wrapTo2Pi(2*pi*opds/lambda + phases);
    ratios(k) = compute_nulling_ratio(N, amplitudes, ...
                shift_phases, positions, theta_star, lambda);
end


if autoplot
    
    style_colors;
    
    figure; hold on;
    plot(lambdas * 1e6, ratios, "LineWidth", 1.5);
    xlabel("Wavelength [micro m]");
    ylabel("Nulling ratio");

    axis tight;   
    set(gca, 'XScale', 'log', 'YScale', 'log');

end

end
