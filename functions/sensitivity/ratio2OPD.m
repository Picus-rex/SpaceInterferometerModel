function OPDs = ratio2OPD(N, amplitudes, phases, positions, ...
                                   theta_star, ratios, lambdas, autoplot)
%RATIO2OPD Compute the minimum optical path distance to have in order to
%have a desired nulling ratio for the given wavelengths interval. 
%
% INPUTS:
%   N[1]                Number of apertures. [-]
%   amplitudes[Nx1]     Amplitudes associated to each aperture.
%   phases[Nx1]         Vector of phase shifts for each aperture.
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambdas[M1x1]       Wavelengths at which verify the analysis. [m]
%   ratios[M2x1]        Desired nulling ratio(s) at which to compute the
%                       response. [-]
%   autoplot[bool]      If given and true, create plots. Default: true
%
% OUTPUTS:
%   OPDs[M1xM2]         Vector of OPDs differences for the second branch 
%                       with respect to the first branch (minimum value,
%                       there is a period as stated in find_minimum_OPD...
%
% NOTES:
%   - OPDs are periodic, therefore each multiple of n * OPD_max, with n
%     integer number will produce the same nulling ratio. 
%
% VERSION HISTORY:
%   2025-03-13 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input management
if nargin < 6
    ratios = logspace(-11, 0, 100);
end
if nargin < 7
    lambdas = linspace(1e-6, 20e-6, 100);
end
if nargin < 8
    autoplot = true;
end

% Allocate output and compute the response by changing the wavelength
OPDs = zeros(length(ratios), length(lambdas));

for i = 1:length(lambdas)
    [OPDs(:, i), ~] = find_minimum_OPD_single_branch(N, amplitudes, ...
               phases, positions, theta_star, lambdas(i), ratios, false);
end

if autoplot
    figure;
    imagesc(ratios, lambdas, log10(OPDs' + eps)); % Add eps to avoid log(0)
    
    style_colors;

    % Set the color limits for the logarithmic scale
    clim([log10(min(OPDs(:) + eps)), log10(max(OPDs(:)))]);
    colormap(darkBlueZero)
    cb = colorbar;
    
    % Update thicks 
    tick_values = get(cb, 'Ticks'); 
    tick_labels = arrayfun(@(x) sprintf('10^{%.0f}', x), tick_values, 'UniformOutput', false);
    set(cb, 'TickLabels', tick_labels);

    % Set axis labels
    xlabel('Nulling Ratios');
    ylabel('Wavelengths [m]');
    title("OPD")
    set(gca, 'XScale', 'log');
end

end