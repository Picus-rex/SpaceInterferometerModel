function ratios = OPDs2ratio(N, amplitudes, phases, positions, theta_star, ...
                                                        lambdas, autoplot)
%OPDS2RATIO By varying the OPDs on the two branches of a double Bracewell
%telescopes, find the resulting nulling ratios for the given interval of
%wavelengths.
%
% INPUTS:
%   N[1]                Number of apertures. [-]
%   amplitudes[Nx1]     Amplitudes associated to each aperture.
%   phases[Nx1]         Vector of phase shifts for each aperture.
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambdas[Mx1]        Wavelengths of observation. By default, a single 
%                       observation at 10e-6 is performed. [m]
%   autoplot[bool]      If given and true, create plots. Default: true
%
% OUTPUTS:
%   ratios[Mx1000x1000] Nulling ratios for each wavelength and for each
%                       combination of shifts on the axes.
%
% NOTES:
%   - OPDs are periodic, therefore each multiple of n * OPD_max, with n
%     integer number will produce the same nulling ratio. 
%
% VERSION HISTORY:
%   2025-03-14 -------- 1.0
%   2025-03-17 -------- 1.1
%                     - Fixed plots consistency across the scripts.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    lambdas = 10e-6;
    autoplot = true;
elseif nargin < 7
    autoplot = true;
end

regen = true;

if isfile('exports/ratios.mat')

    data = load("exports/ratios.mat", "lambdas");
    if ~all(lambdas == data.lambdas)
        fprintf("Saved data existed, but do not match input wavelengths. Regenerating...\n")
    else
        regen = false;
        load("exports/ratios.mat", "lambdas", "ratios", "OPD_max", "opd");
    end

end

if regen
    fprintf("Saved file does not exist. The generation of the map requires few seconds...\n")

    points = 1000;
    ratios = zeros(length(lambdas), points, points);

    opd = zeros(length(lambdas), points);
    
    for k = 1:length(lambdas)

        lambda = lambdas(k);

        % The phases are bound to 2pi, therefore compute the maximum OPD to
        % have the period extension. After OPD_max, everything repeats. 
        % Compute the opd range of points.
        OPD_max = lambda;
        opd(k, :) = logspace(-16, log10(OPD_max), points);
        
        for i = 1:points
            for j = 1:points
                opds = [opd(k, i), opd(k, j), 0, 0];
                shift_phases = wrapTo2Pi(2*pi*opds/lambda + phases);
                ratios(k, i, j) = compute_nulling_ratio(N, amplitudes, ...
                            shift_phases, positions, theta_star, lambda);
            end
        end

    end
        
    fprintf("Saving file...\n")
    save("exports/ratios.mat", "lambdas", "ratios", "OPD_max", "opd")

end

if autoplot
    
    style_colors;

    global_min = min(log10(squeeze(ratios(:, :, :)) + eps), [], 'all');
    global_max = max(log10(squeeze(ratios(:, :, :)) + eps), [], 'all');
    
    for i = 1:length(lambdas)
        figure; hold on;
        imagesc(opd(i, :), opd(i, :), log10(squeeze(ratios(i, :, :)) + eps));
        title(sprintf("OPD [%.000f \\micro m]", lambdas(i)*1e6), "Interpreter", "latex");
        xlabel("OPD difference branch 1-3 [m]");
        ylabel("OPD difference branch 2-4 [m]");
        colormap(darkBlue)
        cb = colorbar();
        clim([global_min global_max]);

        % Update thicks 
        tick_values = get(cb, 'Ticks'); 
        tick_labels = arrayfun(@(x) sprintf('10^{%.0f}', x), tick_values, 'UniformOutput', false);
        set(cb, 'TickLabels', tick_labels);

        axis tight;        
    end

end

end
