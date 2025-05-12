function [U, amplitudes, phases, hs] = compute_optimal_splitting(N, B, ...
    lambda, x_ap, y_ap, theta_max, autoplot, export_settings)
%COMPUTE_OPTIMAL_SPLITTING Following the methodology illustrated in Guyon
%O., this function computes the splitting matrix U, optimally computed
%assuming a certain stellar angular extension.
%
% INPUTS:
%   N[1]                Number of apertures [-]
%   B[1]                Larger baseline of telescope^[1] [m]
%   lambda[1]           Reference wavelength [m]
%   x_ap[N]             Vector with abscissa positions of apertures [m]
%   y_ap[N]             Vector with ordinate positions of apertures [m]
%   theta_max[1]        Angular star extension [rad]
%   autoplot[bool]      If true, produce a plot of the different modes
%   export_setting[struct] Setting for export. 
%
% [^1] The baseline is only used to compute the angular extension of the
% star in accordance with the paper. 
%
% OUTPUTS:
%   U[NxN]              Matrix of the beam combination, with unitary
%                       determinant. It includes the phases delay if
%                       necessary. [-]
%   amplitudes[1xN]     Beam combination amplitudes for optimal nulling
%                       (U's last row). [-]
%   phases[1xN]         Beam combination phases for optimal nulling. [rad]
%
% REFERENCES:
%   Guyon O, Mennesson B, Serabyn E, Martin S. Optimal beam combiner design
%   for nulling interferometers. Publications of the Astronomical Society
%   of the Pacific. 2013 Aug;125(930):951–65.
%
% VERSION HISTORY:
%   2025-02-11 -------- 1.0
%   2025-03-24 -------- 1.1
%                     - Updated to reflect the export settings as in other
%                       functions that create plots.
%   2025-03-25 -------- 1.1.1
%                     - Optimisation of the plotting code.
%                     - Added colours support.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6 || theta_max == 0
    % Compute angular extension of the star as indicated in the referenced
    % paper. This is a simplification that ensure the independency of the
    % solution from the star dimension.
    theta_max = 0.001 * lambda / B;
end
if nargin < 7
    autoplot = false;
end
if nargin < 8
    export_settings = NaN;
    export = false;
else
    export = true;
end

% Definition of vectors and empty placeholders
theta_x = linspace(-theta_max, theta_max, 100);
theta_y = linspace(-theta_max, theta_max, 100);
[THETA_X, THETA_Y] = meshgrid(theta_x, theta_y);
Npt = numel(THETA_X);
A = zeros(Npt, N); 

% For each sampled point compute the response following the extended
% definition without added phases.
for pt = 1:Npt
    alpha = THETA_X(pt);
    beta = THETA_Y(pt);
    A(pt, :) = exp(1i * 2 * pi * (x_ap * alpha + y_ap * beta) / lambda);
end

% Compute SVD in the cheap way.
[U, ~, ~] = svd(A', 'econ');  

% Export amplitude and phase
amplitudes = abs(U(:, end)).';
phases = angle(U(:, end)).';

% Plot options
if autoplot
    
    % Get styling, initialise options
    style_colors;
    hs = [];
    conversion_rad2mas = 1e3 * (3600 * 180) / pi;
    theta_x = theta_x * conversion_rad2mas;
    theta_y = theta_y * conversion_rad2mas;

    % Compute the transformed field at the output and normalised it
    W = U' * A';  
    W_intensity = abs(W).^2;
    W_intensity = W_intensity / max(W_intensity(:));
    
    % PLOT 1: transmission maps
    if ~export
        hs(end+1) = figure;
    end
    for m = 1:min(4, N)
        if export
            hs(end+1) = figure;
        else
            subplot(1, 4, m);
        end
        imagesc(theta_x, theta_y, ...
                reshape(abs(W_intensity(m, :)), sqrt(Npt), sqrt(Npt)));
        colormap(darkBlue)
        colorbar();
        xlabel('\theta_x [mas]');
        ylabel('\theta_y [mas]');
        axis square;
        axis tight;
        if export
            export_figures("font_size", export_settings.sizes.font_size, ...
                "height", export_settings.sizes.height.(export_settings.include.modes_maps.height), ...
                "width", export_settings.sizes.width.(export_settings.include.modes_maps.width), ...
                "name", export_settings.name + "_" + string(length(hs)));
        else
            title(['Mode ', num2str(m)]);
        end
    end

    % PLOT 2: response on horizontal axis
    i_y0 = round(length(theta_y) / 2);
    if ~export
        hs(end+1) = figure;
    end
    for m = 1:min(4, N)
        if export
            hs(end+1) = figure; hold on;
        else

        end
        % Reshape W_intensity into 2D and extract the middle row
        W_reshaped = reshape(W_intensity(m, :), length(theta_y), length(theta_x));
        intensity_profile = W_reshaped(i_y0, :);
        plot(theta_x, intensity_profile, "LineWidth", 1.5, 'DisplayName', ['Mode ', num2str(m)], 'Color', colours(1, :));
        xlabel('\theta_x [mas]');
        ylabel('Intensity profile');
        grid minor;
        axis tight;
        if export
            export_figures("font_size", export_settings.sizes.font_size, ...
                "height", export_settings.sizes.height.(export_settings.include.modes_comparison.height), ...
                "width", export_settings.sizes.width.(export_settings.include.modes_comparison.width), ...
                "name", export_settings.name + "_mode_" + string(length(hs)));
        else
            title(['Mode ', num2str(m)]);
        end
    end
    if export
        hs(end+1) = figure;
    else
        subplot(3, 2, m);
    end
    for m = 1:min(4, N)
        W_reshaped = reshape(W_intensity(m, :), length(theta_y), length(theta_x));
        intensity_profile = W_reshaped(i_y0, :);
        semilogy(theta_x, intensity_profile, 'DisplayName', ['Mode ', num2str(m)], "LineWidth", 1.5, 'Color', ui_colours(m, :));
        hold on;
    end
    xlabel('\theta_x [mas]');
    ylabel('Normalized Intensity');
    ylim([0, 1.2*max(reshape(W_intensity(1, :), length(theta_y), length(theta_x)))])
    legend("NumColumns", 2);
    grid minor; hold off;
    if export
        export_figures("embedded", export_settings, "name", export_settings.name)
    else
        title("Mode comparison")
    end
end

end