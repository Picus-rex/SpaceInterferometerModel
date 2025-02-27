function [U, amplitudes, phases] = ...
        compute_optimal_splitting(N, B, lambda, x_ap, y_ap, theta_max, autoplot, export)
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
%   export[bool]        If true, save figures
%
% [^1] The baseline is only used to compute the angular extension of the
% star in accordance with the paper. 
%
% OUTPUTS:
%   U[NxN]              Matrix of the beam combination, with unitary
%                       determinant. It includes the phases delay if
%                       necessary. [-]
%   amplitudes[NxN]     Beam combination amplitudes for optimal nulling
%                       (U's last row). [-]
%   phases[NxN]         Beam combination phases for optimal nulling. [rad]
%
% REFERENCES:
%   Guyon O, Mennesson B, Serabyn E, Martin S. Optimal beam combiner design
%   for nulling interferometers. Publications of the Astronomical Society
%   of the Pacific. 2013 Aug;125(930):951–65.
%
% VERSION HISTORY:
%   2025-02-11 -------- 1.0
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
    export = false;
end           

% Definition of vectors and empty placeholders
theta_x = linspace(-theta_max, theta_max, 100);
theta_y = linspace(-theta_max, theta_max, 100);
[THETA_X, THETA_Y] = meshgrid(theta_x, theta_y);
Npt = numel(THETA_X);
A = zeros(Npt, N); 

% For each sampled point...
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

if autoplot
    
    % Compute the transformed field at the output and normalised ity
    W = U' * A';  
    W_intensity = abs(W).^2;
    W_intensity = W_intensity / max(W_intensity(:));
    
    % Plot transmission maps
    figure;
    for m = 1:min(4, N)
        subplot(1, 4, m);
        imagesc(theta_x, theta_y, ...
            reshape(abs(W_intensity(m, :)), sqrt(Npt), sqrt(Npt)));
        colorbar; styling; % colormap winter;
        xlabel('\theta_x [rad]');
        ylabel('\theta_y [rad]');
        title(['Mode ', num2str(m)]);
    end
    
    if export
        for m = 1:min(4, N)
            h = figure;
            imagesc(theta_x, theta_y, ...
                reshape(abs(W_intensity(m, :)), sqrt(Npt), sqrt(Npt)));
            colorbar; styling(true, 6, 6, sprintf("exports/optimal_modes_%.0f.pdf", m+4)); % colormap winter;
            xlabel('\theta_x [rad]');
            ylabel('\theta_y [rad]');

            % set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 6, 6]);
            % set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
            % set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
            % set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [6, 6]);
            % set(gcf, 'PaperPositionMode', 'auto');
            % print(gcf, sprintf("exports/optimal_modes_%.0f.pdf", m+4), '-dpdf', '-r300');
            close(h)
        end
    end

    % Plot intensity waves
    figure;
    i_y0 = round(length(theta_y) / 2);
    for m = 1:min(4, N)
        subplot(3, 2, m);
        % Reshape W_intensity into 2D and extract the middle row
        W_reshaped = reshape(W_intensity(m, :), length(theta_y), length(theta_x));
        intensity_profile = W_reshaped(i_y0, :);
        plot(theta_x, intensity_profile, "LineWidth", 1.5, 'DisplayName', ['Mode ', num2str(m)]);
        xlabel('\theta_x (\lambda/B)');
        ylabel('Intensity profile');
        title(['Mode ', num2str(m)]);
    end

    subplot(3, 2, [5, 6]);
    for m = 1:min(4, N)
        W_reshaped = reshape(W_intensity(m, :), length(theta_y), length(theta_x));
        intensity_profile = W_reshaped(i_y0, :);
        semilogy(theta_x, intensity_profile, 'DisplayName', ['Mode ', num2str(m)], "LineWidth", 1.5);
        hold on;
    end

    xlabel('\theta_x (\lambda/B)');
    ylabel('Normalized Intensity');
    legend;
    grid minor; hold off;
end

end