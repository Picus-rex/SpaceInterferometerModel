function [D, a, h] = plot_apertures(positions, A, eta_opt, eta_bc, export_settings)
%PLOT_APERTURES Plot a 2D representation of apertures in a nulling 
% interferometer for visualisation.
%
% INPUTS:
%   positions[Nx2]      Matrix of (x, y) coordinates for N apertures. [m]
%   A[Nx1]              Vector of amplitudes associated to each aperture.
%                       [V/m]
%   eta_opt[1]          Optical transmission efficiency. [-]
%   eta_bc[1]           Beam combiner efficiency. [-]
%   export_setting[struct] Setting for export. 
%
% OUTPUTS:
%   D[1xN]              Possible diameters that correspond to the proposed
%                       amplitudes [m]
%   a[1xN]              Vector of the possible areas [m]
%   h[1x1]              Figure handle.
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005; 
%
% NOTES:
%   - Efficiencies have been assumed from Lay 2004.
%
% VERSION HISTORY:
%   2025-02-14 -------- 1.0
%   2025-03-17 -------- 1.1
%                     - User given efficiencies added as inputs.
%                     - export false does not create a plot.
%   2025-03-26 -------- 1.2
%                     - Uniformated to the plot structures.
%   2025-05-09 -------- 1.2.1
%                     - Fix: error in the surface computation. 
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
if nargin < 3
    eta_opt = 0.2;          % Optical transmission efficiency
    eta_bc = 0.25;          % Beam combiner efficiency
end
if nargin < 5
    export_settings = NaN;
    h = NaN;
end

% Aperture sizes relation
D = sqrt((A.^2) / (eta_opt * eta_bc * pi));
a = pi * (D / 2).^2;

if isstruct(export_settings)
    h = figure; hold on; axis equal;
    style_colors;
    
    % Plot connecting elements
    plot([positions(1,1), positions(3,1)], [positions(1,2), positions(3,2)], "-", "Color", ui_colours(1, :));
    plot([positions(2,1), positions(4,1)], [positions(2,2), positions(4,2)], "-", "Color", ui_colours(1, :));

    % For each aperture...
    for k = 1:length(A)
    
        x = positions(k,1);
        y = positions(k,2);
        r = D(k) / 2; 
        
        % Create circular aperture
        theta = linspace(0, 2*pi, 100);
        X_circle = x + r * cos(theta);
        Y_circle = y + r * sin(theta);
        
        fill(X_circle, Y_circle, colours(1, :), 'FaceAlpha', 0.5, 'EdgeColor', 'k');
        
        text(x, y, sprintf('%d', k), 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 11, 'FontName', 'Times New Roman');
    end

    plot(0, 0, '.', "Color", colours(1, :));
    text(0, 0, "Combiner", 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'FontSize', 11, 'FontName', 'Times New Roman');
    
    
    xlabel('[m]');
    ylabel('[m]');
    grid minor;
    hold off;

    export_figures("embedded", export_settings)

    for i = 1:length(D)
        fprintf("\nSize of aperture %.0f:\t\t%.3f m", i, D(i))
    end
    fprintf("\n");
end

end
