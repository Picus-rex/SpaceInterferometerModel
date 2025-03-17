function [D, a] = plot_apertures(positions, A, eta_opt, eta_bc, export)
%PLOT_APERTURES Plot a 2D representation of apertures in a nulling 
% interferometer for visualisation.
%
% INPUTS:
%   positions[Nx2]      Matrix of (x, y) coordinates for N apertures. [m]
%   A[Nx1]              Vector of amplitudes associated to each aperture.
%                       [V/m]
%   eta_opt[1]          Optical transmission efficiency. [-]
%   eta_bc[1]           Beam combiner efficiency. [-]
%   export[bool]        If true, export plot and results.
%
% OUTPUTS:
%   D[1xN]              Possible diameters that correspond to the proposed
%                       amplitudes [m]
%   a[1xN]              Vector of the possible areas [m]
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
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
if nargin < 3
    eta_opt = 0.2;          % Optical transmission efficiency
    eta_bc = 0.25;          % Beam combiner efficiency
end

% Aperture sizes relation
a = (A.^2) / (eta_opt * eta_bc);
D = 2 * sqrt(a / pi);

if export
    figure; hold on; axis equal;
    colours = styling; colorbar("off");
    
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

    for i = 1:length(D)
        fprintf("\nSize of aperture %.0f:\t\t%.3f m", i, D(i))
    end
    fprintf("\n");
end

end
