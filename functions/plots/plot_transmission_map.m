function h = plot_transmission_map(theta_range, map, export_settings, planets)
%PLOT_TRANSMISSION_MAP Create the transmission map image given the map as
%computed by compute_response_function with settings for exporting.
%
% INPUTS:
%   theta_range[Nx1]    Range of coordinates to consider. [rad]
%   map[NxN]            Transmission map resulting from function. [-]
%  optional:
%   export_setting[struct] Setting for export. 
%   planets[Npx1]       Angular positions of the planets to plots [mas]
%
% OUTPUTS:
%   h                   Handle of the figure.
%
% VERSION HISTORY:
%   2025-03-24 -------- 1.0
%   2025-05-05 -------- 1.1
%                     - Integrated planet plotting if desired.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    export_settings = NaN;
    planets = [];
end

% Call styling
style_colors;

% For planet plotting
Ns = length(planets);
angles = linspace(0, 2*pi, 10000); 

% Conversion factor
conversion_rad2mas = 1e3 * (3600 * 180) / pi;

% Convert data and generate empty arrays
theta_range = theta_range * conversion_rad2mas;

h = figure; 
hold on;
imagesc(theta_range, theta_range, map);
xlabel('\theta_x [mas]');
ylabel('\theta_y [mas]');
colormap(darkBlue)
colorbar();
axis xy; 
axis square;
axis tight;

for p = 1:Ns
    planet_x = planets(p) * cos(angles);
    planet_y = planets(p) * sin(angles);
    plot(planet_x, planet_y, 'w--', 'LineWidth', 1.5);
end


if isstruct(export_settings)
    export_figures("embedded", export_settings)
end

end