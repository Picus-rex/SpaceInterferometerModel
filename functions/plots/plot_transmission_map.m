function h = plot_transmission_map(theta_range, map, export_settings)
%PLOT_TRANSMISSION_MAP Create the transmission map image given the map as
%computed by compute_response_function with settings for exporting.
%
% INPUTS:
%   theta_range[Nx1]    Range of coordinates to consider. [rad]
%   map[NxN]            Transmission map resulting from function. [-]
%  optional:
%   export_setting[struct] Setting for export. 
%
% OUTPUTS:
%   h                   Handle of the figure.
%
% VERSION HISTORY:
%   2025-03-24 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    export_settings = NaN;
end

% Call styling
style_colors;

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

if isstruct(export_settings)
    export_figures("embedded", export_settings)
end

end