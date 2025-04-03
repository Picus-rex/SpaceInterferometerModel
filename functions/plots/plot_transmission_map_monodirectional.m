function h = plot_transmission_map_monodirectional(theta_range, maps, names, embedded, export_settings)
%PLOT_TRANSMISSION_MAP_MONODIRECTIONAL Create the transmission map along
%the horizontal axis given the map as computed by compute_response_function 
%with settings for exporting.
%
% INPUTS:
%   theta_range[Nx1]    Range of coordinates to consider. [rad]
%   maps[NxN]           Selection of desired transmission maps resulting 
%                       from function. [-]
%  optional:
%   export_setting[struct] Setting for export. 
%
% OUTPUTS:
%   h                   Handle of the figure.
%
% VERSION HISTORY:
%   2025-03-24 -------- 1.0
%   2025-04-03 -------- 1.1
%                     - Added arguments and possibility to change the
%                       number of tendencies lines within the plot. 
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    theta_range (:, :)
    maps {mustBeA(maps, "table")}
    names {mustBeA(names, "string")}
    embedded = NaN
    export_settings.tendencies = [2, 4, 6];
end

% Call styling
style_colors;
cols = get_colours(length(names));

% Conversion factor and definition
conversion_rad2mas = 1e3 * (3600 * 180) / pi;
theta_range = theta_range * conversion_rad2mas;
maps_names = fieldnames(maps);
T = maps.(maps_names{1});

h = figure; 
hold on;
for i = 1:length(names)
    semilogy(theta_range, maps.(maps_names{i})(floor(size(T, 1)/2), :), ...
        "LineWidth", 1.5, "DisplayName", names{i}, "Color", cols(i, :)); 
end

% Scaling factor C for theta^2 and theta^4
for i = 1:length(export_settings.tendencies)
    C = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / ...
        theta_range(floor(size(T, 1)/2))^export_settings.tendencies(i);
    name = sprintf("\\theta^{%d}", export_settings.tendencies(i));
    semilogy(theta_range, C * theta_range.^(export_settings.tendencies(i)), ...
     "--", 'DisplayName', name, "LineWidth", 1.5, "Color", ui_colours(i, :));
end

xlabel('\theta_x [mas]');
ylabel('Normalized Intensity');
legend;
grid minor; 
axis tight;
set(gca, 'YScale', 'log');
hold off;

if isstruct(embedded)
    export_figures("embedded", embedded)
end

end