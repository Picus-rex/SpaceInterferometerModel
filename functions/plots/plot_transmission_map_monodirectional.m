function h = plot_transmission_map_monodirectional(theta_range, maps, names, export_settings)
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
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 3
    export_settings = NaN;
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
for i = 1:length(maps_names)-3
    semilogy(theta_range, maps.(maps_names{i})(floor(size(T, 1)/2), :), "LineWidth", 1.5, "DisplayName", names{i}, "Color", cols(i, :)); 
end

% Scaling factor C for theta^2 and theta^4
C2 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^2;
C4 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^4;
C6 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^6;

semilogy(theta_range, C2 * theta_range.^(2), "--", 'DisplayName', "\theta^{2}", "LineWidth", 1.5, "Color", ui_colours(1, :));
semilogy(theta_range, C4 * theta_range.^(4), "--", 'DisplayName', "\theta^{4}", "LineWidth", 1.5, "Color", ui_colours(2, :));
semilogy(theta_range, C6 * theta_range.^(6), "--", 'DisplayName', "\theta^{6}", "LineWidth", 1.5, "Color", ui_colours(3, :));

xlabel('\theta_x [mas]');
ylabel('Normalized Intensity');
legend;
grid minor; 
axis tight;
set(gca, 'YScale', 'log');
hold off;

if isstruct(export_settings)
    export_figures("embedded", export_settings)
end

end