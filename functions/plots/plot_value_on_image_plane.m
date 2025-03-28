function h = plot_value_on_image_plane(values, x_coords, y_coords, ...
    color_intervals, export_settings)
%PLOT_VALUE_ON_IMAGE_PLANE Create an interpolated colormap plot for given 
% values at (x, y) coordinates.
%
% INPUTS:
%   values[Nx1]          Values to be plotted using a colormap.
%   x_coords[Nx1]        X-coordinates of the values.
%   y_coords[Nx1]        Y-coordinates of the values.
%  optional:
%   color_intervals[1xM] Intervals for color mapping (default: auto scaling).
%   export_settings[struct] Setting for export.
%
% OUTPUTS:
%   h                    Handle of the figure.
%
% VERSION HISTORY:
%   2025-03-28 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 || isempty(color_intervals)
    color_intervals = []; % Use default scaling if not provided
end
if nargin < 5
    export_settings = NaN;
end

% Call styling
style_colors;

% Create scattered interpolant
F = scatteredInterpolant(x_coords, y_coords, values, 'natural', 'none');

% Define grid for interpolation
x_range = linspace(min(x_coords), max(x_coords), 100);
y_range = linspace(min(y_coords), max(y_coords), 100);
[X, Y] = meshgrid(x_range, y_range);
Z = F(X, Y);

% Create a circular mask
R = max(sqrt(x_coords.^2 + y_coords.^2)); % Radius of circular aperture
mask = sqrt(X.^2 + Y.^2) <= R;
Z(~mask) = NaN; % Mask out values outside the circular aperture

h = figure;
hold on;
pcolor(X, Y, Z);
shading interp; % Smooth interpolation
set(gca, 'YDir', 'normal'); % Correct axis orientation
xlabel('X Coordinate');
ylabel('Y Coordinate');
colormap(darkBlue);
colorbar();
axis tight;
axis square;

% Apply color intervals if provided
if ~isempty(color_intervals)
    clim([min(color_intervals), max(color_intervals)]);
end

if isstruct(export_settings)
    export_figures("embedded", export_settings);
end

end