function h = plot_value_on_image_plane(values, x_coords, y_coords, ...
    color_intervals, export_settings)
%PLOT_VALUE_ON_IMAGE_PLANE Create an interpolated colormap plot for given 
% values at (x, y) coordinates.
%
% INPUTS:
%   values[Nx1]         Values to be plotted using a colormap.
%   x_coords[Nx1]       X-coordinates of the values.
%   y_coords[Nx1]       Y-coordinates of the values.
%  optional:
%   color_intervals[1xM]Intervals for color mapping (default: auto scaling).
%   embedded[struct]    Setting for export.
%   type[string]        Can be either "angles" to convert the values to
%                       degrees or a string following "linear_%de(+-)%d",
%                       where the exponential number is reversed and
%                       recognised as a scale (if 1e-6 is given, for
%                       example, the string is converted to micro m). 
%   distance_units[string] Conversion to apply to the X, Y values. Can be
%                       either [m] or [cm].
%   title[string]       Title of the plot. The unit of the plotted value is
%                       automatically added. 
%
% OUTPUTS:
%   h                   Handle of the figure.
%
% VERSION HISTORY:
%   2025-03-28 -------- 1.0
%   2025-03-31 -------- 1.1
%                     - Added several plot settings, including type and
%                       scale of values and titles. 
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    values (:, :)
    x_coords (:, :)
    y_coords (:, :)
    color_intervals (1, :) = [] 
    export_settings.embedded = NaN;
    export_settings.type = "linear_1e-6";
    export_settings.distance_units = "cm";
    export_settings.title = "";
end

% Call styling
style_colors;

% Create scattered interpolant
F = scatteredInterpolant(x_coords, y_coords, values, 'natural', 'none');

% Define grid for interpolation
x_range = linspace(min(x_coords), max(x_coords), 1000);
y_range = linspace(min(y_coords), max(y_coords), 1000);
[X, Y] = meshgrid(x_range, y_range);
Z = F(X, Y);

% If the input is an angle, change the visualisation
if strcmp(export_settings.type, "angles")
    Z = rad2deg(Z);
    selected_colormap = angle_darkBlue;
elseif startsWith(export_settings.type, "linear")
    matches = regexp(export_settings.type, '_(\d+e[+-]?\d+)', 'tokens');
    if ~isempty(matches)
        scale = str2double(matches{1}{1});
    else
        scale = 1;
    end
    Z = Z * 1/scale;
    selected_colormap = darkBlue;
end

% Create a circular mask
R = max(sqrt(x_coords.^2 + y_coords.^2)); % Radius of circular aperture
mask = sqrt(X.^2 + Y.^2) <= R;
Z(~mask) = NaN; % Mask out values outside the circular aperture

if strcmp(export_settings.distance_units, "cm")
    X = X * 1e2;
    Y = Y * 1e2;
    unit_string = "[cm]";
elseif strcmp(export_settings.distance_units, "m")
    unit_string = "[m]";
else
    warning("This distance_units value is unknown. Defaulting to metres.")
    unit_string = "[m]";
end

h = figure;
hold on;
pcolor(X, Y, Z);
shading interp; % Smooth interpolation
set(gca, 'YDir', 'normal'); % Correct axis orientation
xlabel('X Coordinate ' + unit_string);
ylabel('Y Coordinate '+ unit_string);
colormap(selected_colormap);
colorbar();
axis tight;
axis square;

if ~strcmp(export_settings.title, "")
    if startsWith(export_settings.type, "linear")
        switch scale
            case 1e-6
                scale_name = " [µm]";
            otherwise
                scale_name = sprintf(" [%.0d m]", scale);
        end
    elseif strcmp(export_settings.type, "angles")
        scale_name = " [deg]";
    else
        scale_name = " [m]";
    end
    title(export_settings.title + scale_name);
end

% Apply color intervals if provided
if ~isempty(color_intervals)
    clim([min(color_intervals), max(color_intervals)]);
end

if isstruct(export_settings) && isstruct(export_settings.embedded)
    export_figures("embedded", export_settings);
end

end