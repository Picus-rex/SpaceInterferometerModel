function [modulation, hs] = planet_modulation(data, export_settings)
%PLANET_MODULATION Computes the modulation of the planet signal by
% rotating the array and measuring the interpolated response function at
% the planet's position.
%
% INPUTS:
%   data[struct]        For integrated development, all the inputs can be
%                       grouped into a single struct
%  optional:
%   export_setting[struct] Setting for export. 
%
% OUTPUTS:
%   modulation[1x n_rotation] Modulated signal values at each rotation step
%   hs                  Handles of figures.
%
% NOTES:
%  - If field planet_modulation_positions is not specified in simulation,
%    then the function acts like compute_response_function. 
%  - Only use a square meshgrid for the analysis.
%
% VERSION HISTORY:
%   2025-03-07 -------- 1.0
%   2025-03-10 -------- 1.1
%                     - If planet_modulation_positions = 0, skip function
%                     - Typo in error function.
%   2025-03-20 -------- 1.2
%                     - Modulation removed because incorrect. Use
%                       compute_modulation_efficiency to extract modulation 
%                       efficiency.
%   2025-03-24 -------- 1.3
%                     - Added export settings to create coherent plots.
%                     - Print angles.
%   2025-03-29 -------- 1.3.1
%                     - Fixed extraction of map following changes in
%                       compute_response_function and fixed problem with
%                       exporting of figures.
%                     - Planet signal modulation capped at 360° along the x
%                       axis on the plot. 
%   2025-05-09 -------- 1.3.2
%                     - Fixed missing export of hs when 
%                       planet_modulation_positions was not given. 
%   2025-05-12 -------- 1.3.3
%                     - Allocate outputs before functions to avoid errors.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    export_settings = NaN;
    exp = false;
else
    exp = true;
end

% Conversion factor and styling for the plots
conversion_rad2mas = 1e3 * (3600 * 180) / pi;
style_colors;

% Allocate outputs
hs = NaN;

% Planet position
planet_x = data.environment.exoplanet_position(1);
planet_y = data.environment.exoplanet_position(2);

% Initialize modulation vector; leave function if no rotations
if isfield(data.simulation, "planet_modulation_positions")
    n_rotation = data.simulation.planet_modulation_positions;
    if n_rotation == 0
        modulation = NaN;
        return
    end
else
    n_rotation = 1;
end
modulation = zeros(1, n_rotation);
if isfield(data.simulation, "monte_carlo_iterations") && data.simulation.monte_carlo_iterations > 0
    modulation_err = zeros(1, n_rotation);
end

% Theta grid
theta_x = data.simulation.theta_x;
theta_y = data.simulation.theta_y;

if planet_x > max(theta_x, [], "all") || planet_x < min(theta_x, [], "all") || ...
        planet_y > max(theta_y, [], "all") || planet_y < min(theta_y, [], "all")
    error("The exoplanet is outside the angular extension of the simulation")
end

% Number of plots to do
div_plot = floor(n_rotation/5);

% Rotation angles
angles = linspace(0, 2*pi, n_rotation);
original_positions = data.instrument.positions;

% Start debug plots if needed
if isfield(data.outputs, "plot_intermediate_rotations") && ...
            data.outputs.plot_intermediate_rotations
    L = ceil(floor(n_rotation/div_plot).^0.5);
    if exp
        hs = [];
    else
        h1 = figure;
        h2 = figure;
    end
    j = 1;
end

% For every rotation...
for i = 1:n_rotation

    if isfield(data.outputs, "verbose") && data.outputs.verbose && ~mod(i, div_plot)
        fprintf("Rotation %.0f of %.0f at angle %.2f\n", i, n_rotation, rad2deg(angles(i)))
    end

    % Rotate the apertures
    rotation_matrix = [cos(angles(i)), -sin(angles(i)); 
                       sin(angles(i)), cos(angles(i))];

    data.instrument.positions = (rotation_matrix * original_positions')';

    % Compute the new response function
    [Maps] = compute_response_function("data", data);
    T_rotated = Maps.T_standard;
    if isfield(data.simulation, "monte_carlo_iterations") && data.simulation.monte_carlo_iterations > 0
        T_real_rotated = Maps.T_real;
    end

    % Plot intermediate maps
    if isfield(data.outputs, "plot_intermediate_rotations") && ...
            data.outputs.plot_intermediate_rotations
        if ~mod(i, div_plot) || i == 1
            
            % Figure 1: Transmission maps
            if exp
                hs(j) = figure(); 
                j = j + 1;
            else
                figure(h1)
                subplot(L, L, j);
            end
            hold on; axis equal;
            imagesc(theta_x(1, :) * conversion_rad2mas, theta_y(:, 1) * conversion_rad2mas, T_rotated);
            colormap(darkBlue)
            colorbar();
            plot(planet_x * conversion_rad2mas, planet_y * conversion_rad2mas, "*", "MarkerSize", 10, "Color", ui_colours(5, :));
            xlabel('\theta_x [mas]');
            ylabel('\theta_y [mas]');
            axis square;
            axis tight;
            if exp
                export_figures("width", export_settings.sizes.width.(export_settings.include.array_plots.width), ...
                    "height", export_settings.sizes.height.(export_settings.include.array_plots.height), ...
                    "name", export_settings.name + "_map_" + string(j))
            else
                title(sprintf("Rotation: %.2f deg", rad2deg(angles(i))))
            end
            hold off;
            
            % Figure 2: rotated arrays   
            if exp
                hs(j) = figure();
            else
                figure(h2)
                subplot(L, L, j); 
            end
            hold on; axis equal;
            for k = 1:size(data.instrument.positions, 1)
                plot(data.instrument.positions(k, 1), data.instrument.positions(k, 2), ".", "MarkerSize", 10, "Color", colours(1, :))
            end
            xlabel('[m]');
            ylabel('[m]');
            grid minor;
            axis tight;
            if strcmp(data.instrument.array, "X-Array")
                plot([data.instrument.positions(1, 1), data.instrument.positions(3, 1)], [data.instrument.positions(1, 2), data.instrument.positions(3, 2)], "-", "LineWidth", 1.5, "Color", colours(3, :))
                plot([data.instrument.positions(2, 1), data.instrument.positions(4, 1)], [data.instrument.positions(2, 2), data.instrument.positions(4, 2)], "-", "LineWidth", 1.5, "Color", colours(3, :))
            else
                plot([data.instrument.positions(1, 1), data.instrument.positions(4, 1)], [data.instrument.positions(1, 2), data.instrument.positions(4, 2)], "-", "LineWidth", 1.5, "Color", colours(3, :))
            end
            if exp
                export_figures("width", export_settings.sizes.width.(export_settings.include.array_plots.width), ...
                    "height", export_settings.sizes.height.(export_settings.include.array_plots.height), ...
                    "name", export_settings.name + "_array_" + string(j))
            else
                title(sprintf("Rotation: %.2f deg", rad2deg(angles(i))))
            end

            hold off;

            j = j + 1;
        end
    end
    
    % Interpolate the signal at the planet position
    modulation(i) = interp2(theta_x, theta_y, T_rotated, planet_x, planet_y, 'linear', 0);

    if isfield(data.simulation, "monte_carlo_iterations") && data.simulation.monte_carlo_iterations > 0
        modulation_err(i) = interp2(theta_x, theta_y, T_real_rotated, planet_x, planet_y, 'linear', 0);
    end
end

% Plot the modulation curve
if isfield(data.outputs, "plot_modulation") && data.outputs.plot_modulation
    if exp
        hs(end+1) = figure();
    else
        h3 = figure();
    end

    hold on;
    plot(rad2deg(angles), modulation, '-', 'LineWidth', 2, "Color", colours(1, :));

    if isfield(data.simulation, "monte_carlo_iterations") && data.simulation.monte_carlo_iterations > 0
        plot(rad2deg(angles), modulation_err, '-', 'LineWidth', 2, "Color", Colours(2, :));
    end
    
    xlim([0, 360]);
    xlabel('Rotation Angle [deg]');
    ylabel('Interpolated Signal at Planet Position');
    grid minor;

    if exp
        export_figures("embedded", export_settings)
    else
        title('Planet Signal Modulation');
        hs = [h1, h2, h3];
    end
end

end
