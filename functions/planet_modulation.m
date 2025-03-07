function [modulation, efficiency] = planet_modulation(data)
%PLANET_MODULATION Computes the modulation of the planet signal by
% rotating the array and measuring the interpolated response function at
% the planet's position.
%
% INPUTS:
%   data[struct]        For integrated development, all the inputs can be
%                       grouped into a single struct
%
% OUTPUTS:
%   modulation[1x n_rotation] Modulated signal values at each rotation step
%   efficiency[1]       Modulation efficiency.
%
% NOTES:
%  - If field planet_modulation_positions is not specified in simulation,
%    then the function acts like compute_response_function. 
%  - Only use a square meshgrid for the analysis.
%
% VERSION HISTORY:
%   2025-03-07 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Planet position
planet_x = data.environment.exoplanet_position{1}(1);
planet_y = data.environment.exoplanet_position{2}(1);

% Theta grid
theta_x = data.simulation.theta_x;
theta_y = data.simulation.theta_y;

if planet_x > max(theta_x, [], "all") || planet_x < min(theta_x, [], "all") || ...
        planet_y > max(theta_y, [], "all") || planet_y < min(theta_y, [], "all")
    eror("The exoplanet is outside the angular extension of the simulation")
end

% Initialize modulation vector
if isfield(data.simulation, "planet_modulation_positions")
    n_rotation = data.simulation.planet_modulation_positions;
else
    n_rotation = 1;
end
modulation = zeros(1, n_rotation);
if isfield(data.simulation, "monte_carlo_iterations") && data.simulation.monte_carlo_iterations > 0
    modulation_err = zeros(1, n_rotation);
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
    h1 = figure;
    h2 = figure;
    j = 1;
end

% For every rotation...
for i = 1:n_rotation

    if data.outputs.verbose && ~mod(i, div_plot)
        fprintf("Rotation %.0f of %.0f\n", i, n_rotation)
    end

    % Rotate the apertures
    rotation_matrix = [cos(angles(i)), -sin(angles(i)); 
                       sin(angles(i)), cos(angles(i))];

    data.instrument.positions = (rotation_matrix * original_positions')';

    % Compute the new response function
    if isfield(data.simulation, "monte_carlo_iterations") && data.simulation.monte_carlo_iterations > 0
        [T_rotated, ~, T_real_rotated] = compute_response_function("data", data);
    else
        T_rotated = compute_response_function("data", data);
    end

    % Plot intermediate maps
    if isfield(data.outputs, "plot_intermediate_rotations") && ...
            data.outputs.plot_intermediate_rotations
        if ~mod(i, div_plot) || i == 1
            
            figure(h1)
            subplot(L, L, j); hold on;
            imagesc(theta_x(1, :), theta_y(:, 1), T_rotated);
            plot(planet_x, planet_y, "*", "MarkerSize", 10);
            title(sprintf("Rotation: %.2f deg", rad2deg(angles(i))))
            axis equal
            hold off;

            figure(h2)
            subplot(L, L, j); hold on;
            for k = 1:size(data.instrument.positions, 1)
                plot(data.instrument.positions(k, 1), data.instrument.positions(k, 2), ".", "MarkerSize", 10)
            end
            title(sprintf("Rotation: %.2f deg", rad2deg(angles(i))))
            axis equal
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

% Modulation efficiency
efficiency = mean(modulation) / max(modulation);

% Plot the modulation curve
if isfield(data.outputs, "plot_modulation") && data.outputs.plot_modulation
    figure; hold on;
    plot(rad2deg(angles), modulation, '-', 'LineWidth', 2);

    if isfield(data.simulation, "monte_carlo_iterations") && data.simulation.monte_carlo_iterations > 0
        plot(rad2deg(angles), modulation_err, '-', 'LineWidth', 2);
    end

    xlabel('Rotation Angle [deg]');
    ylabel('Interpolated Signal at Planet Position');
    title('Planet Signal Modulation');
    grid minor;
end

end
