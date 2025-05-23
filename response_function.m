%% Response Function for a 4-Aperture Nulling Interferometer
% Following the definition of the response function, this script performs
% the optimal splitting for the desired configuration and computes the
% fringes images in the ideal case or with simulated disturbances.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked')

% Change here the name of the file to run
data = ReadYaml('config/linear_array.yml');
data = convert_data(data);

%% Analysis

% Compute optimal splitting
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, true);

% Classify baselines
[data.simulation.baselines, data.simulation.unique_baselines] = ...
    classify_baselines(data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, true);

% Complex Field and Response Function
[Maps, data] = compute_response_function("data", data);
T = Maps.T_standard;
T_chopped = Maps.T_chopped; 
T_real = Maps.T_real; 
T_real_chopped = Maps.T_real_chopped;

% Verify modulation of the signal
[data.simulation.planet_modulation] = planet_modulation(data);

% Find nulling ratio
[data.simulation.nulling_ratio, data.simulation.rejection_factor] = ...
    compute_nulling_ratio(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius, data.instrument.wavelength);

% Solve OPDs
find_minimum_OPD_single_branch(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius, data.instrument.wavelength);

ratio2OPD(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius);

OPDs2ratio(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius);

OPD2ratio_point(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius);

plot_phase_sensitivity(data.instrument.surfaces, data.instrument.intensities, ...
    data.instrument.positions, data.instrument.phase_shifts, ...
    data.instrument.efficiencies.optical_line, data.instrument.wavelength, ...
    data.environment.stellar_angular_radius);

%% Plot

plot_planets = data.outputs.plot_planets;
if plot_planets == 1
    theta_planets = cell2mat(data.outputs.planets_angular_separations);
    Np = length(theta_planets);                 % Simulated planets
elseif plot_planets == 2
    [ang_sep, stellar_ang_radius] = select_random_exoplanet("others/TestPlanetPopulation.txt");
end

plot_array = data.outputs.plot_array;
plot_star = data.outputs.plot_star;

angles = linspace(0, 2*pi, 10000); 
conversion_rad2mas = 1e3 * (3600 * 180) / pi;
theta_range = data.simulation.theta_range * conversion_rad2mas;
%theta_planets = theta_planets * conversion_rad2mas;

figure; hold on;
imagesc(theta_range, theta_range, T);
plot(0, 0, "r*", "LineWidth", 1.5)
if plot_planets == 1
    for p = 1:Np
        planet_x = theta_planets(p) * cos(angles);
        planet_y = theta_planets(p) * sin(angles);
        plot(planet_x, planet_y, 'w--', 'LineWidth', 1.5);
    end
elseif plot_planets == 2
    planet_x = ang_sep * cos(angles);
    planet_y = ang_sep * sin(angles);
    plot(planet_x, planet_y, 'w--', 'LineWidth', 1.5);
end
if plot_star
    star_x = data.environment.stellar_angular_radius * cos(angles);
    star_y = data.environment.stellar_angular_radius * sin(angles);
    plot(star_x, star_y, 'r--', 'LineWidth', 1.5);
end
xlabel('\theta_x [mas]');
ylabel('\theta_y [mas]');
title('Response Function');
colorbar; styling; % colormap winter;
axis xy; axis equal;

figure; 
semilogy(theta_range, T(floor(size(T, 1)/2), :), "LineWidth", 1.5, "DisplayName", "Normal response"); hold on;
semilogy(theta_range, T_chopped(floor(size(T, 1)/2), :), "LineWidth", 1.5, "DisplayName", "Phase chopping response");
if data.simulation.consider_non_ideal
    semilogy(theta_range, T_real(floor(size(T, 1)/2), :), "LineWidth", 1.5, "DisplayName", "Error response");
    semilogy(theta_range, T_real_chopped(floor(size(T, 1)/2), :), "LineWidth", 1.5, "DisplayName", "Error response chopped");
end

% Scaling factor C for theta^2 and theta^4
C2 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^2;
C4 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^4;
C6 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^6;

semilogy(theta_range, C2 * theta_range.^(2), "--", 'DisplayName', "\theta^{2}", "LineWidth", 1.5);
semilogy(theta_range, C4 * theta_range.^(4), "--", 'DisplayName', "\theta^{4}", "LineWidth", 1.5);
semilogy(theta_range, C6 * theta_range.^(6), "--", 'DisplayName', "\theta^{6}", "LineWidth", 1.5);

xlabel('\theta_x [mas]');
ylabel('Normalized Intensity');
legend;
grid minor; hold off;

if plot_planets == 1
    figure;
    for p = 1:Np
        % Simulate planet motion
        theta_p = theta_planets(p);
        planet_response = NaN(size(data.simulation.theta_x));
        
        x_points = [];
        y_points = [];
        intensity_values = [];
        
        for j = 1:length(angles)
            theta_xp = theta_p * cos(angles(j));
            theta_yp = theta_p * sin(angles(j));
            
            % Snap to grid to see response
            [~, i_x] = min(abs(theta_range - theta_xp));
            [~, i_y] = min(abs(theta_range - theta_yp));

            % Store values for scatter plot
            x_points = [x_points, theta_xp];
            y_points = [y_points, theta_yp];
            intensity_values = [intensity_values, T(i_y, i_x)];
        end

        % Plot response modulation for each planet
        subplot(1, Np, p);
        scatter3(x_points, y_points, intensity_values, 20, intensity_values, 'o', 'filled'); 
        styling; % colormap winter; 
        colorbar;
        
        xlabel('\theta_x [mas]');
        ylabel('\theta_y [mas]');
        zlabel('Normalized Intensity');
        title(['Planet ', num2str(p)]);
        grid on; view([0, 0]);
    end
end
