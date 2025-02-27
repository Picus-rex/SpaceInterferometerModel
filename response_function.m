%% Response Function for a 4-Aperture Nulling Interferometer
% Following the definition of the response function, this script performs
% the optimal splitting for the desired configuration and computes the
% fringes images.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') % Change to NORMAL to export

%% Parameters

N = 4;                                          % Apertures
lambda = 10e-6;                                 % Wavelength [m] 
B = 10;                                         % Baseline scaling factor
Ratio = 2;                                      % Apertures position ratio
stellar_ang_radius = 0.001 * lambda / B;        % Default value

plot_planets = 1;                               % 0: No 1: below 2: P-Pop

if plot_planets == 1
    theta_planets = [3e-6, 10e-6, 15e-6];       % Plan. Angular separations
    Np = length(theta_planets);                 % Simulated planets
elseif plot_planets == 2
    [ang_sep, stellar_ang_radius] = select_random_exoplanet("others/TestPlanetPopulation.txt");
end

plot_array = true;                              % Switch array plot
plot_star = true;                               % Plot star

theta_range = linspace(-20e-6, 20e-6, 2000);    % Angular grid [rad]
[theta_x, theta_y] = meshgrid(theta_range, theta_range);

% Positions (x, y) of the apertures 
positions_xarray = define_array("X-Array", B, Ratio);

positions_linear = define_array("Linear", B, Ratio*1.5);

positions = positions_xarray;

% Amplitudes for all apertures.
A = ones(4, 1);

% Plot array
if plot_array
    D = plot_apertures(positions, A, true);
end


% Compute optimal splitting
[~, combination, phase_shifts] = compute_optimal_splitting(N, B, lambda,...
    positions(:, 1), positions(:, 2), stellar_ang_radius, true);

classify_baselines(positions, phase_shifts, true);

%% Complex Field and Response Function

[T, T_chopped] = compute_response_function(lambda, N, ...
                positions, A, phase_shifts, combination, theta_x, theta_y);

%% Plot

angles = linspace(0, 2*pi, 10000); 

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
    star_x = stellar_ang_radius * cos(angles);
    star_y = stellar_ang_radius * sin(angles);
    plot(star_x, star_y, 'r--', 'LineWidth', 1.5);
end
xlabel('\theta_x (rad)');
ylabel('\theta_y (rad)');
title('Response Function');
colorbar; styling; % colormap winter;
axis xy; axis equal;

figure; 
semilogy(theta_range, T(floor(size(T, 1)/2), :), "LineWidth", 1.5, "DisplayName", "Normal response"); hold on;
semilogy(theta_range, T_chopped(floor(size(T, 1)/2), :), "LineWidth", 1.5, "DisplayName", "Phase chopping response");

% Scaling factor C for theta^2 and theta^4
C2 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^2;
C4 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^4;
C6 = T(floor(size(T, 1)/2), floor(size(T, 1)/2)) / theta_range(floor(size(T, 1)/2))^6;

semilogy(theta_range, C2 * theta_range.^(2), "--", 'DisplayName', "\theta^{2}", "LineWidth", 1.5);
semilogy(theta_range, C4 * theta_range.^(4), "--", 'DisplayName', "\theta^{4}", "LineWidth", 1.5);
semilogy(theta_range, C6 * theta_range.^(6), "--", 'DisplayName', "\theta^{6}", "LineWidth", 1.5);

xlabel('\theta_x (\lambda/B)');
ylabel('Normalized Intensity');
legend;
grid minor; hold off;

if plot_planets == 1
    figure;
    for p = 1:Np
        % Simulate planet motion
        theta_p = theta_planets(p);
        planet_response = NaN(size(theta_x));
        
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
        
        xlabel('\theta_x (rad)');
        ylabel('\theta_y (rad)');
        zlabel('Normalized Intensity');
        title(['Planet ', num2str(p)]);
        grid on; view([0, 0]);
    end
end
