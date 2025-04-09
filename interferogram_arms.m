%% INTERFEROGRAM ARMS
% See how the response depends on the phases on the x, y found rays
% positions from CODE V response. This script offers a way to see how the
% interfogram occurs on the lenses space and is auxiliary to the other
% processes.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') 

% Loading of elements from configurations files
data = ReadYaml('config/linear_array.yml');
data = convert_data(data);

[optical_path, OPD, x_coords, y_coords] = load_opd("code_v/perturbed_21804_100sims_0804.txt");
[optical_path_n, OPD_n, x_coords_n, y_coords_n] = load_opd("code_v/nominal_100sim_21804points.txt");

delta_phi = opd2phase(OPD(:, 1), data.instrument.wavelength);
delta_phi_n = opd2phase(OPD_n(:, 1), data.instrument.wavelength);

% Computing of optimal splitting for analysis
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

% Define custom interval for this script
theta = linspace(-700, 700, 1000) * ((pi/180)/3600 / 1e3);

% Allocate space for perturbations
phase_perturbations = zeros(length(x_coords), length(theta));

% Perturb first arm, leave others unchanged
for k = 1:data.instrument.apertures
    if k == 1
        phase_perturbations(:, k) = data.instrument.phase_shifts(k) + delta_phi;
    else
        phase_perturbations(:, k) = data.instrument.phase_shifts(k) + delta_phi_n;
    end
end

% Compute the response function for all the points in the analysis
[R, RF] = create_interferogram(data.instrument.positions, ...
    data.instrument.intensities, data.instrument.combination, ...
    phase_perturbations, data.instrument.wavelength, theta);

% Compute nulling ratio for each point on the screen
ratio = compute_nulling_ratio(data.instrument.apertures, ...
    data.instrument.intensities, phase_perturbations, ...
    data.instrument.positions, data.environment.stellar_angular_radius, ...
    data.instrument.wavelength);

%% Representation of the results

for i = [1, 10, 100, 200, 250]
    label = sprintf("Intensity response at observation angle of %.2f mas", rad2mas(theta(i)));
    h = plot_value_on_image_plane(R(:, i), x_coords(:, 1), y_coords(:, 1), title=label, type="_1e0");
end

h = plot_value_on_image_plane(ratio, x_coords(:, 1), y_coords(:, 1), title="Nulling ratio", type="log");

h = plot_response_function_theta(theta, RF);