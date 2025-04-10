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

% Extract coordinates for plot
[~, ~, x_coords, y_coords] = load_opd("code_v/perturbed_21804_100sims_0804.txt");

% Computing of optimal splitting for analysis
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

% Define custom interval for this script
theta = linspace(-700, 700, 1000) * ((pi/180)/3600 / 1e3);

[RFs, RFn, R, ratio] = interferogram_sensitivity("code_v/nominal_100sim_21804points.txt", ...
    "code_v/perturbed_21804_100sims_0804.txt", data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.instrument.combination, data.instrument.wavelength);

% Representation of the results
maps_to_compute = 1 : floor(length(theta) / 4) : length(theta);

for i = 1:size(R, 3)
    label = sprintf("Intensity response at observation angle of %.2f mas", rad2mas(theta(maps_to_compute(i))));
    h = plot_value_on_image_plane(R(:, i, 1), x_coords(:, 1), y_coords(:, 1), title=label, type="_1e0");
end

h = plot_value_on_image_plane(ratio(:, 2), x_coords(:, 1), y_coords(:, 1), title="Nulling ratio", type="log");

RFp = [RFn; RFs];

h = plot_response_function_theta(theta, RFp, "Normalize", true);

%% Perturbate multiple branches with random picks

[RFs, RFn, R, ratio] = interferogram_sensitivity_multiple_branches("code_v/nominal_100sim_21804points.txt", ...
    "code_v/perturbed_21804_100sims_0804.txt", data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.instrument.combination, data.instrument.wavelength);

% Representation of the results
maps_to_compute = 1 : floor(length(theta) / 4) : length(theta);

for i = 1:size(R, 3)
    label = sprintf("Intensity response at observation angle of %.2f mas", rad2mas(theta(maps_to_compute(i))));
    h = plot_value_on_image_plane(R(:, i, 1), x_coords(:, 1), y_coords(:, 1), title=label, type="_1e0");
end

h = plot_value_on_image_plane(ratio(:, 2), x_coords(:, 1), y_coords(:, 1), title="Nulling ratio", type="log");

RFp = [RFn; RFs];

h = plot_response_function_theta(theta, RFp, "Normalize", true);