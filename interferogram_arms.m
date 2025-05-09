%% INTERFEROGRAM ARMS
% See how the response depends on the phases on the x, y found rays
% positions from CODE V response. This script offers a way to see how the
% interfogram occurs on the lenses space and is auxiliary to the other
% processes.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') 

tic;

% Loading of elements from configurations files
data = convert_data('config/linear_array.yml');

% Overwrite the angular interval for this script and of maps to compute
% within the function.
theta = mas2rad(linspace(-100, 100, 2000));
maps_to_compute = 1 : floor(length(theta) / 4) : length(theta);

% Fix for X Array
%data.instrument.phase_shifts = [0, -pi, 0, -pi];
%data.instrument.combination = [0.2236    0.6708    0.6708    0.2236];

% Compute the response at the pupil screen with the defined function for
% the perturbed and corrected (if available) system
[data.outputs.response_function.perturbed, data.outputs.response_function.nominal, ...
 data.outputs.response_function.maps.perturbed, data.simulation.nulling_ratio_interferogram.perturbed, ...
 data.simulation.modulation_efficiency_interferogram.perturbed] = ...
    interferogram_sensitivity(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.perturbed.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, data.instrument.combination, ...
    data.instrument.surfaces, data.instrument.wavelength, theta, maps_to_compute);

[data.outputs.response_function.corrected, ~, data.outputs.response_function.maps.corrected, ...
    data.simulation.nulling_ratio_interferogram.corrected, ...
    data.simulation.modulation_efficiency_interferogram.corrected] = ...
    interferogram_sensitivity(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.corrected.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, data.instrument.combination, ...
    data.instrument.surfaces, data.instrument.wavelength, theta, maps_to_compute);

% The analysis has been moved to a specific function to call it with
% perturbed and corrected data from the compensator (if available)
plot_interferogram_results(data.outputs.response_function.nominal, data.outputs.response_function.perturbed, ...
    data.outputs.response_function.maps.perturbed, data.simulation.nulling_ratio_interferogram.perturbed, ...
    data.simulation.modulation_efficiency_interferogram.perturbed, data.simulation.code_v.nominal.x, ...
    data.simulation.code_v.perturbed.y, theta, maps_to_compute);
plot_interferogram_results(data.outputs.response_function.nominal, data.outputs.response_function.corrected, ...
    data.outputs.response_function.maps.corrected, data.simulation.nulling_ratio_interferogram.corrected, ...
    data.simulation.modulation_efficiency_interferogram.corrected, data.simulation.code_v.nominal.x, ...
    data.simulation.code_v.nominal.y, theta, maps_to_compute);

% Verify results with compensator.
compare_compensator(data.simulation.code_v.perturbed.opd, ... 
    data.simulation.nulling_ratio_interferogram.perturbed, data.simulation.code_v.corrected.opd, ...
    data.simulation.nulling_ratio_interferogram.corrected, data.simulation.code_v.corrected.x(:, 1), ...
    data.simulation.code_v.corrected.y(:, 1));

toc;

%% Perturbate multiple branches with random picks

[RFs, RFn, R, ratio] = interferogram_sensitivity_multiple_branches(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.perturbed.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.instrument.combination, data.instrument.wavelength, theta, maps_to_compute);

for i = 1:size(R, 3)
    label = sprintf("Intensity response at observation angle of %.2f mas", rad2mas(theta(maps_to_compute(i))));
    h = plot_value_on_image_plane(R(:, i, 1),  data.simulation.code_v.perturbed.x(:, 1),  data.simulation.code_v.perturbed.y(:, 1), title=label, type="_1e0");
end

h = plot_value_on_image_plane(ratio(:, 2),  data.simulation.code_v.perturbed.x(:, 1),  data.simulation.code_v.perturbed.y(:, 1), title="Nulling ratio", type="log");

RFp = [RFn; RFs];

h = plot_response_function_theta(theta, RFp, "Normalize", true);


