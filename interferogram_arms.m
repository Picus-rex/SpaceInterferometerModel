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
 data.outputs.response_function.maps, data.simulation.nulling_ratio_interferogram, ...
 data.simulation.modulation_efficiency_interferogram] = ...
    interferogram_sensitivity(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.perturbed.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, data.instrument.combination, ...
    data.instrument.surfaces, data.instrument.wavelength, theta, maps_to_compute);

[data.outputs.response_function.corrected, ~, data.outputs.response_function.maps_corrected, ...
    data.simulation.nulling_ratio_interferogram_corrected, ...
    data.simulation.modulation_efficiency_interferogram_corrected] = ...
    interferogram_sensitivity(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.corrected.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, data.instrument.combination, ...
    data.instrument.surfaces, data.instrument.wavelength, theta, maps_to_compute);

% The analysis has been moved to a specific function to call it with
% perturbed and corrected data from the compensator (if available)
interferogram_analysis(data, data.simulation.code_v.perturbed, theta, maps_to_compute)
interferogram_analysis(data, data.simulation.code_v.corrected, theta, maps_to_compute)

% Verify results with compensator.
compare_compensator(data.simulation.code_v.perturbed.opd, ... 
    data.simulation.nulling_ratio_interferogram, data.simulation.code_v.corrected.opd, ...
    data.simulation.nulling_ratio_interferogram_corrected, data.simulation.code_v.corrected.x(:, 1), ...
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

%% FUNCTIONS

function interferogram_analysis(data, perturb_data, theta, maps_to_compute)

    % For all results requested with maps_to_compute, plots them. 
    R = data.outputs.response_function.maps;
    for i = 1:size(R, 3)
        label = sprintf("Intensity response at observation angle of %.2f mas", rad2mas(theta(maps_to_compute(i))));
        plot_value_on_image_plane(R(:, i, 1), perturb_data.x(:, 1), perturb_data.y(:, 1), title=label, type="_1e0");
    end
    
    % Group results for plotting the reponse function
    RFp = [data.outputs.response_function.nominal; data.outputs.response_function.perturbed];
    plot_response_function_theta(theta, RFp, "Normalize", true);
    
    % Plot also the nulling ratio (which depends only on the simulation)
    ratio = data.simulation.nulling_ratio_interferogram;
    plot_value_on_image_plane(ratio(:, 1), perturb_data.x(:, 1), perturb_data.y(:, 1), title="Nulling ratio", type="log");
    
    % Plot also the modulation efficieny
    eff = data.simulation.modulation_efficiency_interferogram;
    plot_value_on_image_plane(eff(:, 1), perturb_data.x(:, 1), perturb_data.y(:, 1), title="Modulation efficiency", type="_1e0");
    
    % Mean, Median, Max across Ns simulations
    log_mean = mean(ratio, 2);
    eff_mean = mean(eff, 2);
    log_median = median(ratio, 2);
    log_max = max(ratio, [], 2); 
    
    % eff_mean = mean(eff, 2);
    plot_value_on_image_plane(eff_mean, perturb_data.x(:, 1), perturb_data.y(:, 1), ...
        title="Mean Modulation Efficiency", type="1e0");
    
    plot_value_on_image_plane(log_mean, perturb_data.x(:, 1), perturb_data.y(:, 1), ...
        title="Mean Log Nulling Ratio", type="log");
    
    plot_value_on_image_plane(log_median, perturb_data.x(:, 1), perturb_data.y(:, 1), ...
        title="Median Log Nulling Ratio", type="log");
    
    plot_value_on_image_plane(log_max, perturb_data.x(:, 1), perturb_data.y(:, 1), ...
        title="Worst-case Log Nulling Ratio", type="log");
    
    % Threshold for "good" nulling in log10 scale
    threshold = 1e-4;
    Ns = size(ratio, 2);
    Np = size(ratio, 1);
    
    % Count good cases for each pupil point
    good_counts = sum(ratio < threshold, 2);
    fraction_good = good_counts / Ns;
    
    % Plot how frequently each position is "good"
    plot_value_on_image_plane(fraction_good, perturb_data.x(:, 1), perturb_data.y(:, 1), ...
        title="Fraction of Good Nulling Simulations", type="linear");
    
    % What % of the pupil is reliably good
    good_pupil_fraction = sum(fraction_good > 0.9) / Np;  % e.g., 90% of sims are good
    fprintf("Fraction of pupil with consistently good nulling: %.2f%%\n", 100 * good_pupil_fraction);

end
