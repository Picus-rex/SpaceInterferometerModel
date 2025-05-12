%% COMPLETE ANALYISIS
% This script is used to generate the file that is loading in the image
% exporting files, therefore no visual output is provided. 

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked')

tic;

% CHANGE HERE THE ARRAY
data = convert_data('config/linear_array.yml');

type others/header.txt

% Compute optimal splitting
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

fprintf("\n✓ Optimal splitting\n")

% Classify baselines
[data.instrument.baselines, data.instrument.unique_baselines] = ...
    classify_baselines(data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, false);

fprintf("✓ Baseline classification\n")

% Complex Field and Response Function
[~, data] = compute_response_function("data", data);

fprintf("✓ Response function\n")

% Verify modulation of the signal
[data.simulation.planet_modulation] = planet_modulation(data);

fprintf("✓ Planet modulation\n")

% Find nulling ratio
[data.simulation.nulling_ratio, data.simulation.rejection_factor] = ...
    compute_nulling_ratio(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius, data.instrument.wavelength);

fprintf("✓ Nulling ratio\n")

% Point Spread Function
[data.simulation.PSF, ~, ~] = compute_psf(data.instrument.wavelength, ...
    data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, ...
    data.environment.exoplanet_position, data.simulation.theta_range);

fprintf("✓ Point spread function\n")

% Interferometer sensitivity - single arm
[data.outputs.response_function.perturbed, data.outputs.response_function.nominal, ...
 data.outputs.response_function.maps.perturbed, data.simulation.nulling_ratio_interferogram.perturbed, ...
 data.simulation.modulation_efficiency_interferogram.perturbed] = ...
    interferogram_sensitivity(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.perturbed.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, data.instrument.combination, ...
    data.instrument.surfaces, data.instrument.wavelength);

fprintf("✓ Interferometry sensitivity (perturbed)\n")

[data.outputs.response_function.corrected, ~, data.outputs.response_function.maps.corrected, ...
    data.simulation.nulling_ratio_interferogram.corrected, ...
    data.simulation.modulation_efficiency_interferogram.corrected] = ...
    interferogram_sensitivity(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.corrected.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, data.instrument.combination, ...
    data.instrument.surfaces, data.instrument.wavelength);

fprintf("✓ Interferometry sensitivity (corrected)\n")

% Interferometer sensitivity - multiple arms
[data.simulation.interferogram_multi.RFp, data.simulation.interferogram_multi.RFn, ...
 data.simulation.interferogram_multi.R, data.simulation.interferogram_multi.nulling] = ...
 interferogram_sensitivity_multiple_branches(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.corrected.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.instrument.combination, data.instrument.wavelength);

fprintf("✓ Interferometry sensitivity (multiple arms)\n")

% PSF and Transmission maps of perturbed systems
[data.simulation.perturbation.ratio, data.simulation.perturbation.transmission_maps, ...
    data.simulation.perturbation.nominal_map, data.simulation.perturbation.PSFs, ...
    data.simulation.perturbation.nominal_PSF]  = perturbate_system(data.instrument.apertures, ...
                    data.instrument.intensities, data.instrument.phase_shifts, ...
                    data.instrument.positions, data.environment.stellar_angular_radius, ...
                    data.instrument.wavelength, data.simulation.code_v.corrected.phase, data.instrument.combination, ...
                    perturbed_map_plotting_number=0, compute_PSF=true, create_plots=false);

fprintf("✓ PSF and transmission maps of perturbations\n")

% PPOP yield
data.instrument.IWA = data.instrument.wavelength / (2 * data.instrument.baseline);
data.instrument.OWA = data.instrument.wavelength / data.instrument.diameter(1);
[data.simulation.yield.ppop_table, data.simulation.yield.ppop_matrix] = get_ppop_yield(data.instrument.IWA, ...
    data.instrument.OWA, data.simulation.perturbation.ratio, data.instrument.wavelength, "create_plots", false);

fprintf("✓ PPOP yield\n")

save("exports/data_"+data.instrument.name+".mat", "data", "-v7.3");

fprintf("✓ Output saved \n")
toc;

%                      .
%                     / V\
%                   / `  /
%                  <<   |
%                  /    |
%                /      |
%              /        |
%            /    \  \ /
%           (      ) | |
%   ________|   _/_  | |
% <__________\______)\__)