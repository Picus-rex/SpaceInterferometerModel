%% TOLERANCING ANALYSIS
% Following pages 54-... of Viseur, L., the intensity of a 4 interferometer
% is derived, assuming a perfect signal). The procedure is generalisable to
% any N-apertures interferometer.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') 

% Loading of elements from configurations files
data = ReadYaml('config/linear_array.yml');
data = convert_data(data);

[optical_path, OPD, x_coords, y_coords] = load_opd("code_v/perturbed_0424_100sims.txt");
[optical_path_n, OPD_n, x_coords_n, y_coords_n] = load_opd("code_v/nominal_0424.txt");

phases = opd2phase(optical_path, data.instrument.wavelength);
phases_n = opd2phase(optical_path_n, data.instrument.wavelength);

% Computing of optimal splitting for analysis
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

% Computing of nulling ratios and statistical plots
[ratios, ~, ~]  = perturbate_system(data.instrument.apertures, ...
    data.instrument.intensities, data.instrument.phase_shifts, ...
    data.instrument.positions, data.environment.stellar_angular_radius, ...
    data.instrument.wavelength, phases, data.instrument.combination, ...
    perturbed_map_plotting_number=1, compute_PSF=false, create_plots=true);

% Computing of nulling ratios and statistical plots for nominal system
[ratios_n, ~, ~]  = perturbate_system(data.instrument.apertures, ...
    data.instrument.intensities, data.instrument.phase_shifts, ...
    data.instrument.positions, data.environment.stellar_angular_radius, ...
    data.instrument.wavelength, phases_n, data.instrument.combination, ...
    perturbed_map_plotting_number=1, compute_PSF=false, create_plots=false);

%% Create few plots and statistic

for i = 1:3
    h = plot_value_on_image_plane(OPD(:, i), x_coords(:, i), y_coords(:, i), title="OPD");
    h = plot_value_on_image_plane(phases(:, i), x_coords(:, i), y_coords(:, i), [], type="angles", title="Phase");
end

perform_statistics(OPD);

%% Exoplanet yield from NASA 

exotable2 = get_ppop_yield(IWAs, OWA, ratios, data.instrument.wavelength, "verbose", true, "population", "NASA");

% Get yield from PPOP for the nominal case
perturbed_vect = zeros(size(data.instrument.phase_shifts));
perturbed_vect(1) = rms(phases_n);
pp = data.instrument.phase_shifts + perturbed_vect;
aa = data.instrument.intensities .* data.instrument.combination;
[~, unique_baselines] = classify_baselines(aa, data.instrument.positions, pp, false);
IWA = compute_IWA(unique_baselines, data.instrument.wavelength);
exotable_n = get_ppop_yield(IWA, OWA, ratios_n, data.instrument.wavelength, "verbose", true);
