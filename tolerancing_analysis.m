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

[optical_path, OPD, x_coords, y_coords] = load_opd("code_v/perturbed_100sim_1961points.txt");
[optical_path_n, OPD_n, x_coords_n, y_coords_n] = load_opd("code_v/nominal_100sim_1961points.txt");

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

% Computing of nulling ratios and statistical plots
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

OPD_all = reshape(OPD, [], 1);

perform_statistics(OPD);
perform_statistics(OPD_all, create_plots=false);
perform_statistics(ratios, "label", "Nulling ratio", "type", "_1e0", "scale", "log")

%% Verify exoplanet yield

Ns = size(ratios, 2);

% Allocate space
perturbed_vect = zeros(size(data.instrument.phase_shifts));
IWAs = zeros(Ns, 1);

% For every simulation compute the IWA
for i = 1:Ns
    perturbed_vect(1) = rms(phases(:, i));
    pp = data.instrument.phase_shifts + perturbed_vect;
    aa = data.instrument.intensities .* data.instrument.combination;
    [~, unique_baselines] = classify_baselines(aa, data.instrument.positions, pp, false);

    IWAs(i) = compute_IWA(unique_baselines, data.instrument.wavelength);
end

% Use simplified relation for the OWA
OWA = data.instrument.wavelength / data.instrument.diameter(1);

% Get yield from PPOP
% exotable = get_ppop_yield(IWAs, OWA, ratios, data.instrument.wavelength, "verbose", true);
exotable2 = get_ppop_yield(IWAs, OWA, ratios, data.instrument.wavelength, "verbose", true, "population", "NASA");

% Get yield from PPOP for the nominal case
perturbed_vect = zeros(size(data.instrument.phase_shifts));
perturbed_vect(1) = rms(phases_n);
pp = data.instrument.phase_shifts + perturbed_vect;
aa = data.instrument.intensities .* data.instrument.combination;
[~, unique_baselines] = classify_baselines(aa, data.instrument.positions, pp, false);
IWA = compute_IWA(unique_baselines, data.instrument.wavelength);
exotable_n = get_ppop_yield(IWA, OWA, ratios_n, data.instrument.wavelength, "verbose", true);

%% Improved cases

optical_path_b = optical_path(:, [44, 62, 99]);
optical_path_b(:, 4) = optical_path_n;

phases_b = opd2phase(optical_path_b, data.instrument.wavelength);

[ratios_b, ~, ~]  = perturbate_system(data.instrument.apertures, ...
    data.instrument.intensities, data.instrument.phase_shifts, ...
    data.instrument.positions, data.environment.stellar_angular_radius, ...
    data.instrument.wavelength, phases_b, data.instrument.combination, perturbed_map_plotting_number=4, compute_PSF=false, create_plots=true);

perform_statistics(ratios_b, "label", "Nulling ratio", "type", "_1e0", "scale", "log")
