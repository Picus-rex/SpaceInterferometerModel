%% COMPLETE ANALYISIS
% Run the first section to load data, then every section that has a
% dependency on the desired aspect.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked')

% Change here the name of the file to run
data = ReadYaml('config/shifted_lin_array.yml');
data = convert_data(data);

[optical_path, OPD, x_coords, y_coords] = load_opd("code_v/perturbed_100sim_21804points.txt");
[optical_path_n, OPD_n, x_coords_n, y_coords_n] = load_opd("code_v/nominal_100sim_21804points.txt");

%% OPTIMAL SPLITTING AND BASELINES CLASSIFICATIONS

% Compute optimal splitting
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

% Classify baselines
[data.simulation.baselines, data.simulation.unique_baselines] = ...
    classify_baselines(data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, true);

%% TRANSMISSION MAPS

% Complex Field and Response Function
[Maps, data] = compute_response_function("data", data);

plot_transmission_map(data.simulation.theta_range, Maps.T_standard);

%% MODULATION

% Verify modulation of the signal
[data.simulation.planet_modulation] = planet_modulation(data);

%% NULLING RATIOS

% Find nulling ratio
[data.simulation.nulling_ratio, data.simulation.rejection_factor] = ...
    compute_nulling_ratio(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius, data.instrument.wavelength);

%% OPD AND NULLING RATIOS

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

%% PSF

data.instrument.phase_shifts = deg2rad([270, 180, 90, 0]);

[data.simulation.PSF, FWHM, C_i] = compute_psf(data.instrument.wavelength, ...
    data.instrument.intensities, data.instrument.positions, ...
    data.instrument.phase_shifts, ...
    data.environment.exoplanet_position, data.simulation.theta_range);

data.instrument.modulation_efficiency = compute_modulation_efficiency(...
        data.instrument.efficiencies.optical_line, ...
        data.instrument.surfaces, C_i);

[data.instrument.baselines, data.instrument.unique_baselines] = classify_baselines(data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, false);

data.instrument.FWHM = FWHM / (data.instrument.wavelength/data.instrument.baseline);

data.instrument.IWA = compute_IWA(data.instrument.unique_baselines, data.instrument.wavelength) / (data.instrument.wavelength/data.instrument.baseline);


plot_psf_map(data.simulation.theta_range, data.simulation.PSF);

%% Perturbation 

phases = opd2phase(optical_path, data.instrument.wavelength);
phases_n = opd2phase(optical_path_n, data.instrument.wavelength);

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

h = plot_value_on_image_plane(OPD(:, 1), x_coords(:, 1), y_coords(:, 1), title="OPD");
h = plot_value_on_image_plane(phases(:, 1), x_coords(:, 1), y_coords(:, 1), [], type="angles", title="Phase");

h = plot_value_on_image_plane(OPD_n(:, 1), x_coords_n(:, 1), y_coords_n(:, 1), title="OPD");
h = plot_value_on_image_plane(phases_n(:, 1), x_coords_n(:, 1), y_coords_n(:, 1), [], type="angles", title="Phase");

perform_statistics(OPD);

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
exotable = get_ppop_yield(0.5 * IWAs, OWA, ratios, data.instrument.wavelength, "verbose", true);

%% Interferogram arm (single)

% Define custom interval for this script
theta = linspace(-700, 700, 1000) * ((pi/180)/3600 / 1e3);

[RFs, RFn, R, ratio] = interferogram_sensitivity("code_v/nominal_100sim_21804points.txt", ...
    "code_v/perturbed_100sim_21804points.txt", data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.instrument.combination, data.instrument.wavelength);

% Representation of the results
maps_to_compute = 1 : floor(length(theta) / 4) : length(theta);

for i = 1:size(R, 3)
    label = sprintf("Intensity response at observation angle of %.2f mas", rad2mas(theta(maps_to_compute(i))));
    plot_value_on_image_plane(R(:, i, 1), x_coords(:, 1), y_coords(:, 1), title=label, type="_1e0");
end

plot_value_on_image_plane(ratio(:, 1), x_coords(:, 1), y_coords(:, 1), title="Nulling ratio", type="log");

RFp = [RFn; RFs];

h = plot_response_function_theta(theta, RFp, "Normalize", true);

% Mean, Median, Max across Ns simulations
log_mean = mean(ratio, 2);
log_median = median(ratio, 2);
log_max = max(ratio, [], 2); % worst nulling per point

% Plot them
plot_value_on_image_plane(log_mean, x_coords(:, 1), y_coords(:, 1), ...
    title="Mean Log Nulling Ratio", type="log");

plot_value_on_image_plane(log_median, x_coords(:, 1), y_coords(:, 1), ...
    title="Median Log Nulling Ratio", type="log");

plot_value_on_image_plane(log_max, x_coords(:, 1), y_coords(:, 1), ...
    title="Worst-case Log Nulling Ratio", type="log");

% Threshold for "good" nulling in log10 scale
threshold = 1e-4;
Ns = size(ratio, 2);
Np = size(ratio, 1);

% Count good cases for each pupil point
good_counts = sum(ratio < threshold, 2);
fraction_good = good_counts / Ns;

% Plot how frequently each position is "good"
plot_value_on_image_plane(fraction_good, x_coords(:, 1), y_coords(:, 1), ...
    title="Fraction of Good Nulling Simulations", type="linear");

% What % of the pupil is reliably good
good_pupil_fraction = sum(fraction_good > 0.9) / Np;  % e.g., 90% of sims are good
fprintf("Fraction of pupil with consistently good nulling: %.2f%%\n", 100 * good_pupil_fraction);

%% Interferogram arm multiple

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