%% TOLERANCING ANALYSIS
% Following pages 54-... of Viseur, L., the intensity of a 4 interferometer
% is derived, assuming a perfect signal). The procedure is generalisable to
% any N-apertures interferometer.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') 

data = ReadYaml('config/x_array.yml');
data = convert_data(data);

[optical_path, OPD, x_coords, y_coords] = load_opd(data.simulation.code_v_opd_file);

phases = opd2phase(optical_path, data.instrument.wavelength);

[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

%%

OPD_all = reshape(OPD, [], 1);

perform_statistics(OPD);
perform_statistics(OPD_all, create_plots=false);

%%

for i = 1:3
    h = plot_value_on_image_plane(OPD(:, i), x_coords(:, i), y_coords(:, i), title="OPD");
    h = plot_value_on_image_plane(phases(:, i), x_coords(:, i), y_coords(:, i), [], type="angles", title="Phase");
end

%%

[ratios, ~, PSFs]  = perturbate_system(data.instrument.apertures, ...
    data.instrument.intensities, data.instrument.phase_shifts, ...
    data.instrument.positions, data.environment.stellar_angular_radius, ...
    data.instrument.wavelength, phases, data.instrument.combination);

perform_statistics(ratios, "label", "Nulling ratio", "type", "_1e0", "scale", "log")