%% INTENSITY DISTRIBUTION
% Following pages 54-... of Viseur, L., the intensity of a 4 interferometer
% is derived, assuming a perfect signal). The procedure is generalisable to
% any N-apertures interferometer.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') % Change to NORMAL to export

data = ReadYaml('config/x_array.yml');
data = convert_data(data);

[optical_path, OPD, x_coords, y_coords] = load_opd(data.simulation.code_v_opd_file);

phases = opd2phase(optical_path, data.instrument.wavelength);

%%

perform_statistics(OPD * 1e6);

for i = 1:3
    h = plot_value_on_image_plane(OPD(:, i), x_coords(:, i), y_coords(:, i), title="OPD");
    h = plot_value_on_image_plane(phases(:, i), x_coords(:, i), y_coords(:, i), [], type="angles", title="Phase");
end