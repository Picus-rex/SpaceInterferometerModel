%% GENERATE_FIGURES This script generates all the figures within the thesis
% allowing for easiness of reproduction and substitution.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'normal') 

% Specify path of the git to export images following the convention adopted
path = "/Users/francesco/Library/CloudStorage/OneDrive-Personale/Università/SEMESTRE 12/TESI/Thesis/67a4aa2876545e0e504145e0/images/";
suffix = "_xarray";

warning("Always pull the last version from the system!")

export_data = ReadYaml('config/export_figures.yml');
data = load("exports/data_x_array.mat").data;

for figure = export_data.figures

    export_settings = struct("width", export_data.sizes.width.(figure{1}.width), ...
                "height", export_data.sizes.height.(figure{1}.height), ...
                "font_size", export_data.sizes.font_size, ...
                "name", path + export_data.chapters{figure{1}.chapter}.name + "/" + figure{1}.type + suffix);

    switch figure{1}.type

        case "transmission_map"
            h = plot_transmission_map(data.simulation.theta_range, data.simulation.T_standard, export_settings);
        
        case "transmission_map_planet"
            [~, h] = planet_modulation(data, export_settings);

        case "transmission_map_monodirectional"
            maps = struct();
            for i = 1:length(figure{1}.include)
                maps.(figure{1}.include{i}) = data.simulation.(figure{1}.include{i});
            end
            h = plot_transmission_map_monodirectional(data.simulation.theta_range, maps, figure{1}.names, export_settings);
        
        case "psf_map"
            h = plot_psf_map(data.simulation.theta_range, data.simulation.PSF, export_settings);

        case "optimal_modes"
            [~, ~, ~, h] = compute_optimal_splitting(...
                data.instrument.apertures, data.instrument.baseline, ...
                data.instrument.wavelength, data.instrument.positions(:, 1), ...
                data.instrument.positions(:, 2), ...
                data.environment.stellar_angular_radius, true, export_settings);
    end
end