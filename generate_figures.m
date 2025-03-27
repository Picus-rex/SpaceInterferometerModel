%% GENERATE_FIGURES This script generates all the figures within the thesis
% allowing for easiness of reproduction and substitution.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'normal') 
warning("Always pull the last version from the system!")

% Specify path of the git to export images following the convention adopted
path = "/Users/francesco/Library/CloudStorage/OneDrive-Personale/Università/SEMESTRE 12/TESI/Thesis/67a4aa2876545e0e504145e0/";
export_data = ReadYaml('config/export_figures.yml');
data_matrices = ["exports/data_linear", "exports/data_xarray"];
skip_chapt = [2];

for i = 1:length(data_matrices)

    datas(i) = load(data_matrices(i)).data;
    data = datas(i);
    suffix = "_" + data.instrument.name;
    
    for figure = export_data.figures

        if ismember(figure{1}.chapter, skip_chapt) 
            continue
        end
    
        export_settings = struct("width", export_data.sizes.width.(figure{1}.width), ...
                    "height", export_data.sizes.height.(figure{1}.height), ...
                    "font_size", export_data.sizes.font_size, ...
                    "name", path + "images/" + export_data.chapters{figure{1}.chapter}.name + "/" + figure{1}.type + suffix);
    
        switch figure{1}.type
            
            case "configuration"
                [~, ~, h] = plot_apertures(data.instrument.positions, data.instrument.intensities, data.instrument.efficiencies.optical_line, data.instrument.efficiencies.optical_line, export_settings);

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
                export_settings.include = figure{1}.include;
                export_settings.sizes = export_data.sizes;
                [~, ~, ~, h] = compute_optimal_splitting(...
                    data.instrument.apertures, data.instrument.baseline, ...
                    data.instrument.wavelength, data.instrument.positions(:, 1), ...
                    data.instrument.positions(:, 2), ...
                    data.environment.stellar_angular_radius, true, export_settings);
            
            case "nullratio_single_branch_single_wavelenth"
                min_ratio = compute_nulling_ratio(data.instrument.apertures, ...
                    data.instrument.intensities, data.instrument.phase_shifts, ...
                    data.instrument.positions, data.environment.stellar_angular_radius, ...
                    data.instrument.wavelength);

                find_minimum_OPD_single_branch(data.instrument.apertures, data.instrument.intensities, ...
                    data.instrument.phase_shifts, data.instrument.positions, ...
                    data.environment.stellar_angular_radius, data.instrument.wavelength, ...
                    logspace(log10(min_ratio), 0, 10), export_settings);

            case "nullratio_single_branch"
                ratio2OPD(data.instrument.apertures, data.instrument.intensities, ...
                    data.instrument.phase_shifts, data.instrument.positions, ...
                    data.environment.stellar_angular_radius, logspace(-11, 0, 100), ...
                    linspace(1e-6, 20e-6, 100), export_settings);

            case "nullratio_double_branch_discrete"
                OPDs2ratio(data.instrument.apertures, data.instrument.intensities, ...
                    data.instrument.phase_shifts, data.instrument.positions, ...
                    data.environment.stellar_angular_radius, cell2mat(figure{1}.wavelengths), export_settings);
            
            case "nullratio_double_branch"
                OPD2ratio_point(data.instrument.apertures, data.instrument.intensities, ...
                    data.instrument.phase_shifts, data.instrument.positions, ...
                    data.environment.stellar_angular_radius, ...
                    linspace(1e-6, 1e-4, 1000), [2; 1] * 1e-7, export_settings);

        end
    end

end


for table = export_data.tables

    switch table{1}.type
            case "unique_baselines"
                filename = path + "tables/" + export_data.chapters{table{1}.chapter}.name + "/" + table{1}.type;
                export_baselines_table(datas(2).instrument.unique_baselines, ...
                    datas(1).instrument.unique_baselines, filename)
    end

end