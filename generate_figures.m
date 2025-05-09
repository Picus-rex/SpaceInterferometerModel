%% GENERATE_FIGURES This script generates all the figures within the thesis
% allowing for easiness of reproduction and substitution. This script
% relies on the the export_figures.yml configuration file that must be
% modified accordingly.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'normal') 

% Specify path of the git to export images following the convention adopted
export_data = ReadYaml('config/export_figures.yml');

% Verify if files exist
for i = 1:length(export_data.options.data_files)
    if ~isfile(export_data.options.data_files{i})
        error("File %s does not exist. Generate it using complete_analysis first!", export_data.options.data_files{i});
    end
end

% For every configuration to export...
for i = 1:length(export_data.options.data_files)

    data = load(export_data.options.data_files{i}).data;
    suffix = "_" + data.instrument.name;

    theta = mas2rad(linspace(-100, 100, 1000));
    maps_to_compute = 1 : floor(length(theta) / 3) : length(theta);
    
    for figure = export_data.figures

        cur_fig = figure{1};
    
        if ismember(cur_fig.chapter, cell2mat(export_data.options.skip_chapters))
            continue
        end
        
        label_name = export_data.chapters{cur_fig.chapter}.path + "images/" + export_data.chapters{cur_fig.chapter}.name + "/" + cur_fig.type + suffix;

        export_settings = struct("width", export_data.sizes.width.(cur_fig.width), ...
                    "height", export_data.sizes.height.(cur_fig.height), ...
                    "font_size", export_data.sizes.font_size, ...
                    "name", label_name);
    
        switch cur_fig.type
            
            case "configuration"
                [~, ~, h] = plot_apertures(data.instrument.positions, data.instrument.intensities, data.instrument.efficiencies.optical_line, data.instrument.efficiencies.optical_line, export_settings);

            case "transmission_map"
                if cur_fig.planets
                    theta_planets = [300, 750, 1250];
                    h = plot_transmission_map(data.simulation.theta_range, data.simulation.T_standard, export_settings, theta_planets);
                else
                    h = plot_transmission_map(data.simulation.theta_range, data.simulation.T_standard, export_settings);
                end
            
            case "transmission_map_planet"
                export_settings.include = cur_fig.include;
                export_settings.sizes = export_data.sizes;
                [~, h] = planet_modulation(data, export_settings);
    
            case "transmission_map_monodirectional"
                maps = struct();
                for j = 1:length(cur_fig.include)
                    maps.(cur_fig.include{j}) = data.simulation.(cur_fig.include{j});
                end
                h = plot_transmission_map_monodirectional(data.simulation.theta_range, maps, cur_fig.names, export_settings);
            
            case "psf_map"
                h = plot_psf_map(data.simulation.theta_range, data.simulation.PSF, export_settings);
    
            case "optimal_modes"
                export_settings.include = cur_fig.include;
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
                    data.environment.stellar_angular_radius, cell2mat(cur_fig.wavelengths), export_settings);
            
            case "nullratio_double_branch"
                OPD2ratio_point(data.instrument.apertures, data.instrument.intensities, ...
                    data.instrument.phase_shifts, data.instrument.positions, ...
                    data.environment.stellar_angular_radius, ...
                    linspace(1e-6, 1e-4, 1000), [2; 1] * 1e-7, export_settings);
            
            case "opd_pupil_plane_system"
                export_settings.name = label_name + "_nominal";
                plot_value_on_image_plane(data.simulation.code_v.nominal.opd(:, 1), data.simulation.code_v.nominal.x(:, 1), data.simulation.code_v.nominal.y(:, 1), title="OPD", embedded=export_settings);
                
                export_settings.name = label_name + "_pert1";
                plot_value_on_image_plane(data.simulation.code_v.perturbed.opd(:, 1), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), title="OPD", embedded=export_settings);
                export_settings.name = label_name + "_pert2";
                plot_value_on_image_plane(data.simulation.code_v.perturbed.opd(:, 2), data.simulation.code_v.perturbed.x(:, 2), data.simulation.code_v.perturbed.y(:, 2), title="OPD", embedded=export_settings);
            
            case "phase_pupil_plane_system"     
                export_settings.name = label_name + "_nominal";
                plot_value_on_image_plane(data.simulation.code_v.nominal.phase(:, 1), data.simulation.code_v.nominal.x(:, 1), data.simulation.code_v.nominal.y(:, 1), [], type="angles", title="Phase", embedded=export_settings);
                
                export_settings.name = label_name + "_pert1";
                plot_value_on_image_plane(data.simulation.code_v.perturbed.phase(:, 1), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), [], type="angles", title="Phase", embedded=export_settings);
                export_settings.name = label_name + "_pert2";
                plot_value_on_image_plane(data.simulation.code_v.perturbed.phase(:, 2), data.simulation.code_v.perturbed.x(:, 2), data.simulation.code_v.perturbed.y(:, 2), [], type="angles", title="Phase", embedded=export_settings);
            
            case "sensitivity_opd_system"
                
                % Prepare image export
                exp_figs = {[], []};
                if ismember("fitted_gaussian", cur_fig.include)
                    exp_figs{1} = export_settings;
                    exp_figs{1}.name = label_name + "_fitted_gaussian";
                end
                if ismember("rms_vs_std", cur_fig.include)
                    exp_figs{2} = export_settings;
                    exp_figs{2}.name = label_name + "_rms_std";
                end

                perform_statistics(data.simulation.code_v.perturbed.opd, "embedded", exp_figs);
            
            case "interferogram"
                
                if isfield(cur_fig.include, "interferometer_response_angles")
                    for j = 1:size(data.simulation.interferogram.R, 3)
                        export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_response_angles.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_response_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_angles_%d", j));
                        plot_value_on_image_plane(data.outputs.response_function.maps.(cur_fig.simulation)(:, j, 1), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="_1e0", embedded=export_settings);
                    end
                end

                if isfield(cur_fig.include, "interferometer_modulation_angles")
                    for j = 1:3
                        export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_nulling_angles.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_nulling_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_modulation_%d", j));
                        plot_value_on_image_plane(data.simulation.modulation_efficiency_interferogram.(cur_fig.simulation)(:, j), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="_1e0", embedded=export_settings);
                    end
                end

                if isfield(cur_fig.include, "interferometer_nulling_angles")
                    for j = 1:3
                        export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_nulling_angles.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_nulling_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_nulling_%d", j));
                        plot_value_on_image_plane(data.simulation.nulling_ratio_interferogram.(cur_fig.simulation)(:, j), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);
                    end
                end
                
                if isfield(cur_fig.include, "interferometer_nulling_statistics")
                    log_mean = mean(data.simulation.nulling_ratio_interferogram.(cur_fig.simulation), 2);
                    log_median = median(data.simulation.nulling_ratio_interferogram.(cur_fig.simulation), 2);
                    log_max = max(data.simulation.nulling_ratio_interferogram.(cur_fig.simulation), [], 2);
                    fraction_good = sum(data.simulation.nulling_ratio_interferogram.(cur_fig.simulation) < 1e-4, 2) / size(data.simulation.nulling_ratio_interferogram.(cur_fig.simulation), 2);

                    export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_nulling_statistics.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_nulling_statistics.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name);

                    export_settings.name = label_name + sprintf("_nullingmean_%d", j);
                    plot_value_on_image_plane(log_mean, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);

                    export_settings.name = label_name + sprintf("_nullingmedian_%d", j);
                    plot_value_on_image_plane(log_median, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);

                    export_settings.name = label_name + sprintf("_nullingmax_%d", j);
                    plot_value_on_image_plane(log_max, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);

                    export_settings.name = label_name + sprintf("_nullinggood_%d", j);
                    plot_value_on_image_plane(fraction_good, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="linear", embedded=export_settings);
                end


                if isfield(cur_fig.include, "interferometer_response_function")
                    RFp = [data.outputs.response_function.(cur_fig.simulation); data.outputs.response_function.nominal];
                    export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_response_function.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_response_function.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_rf_%d", j));
                    plot_response_function_theta(theta, RFp, "Normalize", true, "embedded", export_settings);

                end

            case "interferogram_multiple"
                
                if isfield(cur_fig.include, "interferometer_response_angles")
                    for j = 1:size(data.simulation.interferogram_multi.R, 3)
                        export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_response_angles.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_response_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_angles_%d", j));
                        plot_value_on_image_plane(data.simulation.interferogram_multi.R(:, j, 1), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="_1e0", embedded=export_settings);
                    end
                end

                if isfield(cur_fig.include, "interferometer_nulling_angles")
                    for j = 1:3
                        export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_nulling_angles.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_nulling_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_nulling_%d", j));
                        plot_value_on_image_plane(data.simulation.interferogram_multi.nulling(:, j), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);
                    end
                end
                
                if isfield(cur_fig.include, "interferometer_nulling_statistics")
                    log_mean = mean(data.simulation.interferogram_multi.nulling, 2);
                    log_median = median(data.simulation.interferogram_multi.nulling, 2);
                    log_max = max(data.simulation.interferogram_multi.nulling, [], 2);
                    fraction_good = sum(data.simulation.interferogram_multi.nulling < 1e-4, 2) / size(data.simulation.interferogram_multi.nulling, 2);

                    export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_nulling_statistics.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_nulling_statistics.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name);

                    export_settings.name = label_name + sprintf("_nullingmean_%d", j);
                    plot_value_on_image_plane(log_mean, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);

                    export_settings.name = label_name + sprintf("_nullingmedian_%d", j);
                    plot_value_on_image_plane(log_median, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);

                    export_settings.name = label_name + sprintf("_nullingmax_%d", j);
                    plot_value_on_image_plane(log_max, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);

                    export_settings.name = label_name + sprintf("_nullinggood_%d", j);
                    plot_value_on_image_plane(fraction_good, data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="linear", embedded=export_settings);
                end


                if isfield(cur_fig.include, "interferometer_response_function")
                    RFp = [data.simulation.interferogram_multi.RFn; data.simulation.interferogram.RFp];
                    export_settings = struct("width", export_data.sizes.width.(cur_fig.include.interferometer_response_function.width), ...
                                "height", export_data.sizes.height.(cur_fig.include.interferometer_response_function.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_rf_%d", j));
                    plot_response_function_theta(theta, RFp, "Normalize", true, "embedded", export_settings);

                end

            case "transmission_map_perturbed"
                
                exp_figures = {};

                for j = 1:length(cur_fig.include)
                    if strcmp(cur_fig.include{j}, "STD")
                        exp_figures{1} = export_settings;
                        exp_figures{1}.name = label_name + "_STD";
                    elseif strcmp(cur_fig.include{j}, "CDF")
                        exp_figures{2} = export_settings;
                        exp_figures{2}.name = label_name + "_CDF";
                    elseif strcmp(cur_fig.include{j}, "PCA")
                        exp_figures{3} = export_settings;
                        exp_figures{3}.name = label_name + "_PCA";
                    elseif strcmp(cur_fig.include{j}, "PC1")
                        exp_figures{4} = export_settings;
                        exp_figures{4}.name = label_name + "_PC1";
                    end
                end

                plot_variation_map(size(data.simulation.code_v.perturbed.op, 2), ...
                    load_grid("x"), load_grid("y"), data.simulation.perturbation.transmission_maps, ...
                    data.simulation.perturbation.nominal_map, exp_figures);

             case "psf_perturbed"
                
                exp_figures = {};

                for j = 1:length(cur_fig.include)
                    if strcmp(cur_fig.include{j}, "STD")
                        exp_figures{1} = export_settings;
                        exp_figures{1}.name = label_name + "_STD";
                    elseif strcmp(cur_fig.include{j}, "CDF")
                        exp_figures{2} = export_settings;
                        exp_figures{2}.name = label_name + "_CDF";
                    elseif strcmp(cur_fig.include{j}, "PCA")
                        exp_figures{3} = export_settings;
                        exp_figures{3}.name = label_name + "_PCA";
                    elseif strcmp(cur_fig.include{j}, "PC1")
                        exp_figures{4} = export_settings;
                        exp_figures{4}.name = label_name + "_PC1";
                    end
                end

                plot_variation_map(size(data.simulation.code_v.perturbed.op, 2), ...
                    load_grid("x"), load_grid("y"), data.simulation.perturbation.PSFs, ...
                    data.simulation.perturbation.nominal_PSF, exp_figures);

            case "ppop_yield"
                exp_figures = {};
                exp_figures.figures.yields = export_settings;
                exp_figures.figures.yields.name = label_name + "_yield";
                exp_figures.figures.mean_yield = [];
                exp_figures.figures.total_yield = export_settings;
                exp_figures.figures.total_yield.name = label_name + "_total_yield";
                exp_figures.figures.boxplot = [];
                exp_figures.figures.trend = [];
                exp_figures.universes_to_plot = 1;

                plot_ppop_yield(data.simulation.yield.ppop_table, data.instrument.IWA, ...
                    data.instrument.OWA, min(rms(data.simulation.perturbation.ratio)), exp_figures);
        end
    end

    close all;

end