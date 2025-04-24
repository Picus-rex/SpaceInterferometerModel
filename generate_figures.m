%% GENERATE_FIGURES This script generates all the figures within the thesis
% allowing for easiness of reproduction and substitution. This script
% relies on the the export_figures.yml configuration file that must be
% modified accordingly. 

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'normal') 

%% PART 1: GENERATE FILES FOR THE ANALYSIS

data = convert_data('config/x_array.yml');

% Compute optimal splitting
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

% Classify baselines
[data.instrument.baselines, data.instrument.unique_baselines] = ...
    classify_baselines(data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, true);

% Complex Field and Response Function
[~, data] = compute_response_function("data", data);

% Verify modulation of the signal
[data.simulation.planet_modulation] = planet_modulation(data);

% Find nulling ratio
[data.simulation.nulling_ratio, data.simulation.rejection_factor] = ...
    compute_nulling_ratio(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius, data.instrument.wavelength);

% Point Spread Function
[data.simulation.PSF, ~, ~] = compute_psf(data.instrument.wavelength, ...
    data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, ...
    data.environment.exoplanet_position, data.simulation.theta_range);

% Interferometer sensitivity - single arm
[data.simulation.interferogram.RFp, data.simulation.interferogram.RFn, ...
 data.simulation.interferogram.R, data.simulation.interferogram.nulling] = ...
    interferogram_sensitivity(data.simulation.code_v.nominal.opd, data.simulation.code_v.perturbed.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.instrument.combination, data.instrument.wavelength);

% Interferometer sensitivity - multiple arms
[data.simulation.interferogram_multi.RFp, data.simulation.interferogram_multi.RFn, ...
 data.simulation.interferogram_multi.R, data.simulation.interferogram_multi.nulling] = ...
 interferogram_sensitivity_multiple_branches(data.simulation.code_v.nominal.opd, ...
    data.simulation.code_v.perturbed.opd, data.instrument.intensities,...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.instrument.combination, data.instrument.wavelength);

clc;
save("exports/data_"+data.instrument.name+".mat", "data", "-v7.3");

%% PART 2: LOAD AND EXPORT IMAGES

warning("Always pull the last version from the system!")

% Specify path of the git to export images following the convention adopted
export_data = ReadYaml('config/export_figures.yml');
data_matrices = ["exports/data_linear", "exports/data_xarray"];
skip_chapt = [1, 2, 3];

for i = 1:length(data_matrices)

    datas(i) = load(data_matrices(i)).data;
    data = datas(i);
    suffix = "_" + data.instrument.name;

    theta = mas2rad(linspace(-100, 100, 1000));
    maps_to_compute = 1 : floor(length(theta) / 4) : length(theta);
    
    for figure = export_data.figures

        if ismember(figure{1}.chapter, skip_chapt) 
            continue
        end
        
        label_name = export_data.chapters{figure{1}.chapter}.path + "images/" + export_data.chapters{figure{1}.chapter}.name + "/" + figure{1}.type + suffix;

        export_settings = struct("width", export_data.sizes.width.(figure{1}.width), ...
                    "height", export_data.sizes.height.(figure{1}.height), ...
                    "font_size", export_data.sizes.font_size, ...
                    "name", label_name);
    
        switch figure{1}.type
            
            case "configuration"
                [~, ~, h] = plot_apertures(data.instrument.positions, data.instrument.intensities, data.instrument.efficiencies.optical_line, data.instrument.efficiencies.optical_line, export_settings);

            case "transmission_map"
                h = plot_transmission_map(data.simulation.theta_range, data.simulation.T_standard, export_settings);
            
            case "transmission_map_planet"
                export_settings.include = figure{1}.include;
                export_settings.sizes = export_data.sizes;
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
                exp_figs = {[], [], [], [], [], []};
                if ismember("fitted_gaussian", figure{1}.include)
                    exp_figs{2} = export_settings;
                    exp_figs{2}.name = label_name + "_fitted_gaussian";
                end
                if ismember("rms_vs_std", figure{1}.include)
                    exp_figs{5} = export_settings;
                    exp_figs{5}.name = label_name + "_rms_std";
                end

                perform_statistics(data.simulation.code_v.perturbed.opd, "embedded", exp_figs);
            
            case "interferogram"
                
                if isfield(figure{1}.include, "interferometer_response_angles")
                    for j = 1:size(data.simulation.interferogram.R, 3)
                        export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_response_angles.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_response_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_angles_%d", j));
                        plot_value_on_image_plane(data.simulation.interferogram.R(:, j, 1), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="_1e0", embedded=export_settings);
                    end
                end

                if isfield(figure{1}.include, "interferometer_nulling_angles")
                    for j = 1:3
                        export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_nulling_angles.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_nulling_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_nulling_%d", j));
                        plot_value_on_image_plane(data.simulation.interferogram.nulling(:, j), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);
                    end
                end
                
                if isfield(figure{1}.include, "interferometer_nulling_statistics")
                    log_mean = mean(data.simulation.interferogram.nulling, 2);
                    log_median = median(data.simulation.interferogram.nulling, 2);
                    log_max = max(data.simulation.interferogram.nulling, [], 2); 
                    fraction_good = sum(data.simulation.interferogram.nulling < 1e-4, 2) / size(data.simulation.interferogram.nulling, 2);

                    export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_nulling_statistics.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_nulling_statistics.height), ...
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


                if isfield(figure{1}.include, "interferometer_response_function")
                    RFp = [data.simulation.interferogram.RFn; data.simulation.interferogram.RFp];
                    export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_response_function.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_response_function.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_rf_%d", j));
                    plot_response_function_theta(theta, RFp, "Normalize", true, "embedded", export_settings);

                end

            case "interferogram_multiple"
                
                if isfield(figure{1}.include, "interferometer_response_angles")
                    for j = 1:size(data.simulation.interferogram_multi.R, 3)
                        export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_response_angles.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_response_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_angles_%d", j));
                        plot_value_on_image_plane(data.simulation.interferogram_multi.R(:, j, 1), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="_1e0", embedded=export_settings);
                    end
                end

                if isfield(figure{1}.include, "interferometer_nulling_angles")
                    for j = 1:3
                        export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_nulling_angles.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_nulling_angles.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_nulling_%d", j));
                        plot_value_on_image_plane(data.simulation.interferogram_multi.nulling(:, j), data.simulation.code_v.perturbed.x(:, 1), data.simulation.code_v.perturbed.y(:, 1), type="log", embedded=export_settings);
                    end
                end
                
                if isfield(figure{1}.include, "interferometer_nulling_statistics")
                    log_mean = mean(data.simulation.interferogram_multi.nulling, 2);
                    log_median = median(data.simulation.interferogram_multi.nulling, 2);
                    log_max = max(data.simulation.interferogram_multi.nulling, [], 2);
                    fraction_good = sum(data.simulation.interferogram_multi.nulling < 1e-4, 2) / size(data.simulation.interferogram_multi.nulling, 2);

                    export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_nulling_statistics.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_nulling_statistics.height), ...
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


                if isfield(figure{1}.include, "interferometer_response_function")
                    RFp = [data.simulation.interferogram_multi.RFn; data.simulation.interferogram.RFp];
                    export_settings = struct("width", export_data.sizes.width.(figure{1}.include.interferometer_response_function.width), ...
                                "height", export_data.sizes.height.(figure{1}.include.interferometer_response_function.height), ...
                                "font_size", export_data.sizes.font_size, ...
                                "name", label_name + sprintf("_rf_%d", j));
                    plot_response_function_theta(theta, RFp, "Normalize", true, "embedded", export_settings);

                end

        end
    end

    close all;

end


for table = export_data.tables

    switch table{1}.type
            case "unique_baselines"
                filename = path + "tables/" + export_data.chapters{table{1}.chapter}.name + "/" + table{1}.type;
                export_baselines_table(datas(2).instrument.unique_baselines, ...
                    datas(1).instrument.unique_baselines, filename)
    end

end
