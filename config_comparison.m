%% CONFIG_COMPARISON Compare different configurations
% definition in terms of size and efficiencies. 
% Positions (normalised by max baseline)
% Modulation efficiency
% IWA, FWHM (normalised by lambda/max baselien)

clc; clear; close;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked')

% Config files to consider in the analysis. CHANGE HERE by selecting the
% type of the analysis to perform:
%   (1) Lay's configurations (Source: Lay 2005)
%   (2) Real past missions (Source: see thesis)
CONFIG_TYPE = 2;

if CONFIG_TYPE == 1
    configurations = ["lay_x_2", "lay_x_3", "lay_linear_A", "lay_linear_B", "lay_diamond"]';
    names = ["X-Array 2:1", "X-Array 3:1", "Linear DCB A", "Linear DCB B", "Diamond DCB"];
    caption = "Comparison table for five example configurations as presented in Lay \\cite{lay_imaging_2005} with normalised results for positions, modulation efficiency (higher is better), inner working angles (lower is better) and resolution angles (lower is better). Ap. is the aperture, $\varphi$ is the assigned phase and G is the nulling ratio (higher exponent is better).";
    label = "tab:modelling:comparison";
    filename = '/Users/francesco/Library/CloudStorage/OneDrive-Personale/Università/SEMESTRE 12/TESI/Thesis/main/tables/modelling/comparison';
elseif CONFIG_TYPE == 2
    configurations = ["darwin", "tpf_i", "life"]';
    names = ["Darwin", "TPF-I", "LIFE"];
    caption = "Comparison table for Darwin, TPF-I and LIFE with normalised results for positions, modulation efficiency (higher is better), inner working angles (lower is better) and resolution angles (lower is better). Ap. is the aperture, $\varphi$ is the assigned phase and G is the nulling ratio (higher exponent is better).";
    label = "tab:validation:comparison";
    filename = 'exports/images/table';
end

% Empty fields to fill
Normalised_positions_x = [];
Normalised_positions_y = [];
Amplitudes = [];
Diameter = [];
Phases = [];
Modulation_efficiency = [];
IWA = [];
FWHM = [];
Ratio = [];
B_opt = [];
AR_opt = [];

% For every configuration...
for i = 1:length(configurations)
    data = ReadYaml(sprintf('config/%s.yml', configurations(i)));
    data = convert_data(data);

    Normalised_positions_x(i, :) = data.instrument.positions(:, 1)' / data.instrument.baseline;
    Normalised_positions_y(i, :) = data.instrument.positions(:, 2)' / data.instrument.baseline;

    Amplitudes(i, :) = data.instrument.intensities;
    Diameter(i, :) = data.instrument.diameter;
    Phases(i, :) = rad2deg(data.instrument.phase_shifts);

    Ratio(i, 1) = ...
    compute_nulling_ratio(data.instrument.apertures, data.instrument.intensities, ...
    data.instrument.phase_shifts, data.instrument.positions, ...
    data.environment.stellar_angular_radius, data.instrument.wavelength);

    [PSF, FWHM(i, 1), C_i] = compute_psf(data.instrument.wavelength, ...
        data.instrument.intensities, data.instrument.positions, ...
        data.instrument.phase_shifts, data.environment.exoplanet_position, data.simulation.theta_range);
    FWHM(i, 1) = FWHM(i, 1) / (data.instrument.wavelength/data.instrument.baseline);
    
    Modulation_efficiency(i, 1) = compute_modulation_efficiency(...
        data.instrument.efficiencies.optical_line, ...
        data.instrument.surfaces, C_i);

    IWA(i, 1) = compute_IWA(data.instrument.unique_baselines, data.instrument.wavelength) / (data.instrument.wavelength/data.instrument.baseline);
    OWA = data.instrument.wavelength / data.instrument.diameter(1);
    
    [B_opt(i, 1), AR_opt(i, 1), NR_opt(i, 1:2), IWA_opt(i, 1:2)] = ...
        plot_array_size_sensitivity(data.instrument.intensities, data.instrument.array, ...
        data.instrument.baseline, data.instrument.apertures_ratio, data.instrument.phase_shifts, ...
        data.environment.stellar_angular_radius, data.instrument.wavelength, [], 0.4);

    [B_opt(i, 2), AR_opt(i, 2), NR_opt(i, 3:4), IWA_opt(i, 3:4)] = ...
        plot_array_size_sensitivity(data.instrument.intensities, data.instrument.array, ...
        data.instrument.baseline, data.instrument.apertures_ratio, data.instrument.phase_shifts, ...
        data.environment.stellar_angular_radius, data.instrument.wavelength, [], 0.5);

    [B_opt(i, 3), AR_opt(i, 3), NR_opt(i, 5:6), IWA_opt(i, 5:6)] = ...
        plot_array_size_sensitivity(data.instrument.intensities, data.instrument.array, ...
        data.instrument.baseline, data.instrument.apertures_ratio, data.instrument.phase_shifts, ...
        data.environment.stellar_angular_radius, data.instrument.wavelength, [], 0.6);

    for j = 1:6
        
        exotable = get_ppop_yield(IWA_opt(i, j), OWA, NR_opt(i, j), data.instrument.wavelength, "verbose", false, "create_plots", false);
        sumtable_ideal = extract_statistics_exoplanets(exotable, "Table");

        detections(i, j) = sumtable_ideal.MeanDetections;
    
    end
    
    % Remove to see sensitivity for each elements!
    close all;

end

% Convert into a table
results = table(Normalised_positions_x, Normalised_positions_y, ...
    Amplitudes, Diameter, Phases, Modulation_efficiency, IWA, FWHM, Ratio, B_opt, AR_opt, NR_opt, IWA_opt, detections, RowNames=configurations);

% save exports/results.mat results;

% Export table on the latex format.
export_comparison_table(results, names, filename, caption, label)
export_sensitivity_table(results, names, 'exports/images/sensitable')


% For the last array, compute some plot about how the dimensions affect
% some metrics.
plot_array_size_sensitivity(data.instrument.intensities, data.instrument.array, ...
    data.instrument.baseline, data.instrument.apertures_ratio, data.instrument.phase_shifts, ...
    data.environment.stellar_angular_radius, data.instrument.wavelength)