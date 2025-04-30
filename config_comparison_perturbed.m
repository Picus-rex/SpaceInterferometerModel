%% CONFIG_COMPARISON_PERTURBE Compare different configurations following 
% Lay definition in terms of size and efficiencies. 
% Positions (normalised by max baseline)
% Modulation efficiency
% IWA, FWHM (normalised by lambda/max baselien)

clc; clear; close;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked')

% Config files to consider in the analysis.
configurations = ["lay_x_2", "lay_x_3", "lay_linear_A", "lay_linear_B", "lay_diamond"]';
configurations = ["lay_linear_A"]';
names = ["X-Array 2:1", "X-Array 3:1", "Linear DCB A", "Linear DCB B", "Diamond DCB"];
filename = '/Users/francesco/Library/CloudStorage/OneDrive-Personale/Università/SEMESTRE 12/TESI/Thesis/67a4aa2876545e0e504145e0/tables/modelling/comparison';


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

% For every configuration...
for i = 1:length(configurations)
    fprintf("Config %d\n", i);
    data = convert_data(sprintf('config/%s.yml', configurations(i)));

    % Extract disturbances on phases
    dist_phases = data.simulation.code_v.perturbed.phase;
    perturbed_vect = zeros(size(data.instrument.phase_shifts));

    Normalised_positions_x(i, :) = data.instrument.positions(:, 1)' / data.instrument.baseline;
    Normalised_positions_y(i, :) = data.instrument.positions(:, 2)' / data.instrument.baseline;

    Amplitudes(i, :) = data.instrument.intensities;
    Diameter(i, :) = data.instrument.diameter;
    Phases(i, :) = rad2deg(data.instrument.phase_shifts);

    Ns = size(dist_phases, 2);

    for j = 1:Ns
        perturbed_vect(1) = rms(dist_phases(:, j));
        phases = data.instrument.phase_shifts + perturbed_vect;

        Ratio(i, j) = ...
        compute_nulling_ratio(data.instrument.apertures, ...
        data.instrument.intensities, phases, data.instrument.positions, ...
        data.environment.stellar_angular_radius, data.instrument.wavelength);

        [PSF, FWHM(i, j), C_i] = compute_psf(data.instrument.wavelength, ...
        data.instrument.intensities, data.instrument.positions, ...
            phases, data.environment.exoplanet_position, data.simulation.theta_range);
        
        FWHM(i, j) = FWHM(i, j) / (data.instrument.wavelength/data.instrument.baseline);
    
        Modulation_efficiency(i, j) = compute_modulation_efficiency(...
            data.instrument.efficiencies.optical_line, ...
            data.instrument.surfaces, C_i);
    end

    IWA(i, 1) = compute_IWA(data.instrument.unique_baselines, data.instrument.wavelength) / (data.instrument.wavelength/data.instrument.baseline);
end

% Convert into a table
results = table(Normalised_positions_x, Normalised_positions_y, ...
    Amplitudes, Diameter, Phases, Modulation_efficiency, IWA, FWHM, Ratio, RowNames=configurations);

% save exports/pert_results.mat results;

% Export table on the latex format. MUST BE CORRECTED (show the average?)
% export_comparison_table(results, names, filename)

%%

nominal = load("exports/results.mat", "results");
perturbed = load("exports/pert_results.mat", "results");
names = ["X-Array 2:1", "X-Array 3:1", "Linear DCB A", "Linear DCB B", "Diamond DCB"];

compare_configurations(nominal.results, perturbed.results, names)