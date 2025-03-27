%% CONFIG_COMPARISON Compare different configurations following Lay 
% definition in terms of size and efficiencies. 
% Positions (normalised by max baseline)
% Modulation efficiency
% IWA, FWHM (normalised by lambda/max baselien)

clc; clear; close;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') % Change to NORMAL to export

% Config files to consider in the analysis.
configurations = ["lay_x_2", "lay_x_3", "lay_linear_A", "lay_linear_B", "lay_diamond"]';
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
end

% Convert into a table
results = table(Normalised_positions_x, Normalised_positions_y, ...
    Amplitudes, Diameter, Phases, Modulation_efficiency, IWA, FWHM, Ratio, RowNames=configurations);


export_comparison_table(results, names, filename)