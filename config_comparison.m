%% CONFIG_COMPARISON Compare different configurations following Lay 
% definition in terms of size and efficiencies. 
% Positions (normalised by max baseline)
% Modulation efficiency
% IWA, FWHM (normalised by lambda/max baselien)

clc; clear; close;

configurations = ["lay_x_2", "lay_x_3", "lay_linear_A", "lay_linear_B", "lay_diamond"]';
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') % Change to NORMAL to export

Normalised_positions_x = [];
Normalised_positions_y = [];
Amplitudes = [];
Diameter = [];
Phases = [];
Modulation_efficiency = [];
IWA = [];
FWHM = [];

for i = 1:length(configurations)
    data = ReadYaml(sprintf('config/%s.yml', configurations(i)));
    data = convert_data(data);

    Normalised_positions_x(i, :) = data.instrument.positions(:, 1)' / data.instrument.baseline;
    Normalised_positions_y(i, :) = data.instrument.positions(:, 2)' / data.instrument.baseline;

    Amplitudes(i, :) = data.instrument.intensities;
    Diameter(i, :) = data.instrument.diameter;
    Phases(i, :) = rad2deg(data.instrument.phase_shifts);

    [PSF, FWHM(i, 1), C_i] = compute_psf(data.instrument.wavelength, ...
        data.instrument.intensities, data.instrument.positions, ...
        data.instrument.phase_shifts, data.environment.exoplanet_position, data.simulation.theta_range);
    FWHM(i, 1) = FWHM(i, 1) / (data.instrument.wavelength/data.instrument.baseline);

    [Modulation_efficiency(i, 1), IWA(i, 1)] = compute_IWA(PSF, ...
        data.simulation.theta_range, ...
        data.instrument.efficiencies.optical_line, ...
        data.instrument.surfaces, C_i, false);
    IWA(i, 1) = IWA(i, 1) / (data.instrument.wavelength/data.instrument.baseline);
end

results = table(Normalised_positions_x, Normalised_positions_y, ...
    Amplitudes, Diameter, Phases, Modulation_efficiency, IWA, FWHM, RowNames=configurations);