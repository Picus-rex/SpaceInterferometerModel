%% EXOPLANETS YIELDS
% Analyse how a system is detecting exoplanets. 

clc; clear; close;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked')

% Loading of elements from configurations files
data = convert_data('config/darwin.yml');

% Compute IWA, OWA
IWA = compute_IWA(data.instrument.unique_baselines, data.instrument.wavelength);
OWA = data.instrument.wavelength / data.instrument.diameter(1);

%% Ideal case

[ratio, ~] = compute_nulling_ratio(data.instrument.apertures, data.instrument.intensities, ...
     data.instrument.phase_shifts, data.instrument.positions, ...
     data.environment.stellar_angular_radius, data.instrument.wavelength);

exotable = get_ppop_yield(IWA, OWA, ratio, data.instrument.wavelength, "verbose", true);

sumtable_ideal = extract_statistics_exoplanets(exotable, "Darwin ideal");

startable = extract_stars(exotable, "compute_HZ", true);
summaryTable = summarise_HZ_by_type(startable);

plot_HZ_limits(startable)
plot_HZ_limits([], summaryTable)

plot_ppop_yield(exotable, IWA, OWA, ratio, "HZ_table", summaryTable, "plot_heatmap", false)
plot_ppop_yield(exotable, IWA, OWA, ratio, "HZ_table", summaryTable, "plot_heatmap", false, "filter_types", ["F", "K"], "show_legend", "hide");
plot_ppop_yield(exotable, IWA, OWA, ratio, "HZ_table", summaryTable, "plot_heatmap", false, "filter_types", ["G", "M"], "show_legend", "hide");

%% Perturbed cases 

Ns = data.simulation.code_v.num_perturbed_simulations;
phases = data.simulation.code_v.corrected.phase;

% Extract nulling ratio
[ratios, ~, ~]  = perturbate_system(data.instrument.apertures, ...
    data.instrument.intensities, data.instrument.phase_shifts, ...
    data.instrument.positions, data.environment.stellar_angular_radius, ...
    data.instrument.wavelength, phases, data.instrument.combination, ...
    perturbed_map_plotting_number=1, compute_PSF=false, create_plots=false);

% Get yield from PPOP
exotable = get_ppop_yield(IWA, OWA, ratios, data.instrument.wavelength, "verbose", true);

sumtable_pert = extract_statistics_exoplanets(exotable, "Darwin perturbed");

%% Config analysis

configurations = ["darwin", "tpf_i", "life"]';
names = ["Darwin", "TPF-I", "LIFE"];

sumtable = table();

% For every configuration...
for i = 1:length(configurations)
    data = convert_data(sprintf('config/%s.yml', configurations(i)));
    
    IWA = compute_IWA(data.instrument.unique_baselines, data.instrument.wavelength);
    OWA = data.instrument.wavelength / data.instrument.diameter(1);

    [ratio, ~] = compute_nulling_ratio(data.instrument.apertures, data.instrument.intensities, ...
     data.instrument.phase_shifts, data.instrument.positions, ...
     data.environment.stellar_angular_radius, data.instrument.wavelength);

    exotable = get_ppop_yield(IWA, OWA, ratio, data.instrument.wavelength, "verbose", false, "create_plots", false);

    sumtable(i, :) = extract_statistics_exoplanets(exotable, names(i));

end

%% Export table

filename = '/Users/francesco/Library/CloudStorage/OneDrive-Personale/Università/SEMESTRE 12/TESI/Thesis/main/tables/validation/exo_data';
export_exoplanets_table(sumtable, filename);