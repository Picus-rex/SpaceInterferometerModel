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