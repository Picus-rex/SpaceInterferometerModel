%% Point Spread Function of a nulling interferometry
% Based on Lay OP. Imaging properties of rotating nulling interferometers. 
% 2005, this script visualises the PSF associated to a generic exoplanet in
% a given position.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') % Change to NORMAL to export

%% Data
% Taken from the paper

data = ReadYaml('config/psf_array.yml');

lambda = data.instrument.wavelength;
B = data.instrument.baseline;
A = cell2mat(data.instrument.intensities);
phase_shifts = cell2mat(data.instrument.phase_shifts);
theta_range = define_range(data.simulation.angular_extension);

positions = define_array(data.instrument.array, B, data.instrument.apertures_ratio);

conversion_rad2mas = 1e3 * (3600 * 180) / pi;

%% Compute the PSF
PSF = compute_psf(lambda, A, positions, phase_shifts, [1.453e-7, 0], theta_range);

% The maximum of the PSF is a first index of the integration time
% measurement that can be taken as a reference.
fprintf("PSF Max Value: \t%.2f\n", max(PSF, [], "all"));


% Compute the angular resolution to plot it on the graph
r = 1.22 * lambda/B;
angles = linspace(0, 2*pi, 100);
x = r * cos(angles);
y = r * sin(angles);

figure; hold on;
imagesc(theta_range * conversion_rad2mas, theta_range * conversion_rad2mas, PSF);
plot(x * conversion_rad2mas, y * conversion_rad2mas, 'r--')
colorbar;
title("PSF");

%% Compute the IWA

[D, a] = plot_apertures(positions, A, ...
    data.instrument.efficiencies.optical_line, ...
    data.instrument.efficiencies.beam_combiner, false);

eta_opt = data.instrument.efficiencies.optical_line;
[eta_mod, IWA] = compute_IWA(PSF, theta_range, eta_opt, a, true);

fprintf("Inner working angle: \t%.10f\n", IWA);