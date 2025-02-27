%% Display PSF
% Based on 
% 1. Lay OP. Imaging properties of rotating nulling interferometers. 2005; 

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') % Change to NORMAL to export

%% Data
% Taken from the paper

lambda = 10e-6;
B = 89;
positions = [-40, 20;
             -40, -20;
             40, 20;
             40, -20];
% B = 10;
% positions = [-5,    0;
%              -2.5, -0;
%               2.5,  0;
%               5,   -0];

A = [0.5, 0.5, 0.5, 0.5];

phase_shifts = [0, pi, pi/2, 3/2*pi];
%phase_shifts = [0, 3/2*pi, pi/2, pi];

theta_range = linspace(-2.90888208665724e-7, 2.90888208665724e-7, 2000);    % Angular grid [rad]

%% Compute the PSF
PSF = compute_psf(lambda, A, positions, phase_shifts, [1.453e-7, 0], theta_range);

% The maximum of the PSF is a first index of the integration time
% measurement that can be taken as a reference.
disp(max(PSF, [], "all"))

% Compute the angular resolution to plot it on the graph
r = 1.22 * lambda/B;
angles = linspace(0, 2*pi, 100);
x = r * cos(angles);
y = r * sin(angles);

figure; hold on;
imagesc(theta_range, theta_range, PSF);
plot(x, y, 'r--')
colorbar;
title("PSF");

%% Compute the IWA

[D, a] = plot_apertures(positions, A, false);

eta_opt = 0.85 * ones(4, 1);
[eta_mod, IWA] = compute_IWA(PSF, theta_range, eta_opt, a, true);

disp(IWA)