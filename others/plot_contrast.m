% Clear workspace and figures
clear; close all; clc;

% Specify the file name
filename = 'others/TestPlanetPopulation.txt';

% Read the catalogue.
% (Assumes whitespace‐delimited text file with header rows.)
% If your file has two header lines, you may need to skip one of them.
% Here we use the first header line (column names)
opts = detectImportOptions(filename, 'Delimiter', {'\t',' '});
% Sometimes the file might have multiple header lines.
% Adjust opts.DataLines if needed (e.g., opts.DataLines = [3, Inf];)
T = readtable(filename, opts);

% Look at the variable names:
disp('Variables in the table:');
disp(T.Properties.VariableNames);

% Filter to select only those planets from desired universe.
% According to your file, the column corresponding to the universe number is named either 'nMC' or 'Nuniverse'.
% Here we assume the first column (nMC) is used.
universeToPlot = 0;
rows = T.nMC == universeToPlot;
Tsub = T(rows, :);

% For the x axis, we use the stellar type.
% In your file the column is named 'stype'.
% Convert to categorical for proper plotting.
stellarTypes = categorical(Tsub.stype);

% Flux contrast
h = 6.626e-34;                  % Planck constant [J s]
c = 3e8;                        % Speed of light [m/s]
k = 1.381e-23;                  % Boltzmann constant [J/K]
R_earth = 6.371e6;              % Earth radius [m]
R_sun   = 6.957e8;              % Solar radius [m]
Baseline = 10;                  % Baseline [m]
Diameter = 2.5;                 % Diameter [m]

lambda1 = 10e-6;                 % Wavelength [m] (10 microns)
lambda2 = 500e-9;                 % Wavelength [m] (10 microns)

B = @(T_val, lam) (2*h*c^2 ./ lam^5) ./ (exp(h*c./(lam*k*T_val)) - 1);

contrast_1 = ((Tsub.rp * R_earth).^2 .* B(Tsub.Tp, lambda1)) ./ ((Tsub.Rs * R_sun).^2 .* B(Tsub.Ts, lambda1));
contrast_2 = ((Tsub.rp * R_earth).^2 .* B(Tsub.Tp, lambda2)) ./ ((Tsub.Rs * R_sun).^2 .* B(Tsub.Ts, lambda2));

set(0, 'DefaultFigureWindowStyle', 'normal') % Change to NORMAL to export
col = styling;

figure;
semilogy(stellarTypes, contrast_1, 'o', 'LineWidth', 1.5, 'MarkerSize', 8, "Color", col(1, :));
xlabel('Stellar Spectral Type');
ylabel('Planet/Star Flux Contrast');
grid on;
colorbar('off')
styling(true, 10, 7, "exports/contrast_1", false);

figure;
semilogy(stellarTypes, contrast_2, 'o', 'LineWidth', 1.5, 'MarkerSize', 8, "Color", col(3, :));
xlabel('Stellar Spectral Type');
ylabel('Planet/Star Flux Contrast');
ylim([1e-20, 1e-5])
grid on;
colorbar('off')
styling(true, 10, 7, "exports/contrast_2", false);