clc; clear; close all;

data = readtable("others/MIR.csv", "VariableNamingRule", "preserve");

names = table2array(data(:, 1));
diameters = table2array(data(:, 2));
wavelengths = table2array(data(:, 3));

% Convert to arcsec
Resolutions = 1.22 * wavelengths ./ diameters * 206265 * 1000;

% figure;
% 
% for i = 1:numel(diameters)
%     semilogy(diameters(i), Resolutions(i), ".", "MarkerSize", 10)
%     hold on;
% end
% yline(0.1, "--", "LineWidth", 1.5)
% 
% legend(names, "Location", "Best")

final = table(names, diameters, wavelengths, Resolutions);