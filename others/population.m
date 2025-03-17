clc; clear; close;

% IWA, OWA from Viseurs

% Read the catalogue
opts = detectImportOptions('others/TestPlanetPopulation.txt', 'NumHeaderLines', 1);
T = readtable('others/TestPlanetPopulation.txt', opts);

universeToPlot = 0;
rows = T.Nuniverse == universeToPlot;
T = T(rows, :);


arcsec2rad = (pi/180)/3600;
ang_sep_rad = T.AngSep * arcsec2rad;

% Flux contrast
h = 6.626e-34;                  % Planck constant [J s]
c = 3e8;                        % Speed of light [m/s]
k = 1.381e-23;                  % Boltzmann constant [J/K]
lambda = 10e-6;                 % Wavelength [m] (10 microns)
R_earth = 6.371e6;              % Earth radius [m]
R_sun   = 6.957e8;              % Solar radius [m]
Baseline = 10;                  % Baseline [m]
Diameter = 2.5;                 % Diameter [m]

resol = 0.22 * lambda/Baseline;
IWA = lambda/Baseline;
OWA = lambda/Diameter;

B = @(T_val) (2*h*c^2 ./ lambda^5) ./ (exp(h*c./(lambda*k*T_val)) - 1);

contrast = ((T.rp * R_earth).^2 .* B(T.Tp)) ./ ((T.Rs * R_sun).^2 .* B(T.Ts));

%% Plot

figure;
hold on; 
set(gca, 'XScale', 'log', 'YScale', 'log'); 

starTypes = unique(T.Stype);

% Define markers and colors for different groups (adjust as desired)
markers = {'o','s','^','d','v','>','<','p','h'};
%colors = lines(numel(starTypes));

style_colors;

for i = 1:length(starTypes)

    % Select indices for the current star type
    idx = strcmp(T.Stype(:), {starTypes{i}});

    % Create a log–log scatter plot for this group
    scatter(ang_sep_rad(idx) * 648000000/pi, contrast(idx), 40, ...
        colours(i,:), 'filled', 'DisplayName', starTypes{i});

end

% xline(resol, '--', "LineWidth", 1.5, "DisplayName", "Angular resolution");
% xline(IWA, 'b--', "LineWidth", 1.5, "DisplayName", "IWA");
% xline(OWA, 'r--', "LineWidth", 1.5, "DisplayName", "OWA");

xlabel('Angular Separation [mas]');
ylabel('Flux Contrast at 10 \mum');
% title('Exoplanet Population: Angular Separation vs. Flux Contrast');
legend('show','Location','best');
grid on;
hold off;

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);


set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 18, 8]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [18, 8]);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, "exports/ppop_population", '-dpdf', '-r300');