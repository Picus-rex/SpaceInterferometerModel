function [T, total_yield_matrix] = get_ppop_yield(IWAs, OWA, ratios, lambda, sim_options, export_setup)
%GET_PPOP_YIELD Verify the yield of the array from comparing the nulling
%ratio of the array over the integration time of the system when subjected
%to perturbations.
%
% INPUTS:
%   IWAs[Ns x 1]        Inner working angle for each simulation [rad]
%   OWA[Ns x 1]         Outer working angle [rad]
%   ratios[Np x Ns]     Nulling ratio for each point of the decentered
%                       ray tracing and for each simulation [-]
%   lambda[1]           Wavelength [m]
% 
% ARGUMENT OPTIONAL INPUTS:
%   create_plots        Flag for plot generations
%   integration_time    Integration time for incoming spectral radiance [s]
%   population          Either "PPOP" or "NASA"
%   universes_to_plot   How many plots to randomly generate from available.
%   verbose             If true, display text.
%
% OUTPUTS:
%   T[table]            Exoplanets table from P-Pop with additional fields,
%                       in particular:
%       yields[Ne x Ns] Number of detected exoplanet for each simulation
%                       and for each universe of the simulation.
%       candidates[Ne x Ns] Exoplanets that could be detected if the
%                       integration time is infinite. 
%   total_yield_matrix[Nu x Ns] Total exoplanet detection for each universe
%                       and simulation. 
%
% NOTE:
%   Assumed 0.2 as AgeomVIS from C. Dandumont when using NASA table.
%
% VERSION HISTORY:
%   2025-04-07 -------- 1.0
%   2025-04-08 -------- 1.1
%                     - Added debug printing information and slightly
%                       changed the way rms ratios are addressed for plots.
%   2025-04-10 -------- 1.2
%                     - Added option for universes_to_plot in
%                       plot_ppop_yield.
%                     - Added total_yield_matrix as export.
%                     - Added compatibility with real exoplanets.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    IWAs (:, :)
    OWA (1, 1)
    ratios (:, :)
    lambda (1, 1)
    sim_options.integration_time = 1 * 24 * 3600;
    sim_options.verbose = false;
    sim_options.population {mustBeMember(sim_options.population,["PPOP", "NASA"])}  = "PPOP"
    export_setup.create_plots = true;
    export_setup.universes_to_plot = 5;
end

% Extract the number of simulations from the system.
Ns = size(IWAs, 1);

% Constants
AU = 1.496e11;                  % Astronomical unit in meters
h = 6.626e-34;                  % Planck constant [J s]
c = 3e8;                        % Speed of light [m/s]
k = 1.381e-23;                  % Boltzmann constant [J/K]
R_earth = 6.371e6;              % Earth radius [m]
R_sun   = 6.957e8;              % Solar radius [m]

if strcmp(sim_options.population, "PPOP")
    % Extract the exoplanets from the table
    opts = detectImportOptions('others/TestPlanetPopulation.txt', 'NumHeaderLines', 1);
    T = readtable('others/TestPlanetPopulation.txt', opts);
    T.ang_sep_rad = arcsec2rad(T.AngSep);
else
    % Fill needed rows
    T = load_NASA_archive();
    T.Nuniverse(:) = 0;
    T.ang_sep_rad = T.AngSep;
    T.AgeomVIS(:) = 0.2;
    export_setup.universes_to_plot = 1;
end

% Extract the number of universes from the table
universes = unique(T.Nuniverse);
Nu = length(universes);

% Spectral radiance
B = @(T_val) (2*h*c^2 ./ lambda^5) ./ (exp(h*c./(lambda*k*T_val)) - 1);

% Compute for all the viable planets the ratio under the given integration time. 
% Thermal emission from the exoplanet and star
thermal_flux_planet = B(T.Tp) .* (pi * (T.Rp * R_earth).^2) * sim_options.integration_time;
thermal_flux_star = B(T.Ts) .* (pi * (T.Rs * R_sun).^2) * sim_options.integration_time;

% Reflected light from the exoplanet
reflected_flux_planet = thermal_flux_star .* T.AgeomVIS .* ((T.Rp * R_earth).^2 ./ (T.ap * AU).^2) * sim_options.integration_time;

% Contrasts
thermal_contrast = thermal_flux_planet ./ thermal_flux_star;
reflected_contrast = reflected_flux_planet ./ thermal_flux_star;
T.contrast = thermal_contrast + reflected_contrast;

% Allocate space
T.yields = zeros(size(T, 1), Ns);
T.candidates = zeros(size(T, 1), Ns);
rms_ratios = zeros(Ns, 1);

% For every simulation...
for i = 1:Ns
    IWA = IWAs(i);
    rms_ratios(i) = rms(ratios(:, i));

    T.candidates(:, i) = T.ang_sep_rad(:) > IWA & T.ang_sep_rad(:) < OWA;
    T.yields(:, i) = T.ang_sep_rad(:) > IWA & T.ang_sep_rad(:) < OWA & T.contrast(:) > rms_ratios(i);
    
    % For debug
    if sim_options.verbose
        fprintf("%d\t%.3e\t%d\n", i, rms_ratios(i), sum(T.yields(:, i)))
    end
end

% Allocate space for additional results
total_yield_matrix = zeros(Nu, Ns); % Nu x Ns

for u = 1:Nu
    idx = T.Nuniverse == u-1;                           % Logical indexing for current universe
    total_yield_matrix(u, :) = sum(T.yields(idx, :));   % Total across planets in the universe    
end

% For debug
if sim_options.verbose
    fprintf("\nAverage for universe in the best case:\t%.2f planets\n", max(mean(total_yield_matrix, 1)));
    fprintf("Average for universe in the mean case:\t%.2f planets\n", mean(mean(total_yield_matrix, 2)));
    fprintf("Average for all universes:\t%.2f planets\n", mean(sum(total_yield_matrix, 1)));
end


% PLOT SECTION
if export_setup.create_plots
    plot_ppop_yield(T, mean(IWAs), OWA, min(rms_ratios), export_setup);
end

end