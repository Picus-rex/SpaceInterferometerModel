function T = get_ppop_yield(IWAs, OWA, ratios, lambda, sim_options, export_setup)
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
% ARGUMENT OPTINAL INPUTS:
%   integration_time    Integration time for incoming spectral radiance [s]
%   create_plots        Flag for plot generations
%
% OUTPUTS:
%   T[table]            Exoplanets table from P-Pop with additional fields,
%                       in particular:
%       yields[Ne x Ns] Number of detected exoplanet for each simulation
%                       and for each universe of the simulation.
%       candidates[Ne x Ns] Exoplanets that could be detected if the
%                       integration time is infinite. 
%
% VERSION HISTORY:
%   2025-04-07 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    IWAs (:, :)
    OWA (1, 1)
    ratios (:, :)
    lambda (1, 1)
    sim_options.integration_time = 1 * 24 * 3600;
    export_setup.create_plots = true;
end

% Extract the number of simulations from the system.
Ns = size(IWAs, 1);

% Extract the exoplanets from the table
opts = detectImportOptions('others/TestPlanetPopulation.txt', 'NumHeaderLines', 1);
T = readtable('others/TestPlanetPopulation.txt', opts);

% Constants
AU = 1.496e11;                  % Astronomical unit in meters
h = 6.626e-34;                  % Planck constant [J s]
c = 3e8;                        % Speed of light [m/s]
k = 1.381e-23;                  % Boltzmann constant [J/K]
R_earth = 6.371e6;              % Earth radius [m]
R_sun   = 6.957e8;              % Solar radius [m]
arcsec2rad = (pi/180)/3600;     % arcsec to radians conversion factor

% Spectral radiance
B = @(T_val) (2*h*c^2 ./ lambda^5) ./ (exp(h*c./(lambda*k*T_val)) - 1);

% Compute for all the viable planets the ratio under the given integration time. 
T.ang_sep_rad = T.AngSep * arcsec2rad;

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

% For every simulation...
for i = 1:Ns
    IWA = IWAs(i);
    ratio = rms(ratios(:, i));

    T.candidates(:, i) = T.ang_sep_rad(:) > IWA & T.ang_sep_rad(:) < OWA;
    T.yields(:, i) = T.ang_sep_rad(:) > IWA & T.ang_sep_rad(:) < OWA & T.contrast(:) > ratio;
end

% PLOT SECTION
if export_setup.create_plots
    plot_ppop_yield(T, mean(IWAs), OWA, export_setup);
end

end