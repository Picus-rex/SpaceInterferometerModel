function [ratios, Maps, Nominal_Map, PSFs, Nominal_PSF] = ...
    perturbate_system(N, amplitudes_nominal, ...
    phases_nominal, positions_nominal, theta_star, lambda, ...
    phases_perturbed, combinations, theta_x, theta_y, ...
    exoplanet_position, export_setup)
%PERTURBATE_SYSTEM From the nominal characteristics of the system, computes
%the variations based on the provided phases shift. 
%
% INPUTS:
%   N[1]                Number of apertures [-]
%   amplitude_nom[Nx1]  Vector of nominal aperture amplitudes [m^2]
%   phases_nom[Nx1]     Vector of nominal phase shifts for each aperture
%   positions_nom[Nx2]  Matrix of (x, y) nominal positions of the apertures
%   theta_star[1]       Angular dimension of the star [rad]
%   lambda[1]           Wavelength [m]
%   phases_pert[NxNs]   Resulting perturbed phases for one aperture over
%                       all the rays (N points) over Ns simulations
%   combinations[Nx1]   Vector of beam combiner coefficients
%
% OPTIONAL INPUTS: 
%   theta_x, theta_y    Meshgrid of angular coordinates
%   exoplanet_position  Position of the exoplanet for PSF computation
% 
% ARGUMENT OPTINAL INPUTS:
%   create_plots[bool]  If true, creates several plots.
%   type[string]        String following "%unit_%de(+-)%d",
%                       where the exponential number is reversed and
%                       recognised as a scale (if m_1e-6 is given, for
%                       example, the string is converted to micro m). 
%   perturbed_map_plotting_number[1] Number of perturbed maps to plot (must
%                       be below Ns.
%
% OUTPUTS:
%   ratios[N x Ns]      Nulling ratios associated to every provided OPD.
%   Maps[TX x TY x Ns]  3D matrix where the first two dimensions are the
%                       standard normalised response function along the
%                       angular coordinates provided.
%   Nominal_Map[TX x TY]Nominal response map.
%   PSFs[TX x TY x Ns]  3D matrix where the first two dimensions are the
%                       PSFs along the angular coordinates provided.
%   Nominal_PSF[TX x TY]Nominal PSF.
%
% VERSION HISTORY:
%   2025-04-02 -------- 1.0
%   2025-04-03 -------- 1.1
%                     - Added support for PSF with integration of inputs
%                       and outputs.
%                     - Moved the plots of variation maps to a dedicated
%                       function.
%   2025-04-24 -------- 1.2
%                     - Added additional exports to speed up the analysis.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    N {mustBeInteger}
    amplitudes_nominal (:, :) {mustBeNumeric}
    phases_nominal (:, :) {mustBeNumeric}
    positions_nominal (:, :) {mustBeNumeric}
    theta_star (1, 1) {mustBeNumeric}
    lambda (1, 1) {mustBeNumeric}
    phases_perturbed (:, :) {mustBeNumeric}
    combinations (:, :) {mustBeNumeric}
    theta_x = load_grid("x")
    theta_y = load_grid("y")
    exoplanet_position = [1.453e-7, 0]
    export_setup.create_plots = true
    export_setup.type = "m_1e-6"
    export_setup.perturbed_map_plotting_number = 3
    export_setup.compute_PSF = true
    export_setup.embedded = []
end

% Number of series, obtain OPDs as well as the nominal map
Ns = size(phases_perturbed, 2);
OPDs = phase2opd(phases_perturbed, lambda);
nom_table = compute_response_function("lambda", lambda, "N", N, ...
        "positions", positions_nominal, "A", amplitudes_nominal, ...
        "phase_shifts", phases_nominal, "combinations", combinations, ...
        "theta_x", theta_x, "theta_y", theta_y);
Nominal_Map = nom_table(:, :).T_standard;
Nominal_PSF = compute_psf(lambda, amplitudes_nominal, positions_nominal, ...
                    phases_nominal, exoplanet_position, theta_x(1, :));

% Allocation and reordering of the vector
phases_perturbed_all = reshape(phases_perturbed, [], 1);
ratios = zeros(length(phases_perturbed_all), 1);
perturbed_vect = zeros(size(phases_nominal));
Maps = zeros(length(theta_x), length(theta_y), Ns);
PSFs = zeros(length(theta_x), length(theta_y), Ns);

% Add the single perturbation to the phases to compute the nulling ratio
for i = 1:length(phases_perturbed_all)
    perturbed_vect(1) = phases_perturbed_all(i);
    phases = phases_nominal + perturbed_vect;

    amplitudes_modified = amplitudes_nominal .* combinations;

    ratios(i) = compute_nulling_ratio(N, amplitudes_modified, phases, ...
                                    positions_nominal, theta_star, lambda);
end

for i = 1:Ns
    perturbed_vect(1) = rms(phases_perturbed(:, i));
    phases = phases_nominal + perturbed_vect;

    map_table = compute_response_function("lambda", lambda, "N", N, ...
        "positions", positions_nominal, "A", amplitudes_nominal, ...
        "phase_shifts", phases, "combinations", combinations, "theta_x", ...
        theta_x, "theta_y", theta_y);

    Maps(:, :, i) = map_table(:, :).T_standard;
    
    if export_setup.compute_PSF
        PSFs(:, :, i) = compute_psf(lambda, amplitudes_nominal, ...
            positions_nominal, phases, exoplanet_position, theta_x(1, :));
    end

end

% Reformat the vector 
ratios = reshape(ratios, [], Ns);

if export_setup.create_plots

    % Get visualisation style
    [scale, scale_tag] = get_scale_plots(export_setup.type);
    elem_label = sprintf("OPD [%s]", scale_tag);
    
    % OPDs and nulling ratios
    figure; hold on;
    cols = get_colours(Ns);
    for i = 1:Ns
        scatter(OPDs(:, i) * scale, ratios(:, i), 36, cols(i, :), ...
                           "filled", "DisplayName", sprintf("Sim. %d", i));
    end
    xlabel(elem_label)
    ylabel("Nulling ratios")
    legend("Location","best");
    grid minor;
    set(gca, "YScale", "log");

    % Transmission maps
    plot_variation_map(Ns, theta_x, theta_y, Maps, Nominal_Map, export_setup);

    % PSF
    plot_variation_map(Ns, theta_x, theta_y, PSFs, Nominal_PSF, export_setup);

    % Monodirectional transmission maps
    maps_named = table();
    names_plot = strings(Ns, 1);
    for i = 1:Ns
        names_plot(i) = sprintf("Simulation %d", i);
        maps_named.(names_plot(i)) = Maps(:, :, i);
    end
    plot_transmission_map_monodirectional(theta_x(1, :), maps_named, names_plot, NaN, "tendencies", 2);

end



end