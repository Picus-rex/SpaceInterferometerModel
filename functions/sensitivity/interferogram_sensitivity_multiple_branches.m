function [RFs, RFn, R, ratio] = ...
    interferogram_sensitivity_multiple_branches(nominal_file, ...
    perturbed_file, intensities, nominal_phases, positions, combination,...
    lambda, theta, stellar_angular_radius, maps_to_compute, simulations)
%INTERFEROGRAM_SENSITIVITY_MULTIPLE_BRANCHESE Computes the nominal and
%perturbed response functions for interferogram sensitivity analysis, by
%applying a random perturbation to each aperture.
%
% ARGUMENT INPUTS:
%   nominal_file[string] Path to the nominal CSV file from CODE V.
%   perturbed_file[string] Path to the perturbed CSV from CODE V.
%   intensities[Nx1]    Vector of aperture intensities [V/m]
%   nominal_phases[Nx1] Vector of nom. phase shifts for each aperture [rad]
%   positions[Nx2]      Matrix of (x, y) positions of the apertures [m]
%   combination[Nx1]    Vector of beam combiner coefficients [-]
%   lambda[1]           Wavelength [m] (default: 1e-5 m)
%   theta[Nt x 1]       Vector of angular coordinates [rad] 
%                       (default: mas2rad(linspace(-700, 700, 1000)))
%   stellar_angular_radius[1] Angular radius of the star [rad] 
%                       (default: Proxima Centauri value)
%   maps_to_compute[1xNm] Indices of the maps to include in the R 3D matrix
%                       (default: [1, 251, 501, 751])
%   simulations[1]      Crossed simulations to perform in the random
%                       selection process.
%
% OUTPUTS:
%   RFs[Ns x Nt]        Matrix of perturbed response functions
%   RFn[1 x Nt]         Vector of nominal response function
%   R[Np x Nt x Ns]     3D matrix of the response map over apertures
%   ratio[Np x Ns]      Nulling ratio for every point for every simulation
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    nominal_file {mustBeFile}
    perturbed_file {mustBeFile}
    intensities 
    nominal_phases 
    positions 
    combination 
    lambda = 1e-5
    theta = mas2rad(linspace(-700, 700, 1000))
    stellar_angular_radius = 1.50444523662918e-09
    maps_to_compute = 1 : floor(length(theta) / 4) : length(theta)
    simulations = 100
end

% Load data from CODE V
[~, OPD, ~, ~] = load_opd(perturbed_file);     % size: Np x Ns
[~, OPD_n, ~, ~] = load_opd(nominal_file);     % size: Np x 1

% Setup
delta_phi_n = opd2phase(OPD_n(:, 1), lambda);  % Np x 1
Np = size(OPD, 1);
N = size(positions, 1);
Nt = length(theta);
Ns = size(OPD, 2);  % Total available precomputed perturbations

% Nominal response
phases_nominal = nominal_phases + delta_phi_n;
[~, RFn] = create_interferogram(positions, intensities, combination, phases_nominal, lambda, theta);

% Allocate space
RFs = zeros(simulations, Nt);
R = zeros(Np, Nt, length(maps_to_compute));
ratio = zeros(Np, simulations);

% Run randomized simulations
jR = 1;
for i = 1:simulations
    % Allocate space for this simulation's phase perturbations
    phase_perturbations = zeros(Np, N);

    % Randomly perturb each aperture independently
    for k = 1:N
        random_column = randi(Ns);  % pick a random simulation
        delta_phi_sample = opd2phase(OPD(:, random_column), lambda);
        phase_perturbations(:, k) = nominal_phases(k) + delta_phi_sample;
    end

    % Compute the response function for this perturbed configuration
    [Rm, RFs(i, :)] = create_interferogram(positions, intensities, combination, ...
        phase_perturbations, lambda, theta);

    if ismember(i, maps_to_compute)
        R(:, :, jR) = Rm;
        jR = jR + 1;
    end

    % Compute nulling ratio for this perturbation
    ratio(:, i) = compute_nulling_ratio(N, intensities, phase_perturbations, ...
        positions, stellar_angular_radius, lambda);
end

end