function [RFs, RFn, R, ratio] = interferogram_sensitivity_multiple_branches(OPD_n, ...
    OPD, intensities, nominal_phases, positions, combination,...
    lambda, theta, maps_to_compute, simulations, stellar_angular_radius)
%INTERFEROGRAM_SENSITIVITY_MULTIPLE_BRANCHESE Computes the nominal and
%perturbed response functions for interferogram sensitivity analysis, by
%applying a random perturbation to each aperture.
%
% ARGUMENT INPUTS:
%   OPD_n[Np x 1]       OPD in the nominal case for all the Np points [m]
%   OPD[Np x Ns]        OPD in the perturbed cases for all points and
%                       simulations [m]
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
%   2025-04-24 -------- 1.1
%                     - Inputs improved: now, instead of the file, the
%                       function requires information that is already
%                       available in the main script. Order updated.
%   2025-05-20 -------- 1.1.1
%                     - Uses default_arguments for parsing.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    OPD_n
    OPD
    intensities 
    nominal_phases 
    positions 
    combination 
    lambda = 1e-5
    theta = default_arguments("theta")
    maps_to_compute = default_arguments("maps_to_compute")
    simulations = 100
    stellar_angular_radius = default_arguments("stellar_angular_radius")
end

% Compute nominal phase
delta_phi_n = opd2phase(OPD_n(:, 1), lambda);  % Np x 1

% Extract number of simulations
Np = size(OPD, 1);
N = size(positions, 1);
Nt = length(theta);
Ns = size(OPD, 2);  

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