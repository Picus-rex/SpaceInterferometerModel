function [RFs, RFn, R, ratio] = interferogram_sensitivity(nominal_file, ...
    perturbed_file, intensities, nominal_phases, positions, combination,...
    lambda, theta, stellar_angular_radius, maps_to_compute)
%INTERFEROGRAM_SENSITIVITY Computes the nominal and perturbed response
%functions for interferogram sensitivity analysis.
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
end

% Load data from CODE V and compute nominal phase
[~, OPD, ~, ~] = load_opd(perturbed_file);
[~, OPD_n, ~, ~] = load_opd(nominal_file);
delta_phi_n = opd2phase(OPD_n(:, 1), lambda);

% Extract number of simulations
Ns = size(OPD, 2);
Np = size(OPD, 1);
N = size(positions, 1);
Nt = length(theta);

% Obtain phase vector nominal and compute response function
phases_nominal = (nominal_phases + delta_phi_n);
[~, RFn] = create_interferogram(positions, intensities, combination, ...
                                phases_nominal, lambda, theta);

% Allocate space
R = zeros(Np, Nt, length(maps_to_compute));
RFs = zeros(Ns, Nt);
ratio = zeros(Np, Ns);

% Perturbate a single branch
jR = 1;
for i = 1:Ns
    
    % Extract the phase of the ith branch
    delta_phi = opd2phase(OPD(:, i), lambda);

    % Allocate space for perturbations
    phase_perturbations = zeros(Np, N);

    % Perturb first arm, leave others unchanged
    for k = 1:N
        if k == 1
            phase_perturbations(:, k) = nominal_phases(k) + delta_phi;
        else
            phase_perturbations(:, k) = nominal_phases(k) + delta_phi_n;
        end
    end

    % Compute the response function for all the points in the analysis
    [Rm, RFs(i, :)] = create_interferogram(positions, intensities, ...
                    combination, phase_perturbations, lambda, theta);
    
    if ismember(i, maps_to_compute)
        R(:, :, jR) = Rm;
        jR = jR + 1;
    end

    % Compute nulling ratio for each point on the screen
    ratio(:, i) = compute_nulling_ratio(N, intensities, phase_perturbations, ...
        positions, stellar_angular_radius, lambda);

end

end