function [RFs, RFn, R, ratio, eta_mod] = interferogram_sensitivity(OPD_n, OPD, ...
    intensities, nominal_phases, positions, combination, areas, ...
    lambda, theta, maps_to_compute, eta_opt, stellar_angular_radius)
%INTERFEROGRAM_SENSITIVITY Computes the nominal and perturbed response
%functions for interferogram sensitivity analysis.
%
% ARGUMENT INPUTS:
%   OPD_n[Np x 1]       OPD in the nominal case for all the Np points [m]
%   OPD[Np x Ns]        OPD in the perturbed cases for all points and
%                       simulations [m]
%   intensities[Nx1]    Vector of aperture intensities [V/m]
%   nominal_phases[Nx1] Vector of nom. phase shifts for each aperture [rad]
%   positions[Nx2]      Matrix of (x, y) positions of the apertures [m]
%   combination[Nx1]    Vector of beam combiner coefficients [-]
%   surfaces[Nx1]       Surfaces associated to each apertures [m^2]
%   lambda[1]           Wavelength [m] (default: 1e-5 m)
%   theta[Nt x 1]       Vector of angular coordinates [rad] 
%                       (default: mas2rad(linspace(-700, 700, 1000)))
%   maps_to_compute[1xNm] Indices of the maps to include in the R 3D matrix
%                       (default: [1, 251, 501, 751])
%   eta_opt[1]          Efficiency of optical line (default: 0.25)
%   stellar_angular_radius[1] Angular radius of the star [rad] 
%                       (default: Proxima Centauri value)
%
% OUTPUTS:
%   RFs[Ns x Nt]        Matrix of perturbed response functions
%   RFn[1 x Nt]         Vector of nominal response function
%   R[Np x Nt x Ns]     3D matrix of the response map over apertures
%   ratio[Np x Ns]      Nulling ratio for every point for every simulation
%   eta_mod[Np x Ns]    Modulation eff. for every point for every sim. 
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%   2025-04-16 -------- 1.1
%                     - Inputs improved: now, instead of the file, the
%                       function requires information that is already
%                       available in the main script.
%   2025-04-22 -------- 1.1.1
%                     - Fix in the order of inputs.
%   2025-04-30 -------- 1.2
%                     - Added modulation efficiency as output.
%                     - Inputs reorganised.
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
    areas
    lambda = 1e-5
    theta = mas2rad(linspace(-700, 700, 1000))
    maps_to_compute = 1 : floor(length(theta) / 4) : length(theta)
    eta_opt = 0.1
    stellar_angular_radius = 1.50444523662918e-09
end

% Compute nominal phase
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
eta_mod = zeros(Np, Ns);

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
    
    if nargout > 4
        % Export unique_baselines for modulation efficiency
        C_i = classify_baselines_matrix(intensities, positions, phase_perturbations);
        eta_mod(:, i) = compute_modulation_efficiency(eta_opt, areas, C_i);
    end

end

end