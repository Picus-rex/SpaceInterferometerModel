function [T_standard, T_chopped, T_real, T_real_std] = ...
    compute_response_function(lambda, N, positions, A, phase_shifts, ...
    combination, theta_x, theta_y, N_MC, sigma_phase, sigma_intensity, sigma_pol)
%COMPUTE_RESPONSE_FUNCTION Computes the response function and the phase 
% chopping response function for a system of N apertures.
%
% INPUTS:
%   lambda[1]           Wavelength [m]
%   N[1]                Number of apertures [-]
%   positions[Nx2]      Matrix of (x, y) positions of the apertures
%   A[Nx1]              Vector of aperture amplitudes [V/m]
%   phase_shifts[Nx1]   Vector of phase shifts for each aperture
%   combination[Nx1]    Vector of beam combiner coefficients
%   theta_x, theta_y    Meshgrid of angular coordinates
%  Optional:
%   N_MC[1]             Number of MC simulations to perform
%   sigma_phase[1]      Standard deviation for the phase errors [rad]
%   sigma_intensity[1]  Standard deviation for relative intensity errors 
%   sigma_pol[1]        Standard deviation for polarisation-induced phase errors [rad]
%
% OUTPUTS:
%   T_standard          Normalised response function 
%   T_chopped           Normalised phase chopping response function
%  If optional arguments have been specified:
%   T_real              Normalized mean intensity response with errors
%   T_real_std          Normalized 1-sigma envelope of the error response
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%   Absil O. Astrophysical studies of extrasolar planetary systems using 
%   infrared interferometric techniques [Internet]. 2006.
%
% VERSION HISTORY:
%   2025-02-11 -------- 1.0
%   2025-02-27 -------- 1.1
%                     - Added a second function to compute the field,
%                       simplifying the code.
%                     - Added compatibility for MC analysis and error
%                       definition. The function can still be used as it
%                       is.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dealing with partial outputs
if nargin <= 8
    N_MC = 0;
end

% If a single struct is passed (in the integrated module), extract
% everything to continue.
if isstruct(lambda)

    data = lambda;
    instrument =  data.instrument;
    simulation =  data.simulation;
    environment = data.environment;

    lambda = instrument.wavelength;
    N = instrument.apertures;
    positions = instrument.positions; 
    A = instrument.intensities;
    phase_shifts = instrument.phase_shifts;
    combination = instrument.combination;
    
    theta_x = simulation.theta_x;
    theta_y = simulation.theta_y;

    if isfield(data.simulation, "monte_carlo_iterations") && simulation.monte_carlo_iterations > 0
        N_MC = simulation.monte_carlo_iterations;
        
        sigma_phase = instrument.instrumental_leakage.phase;
        sigma_intensity = instrument.instrumental_leakage.intensity;
        sigma_pol = instrument.instrumental_leakage.polarisation;
    else
        N_MC = 0;
    end

end

% Grid dimensions
nx = size(theta_x, 1);
ny = size(theta_y, 1);

R_standard = create_field(theta_x, theta_y, N, lambda, positions, A, ...
    combination, phase_shifts);
R_negated = create_field(theta_x, theta_y, N, lambda, positions, A, ...
    combination, -phase_shifts);

% Compute Phase Chopping Response
R_chopped = 0.5 * abs(R_standard - R_negated);

% Normalize to max value
T_standard = R_standard / max(R_standard(:));
T_chopped = R_chopped / max(R_chopped(:));

if N_MC > 0

    % Preallocate for MC analysis 
    R_MC = zeros(nx, ny, N_MC);
    
    % Monte Carlo simulation loop
    for mc = 1:N_MC
    
        % Perturb amplitude
        A_err = A .* (1 + sigma_intensity * randn(size(A)));
        
        % Perturb phase: independent errors for phase and polarization
        phase_err = phase_shifts + sigma_phase * randn(size(phase_shifts)) ...
            + sigma_pol * randn(size(phase_shifts));

        R_MC(:,:,mc) = create_field(theta_x, theta_y, N, lambda, ...
            positions, A_err, combination, phase_err);
    end
    
    % CMean and standard deviation
    R_mean = mean(R_MC, 3);
    R_std  = std(R_MC, 0, 3);
    
    % Normalize the responses
    norm_factor = max(R_mean(:));
    
    T_real = R_mean / norm_factor;
    T_real_std = R_std / norm_factor;

end

end





function R = create_field(theta_x, theta_y, N, lambda, positions, A, ...
    combination, phase_shifts)

E = zeros(size(theta_x));

for k = 1:N

    xk = positions(k,1);
    yk = positions(k,2);
    
    % Phase term
    phase_k = 2 * pi * (xk * theta_x + yk * theta_y) / lambda + phase_shifts(k);
    
    % Add contribution from each aperture
    E = E + combination(k) * A(k) * exp(1i * phase_k);

end

R = abs(E).^2;

end