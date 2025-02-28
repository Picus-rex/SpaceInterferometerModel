function [T_standard, T_chopped, T_real, T_real_std, data] = ...
    compute_response_function(data)
%COMPUTE_RESPONSE_FUNCTION Computes the response function and the phase 
% chopping response function for a system of N apertures.
%
% ARGUMENT INPUTS:
%   lambda[1]           Wavelength [m]
%   N[1]                Number of apertures [-]
%   positions[Nx2]      Matrix of (x, y) positions of the apertures
%   A[Nx1]              Vector of aperture amplitudes [V/m]
%   phase_shifts[Nx1]   Vector of phase shifts for each aperture
%   combination[Nx1]    Vector of beam combiner coefficients
%   theta_x, theta_y    Meshgrid of angular coordinates
%   N_MC[1]             Number of MC simulations to perform
%   sigma_phase[1]      Standard deviation for the phase errors [rad]
%   sigma_intensity[1]  Standard deviation for relative intensity errors 
%   sigma_pol[1]        Standard deviation for polarisation-induced phase errors [rad]
%   environment[struct] External perturbations (see configuration)
%
%   OR
%
%   data[struct]        For integrated development, all the inputs can be
%                       grouped into a single struct
%
% OUTPUTS:
%   T_standard          Normalised response function 
%   T_chopped           Normalised phase chopping response function
%   data[struct]        For integrated development, group all the results                
% 
%   With specific inputs:
%
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
%   2025-02-28 -------- 2.0
%                     - Use of name-value pairs for inputs.
%                     - Add external sensitivity.
%                     - Not retrocompatible with older versions.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data.lambda (1, 1) double = 10e-6
    data.N (1, 1) double = 4
    data.positions (:, 2) {mustBeNumeric}
    data.A (:, 1) double = [1; 1; 1; 1]
    data.phase_shifts (:, 1) double
    data.combinations (:, 1) double
    data.theta_x (:, :) double
    data.theta_y (:, :) double

    data.N_MC (1, 1) double = 0
    data.sigma_phase (1, 1) double = 0
    data.sigma_intensity (1, 1) double = 0
    data.sigma_pol (1, 1) double = 0

    data.environment struct = struct()

    data.data struct
end

% If a single struct is passed (in the integrated module), extract
% everything to continue.
if ~isfield(data, "data")
    
    data.instrument = {
        "wavelength": data.lambda, ...
        "apertures" : data.N, ...
        "positions" : data.positions, ...
        "intensities": data.A, ...
        "phase_shifts": data.phase_shifts, ...
        "combination": data.combination
    };

    data.simulation = {
        "theta_x" : data.theta_x, ...
        "theta_y" : data.theta_y, ...
        "monte_carlo_iterations" : data.N_MC
    };

    if data.N_MC > 0
         
        data.instrument.instrumental_leakage.phase = data.sigma_phase;
        data.instrument.instrumental_leakage.intensity = data.sigma_intensity;
        data.instrument.instrumental_leakage.polarisation = data.sigma_pol;

    end


else
    data = data.data;
end

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
    eps_MC = zeros(N_MC, 1);
    
    % Monte Carlo simulation loop
    for mc = 1:N_MC
    
        % Perturb amplitude
        A_err = A .* (1 + sigma_intensity * randn(size(A)));
        
        % Perturb phase: independent errors for phase and polarization
        phase_err = phase_shifts + sigma_phase * randn(size(phase_shifts)) ...
            + sigma_pol * randn(size(phase_shifts));

        R_MC(:,:,mc) = create_field(theta_x, theta_y, N, lambda, ...
            positions, A_err, combination, phase_err);

        if ~isempty(fieldnames(environment))
            eps_MC(mc) = add_external_sensitivity(instrument, environment);
        end

    end
    
    % CMean and standard deviation
    R_mean = mean(R_MC, 3);
    R_std  = std(R_MC, 0, 3);
    epsilon = mean(eps_MC);
    epsilon_std = std(eps_MC, 0);
    
    % Normalize the responses
    norm_factor = max(R_mean(:));
    
    T_real = R_mean * (1 + epsilon) / norm_factor;
    T_real_std = R_std / norm_factor;

end

% Group results
data.simulation.T_standard = T_standard;
data.simulation.T_chopped = T_chopped;
data.simulation.T_real = T_real;
data.simulation.T_read_std = T_real_std;
data.simulation.epsilon_std = epsilon_std;

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