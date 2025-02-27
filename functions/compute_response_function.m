function [T_standard, T_chopped] = compute_response_function(lambda, N, ...
    positions, A, phase_shifts, combination, theta_x, theta_y)
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
%
% OUTPUTS:
%   T_standard          Normalised response function 
%   T_chopped           Normalised phase chopping response function
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005; 
%
% VERSION HISTORY:
%   2025-02-11 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Standard and Negated Response Function
E_standard = zeros(size(theta_x));
E_negated = zeros(size(theta_x));

for k = 1:N

    xk = positions(k,1);
    yk = positions(k,2);
    
    % Phase term
    phase_k = 2 * pi * (xk * theta_x + yk * theta_y) / lambda + phase_shifts(k);
    phase_k_negated = 2 * pi * (xk * theta_x + yk * theta_y) / lambda - phase_shifts(k);
    
    % Add contribution from each aperture
    E_standard = E_standard + combination(k) * A(k) * exp(1i * phase_k);
    E_negated = E_negated + combination(k) * A(k) * exp(1i * phase_k_negated);

end

% Intensity
R_standard = abs(E_standard).^2;
R_negated = abs(E_negated).^2;

% Compute Phase Chopping Response
R_chopped = 0.5 * abs(R_standard - R_negated);

% Normalize to max value
T_standard = R_standard / max(R_standard(:));
T_chopped = R_chopped / max(R_chopped(:));

end
