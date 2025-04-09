function [R, RF] = create_interferogram(positions, A, combination, ...
            phase_perturbations, lambda, theta_range, export_setup)
%CREATE_INTERFEROGRAM_SINGLE_PERTURBATION Computes the interferogram
%response map and response function for a system of apertures with phase
%perturbations.
%
% INPUTS:
%   positions[Nx2]      Matrix of (x, y) positions of the apertures [m]
%   A[Nx1]              Vector of aperture amplitudes [V/m]
%   combination[Nx1]    Vector of beam combiner coefficients [-]
%   phase_perturbations[Np x N] Matrix of phase perturbations for each 
%                       aperture [rad]
% OPTIONAL INPUTS:
%   lambda[1]           Wavelength [m] (default: 1e-5 m)
%   theta_range[Na x 1] Vector of angular coordinates [rad] 
%                       (default: linspace(-700, 700, 1000) mas converted
%
% OUTPUTS:
%   R[Np x Na]          Response map matrix for every point and angle
%   RF[Na x 1]          Response function vector (integrated over point)
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    positions  
    A 
    combination 
    phase_perturbations
    lambda = 1e-5
    theta_range = mas2rad(linspace(-700, 700, 1000))
    export_setup.create_plot = true
    export_setup.interferograms_to_plot = 3
end

% Compute points size
N = size(positions, 1);
Np = size(phase_perturbations, 1);
Na = length(theta_range);

% Allocate space
R = zeros(Np, Na);

% For every angular position compute the response of the instrument
for i = 1:length(theta_range)
    R(:, i) = create_field(theta_range(i), 0, N, lambda, positions, A, ...
            combination, phase_perturbations);
end

% Integrate the response over all the surface to obtain the response
% function for every angular positions.
RF = sum(R, 1);

end