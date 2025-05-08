function [ratio, rejection] = compute_nulling_ratio(N, amplitudes, ...
                                    phases, positions, theta_star, lambda)
%COMPUTE_NULLING_RATIO Computation of the nulling ratio and rejection
%fraction for a given array.
%
% INPUTS:
%   N[1]                Number of apertures. [-]
%   amplitudes[Nx1]     Amplitudes associated to each aperture.
%   phases  [Nx1]       Vector of phase shifts for each aperture OR
%           [Np x N]    Matrix of phase shifts associated to apertures and
%                       poisitons points on the aperture.
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation. [m]
%
% OUTPUTS:
%   ratio   [1] or [Np] Nulling ratio of the array.
%   rejection[1]or [Np] Rejection fraction of the array.
%
% REFERENCES:
%   Defrère D. Characterizing extrasolar planetary systems 
%   using infrared interferometry and Absil O. Astrophysical studies of 
%   extrasolar planetary systems using infrared interferometric techniques 
%   [Internet]. 2006. Available from: 
%   https://theses.hal.science/tel-00124720v1/file/Absil06_thesis.pdf
%
% NOTES:
%   - Perturbations need to be introduced as phase shifts above.
%
% VERSION HISTORY:
%   2025-02-28 -------- 1.0
%   2025-03-11 -------- 1.1
%                     - Use of the base formula to address the possible
%                       amplitudes interferences.
%   2025-04-09 -------- 1.2
%                     - Introduced possibility to accept matrices for
%                       phases in the case of spatial positions considered.
%                     - rejection computed only if required to speed up.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate space for geometric leakage depending on the shape
if isvector(phases)
    N_g = 0;
else
    N_g = zeros(size(phases, 1), 1);
end

for j = 1:N
    for k = 1:N
        
        if j ~= k
            % Baseline b_jk
            b_jk = norm(positions(j,:) - positions(k,:));
    
            % Stellar brightness B_star_jk
            B_star_jk = 2 * besselj(1, 2 * pi * b_jk * theta_star / lambda) / (2 * pi * b_jk * theta_star / lambda);
        else
            B_star_jk = 1;
        end
        
        if isvector(phases)
            N_g = N_g + amplitudes(j) * amplitudes(k) * cos(phases(j) - phases(k)) * B_star_jk;
        else
            N_g = N_g + amplitudes(j) * amplitudes(k) * cos(phases(:, j) - phases(:, k)) * B_star_jk;
        end
    end
end

% Denominator
sum_N_star = sum(amplitudes.^2);

% Nulling ratio and rejection factor (only if required)
ratio = N_g / sum_N_star;
if nargout == 2
    rejection = ratio.^(-1);
end

end