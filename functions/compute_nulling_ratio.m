function [ratio, rejection] = ...
    compute_nulling_ratio(N, phases, positions, theta_star, lambda)
%COMPUTE_NULLING_RATIO Computation of the nulling ratio and rejection
%fraction for a given array.
%
% INPUTS:
%   N[1]                Number of apertures. [-]
%   phases[Nx1]         Vector of phase shifts for each aperture
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation. [m]
%
% OUTPUTS:
%   ratio[1]            Nulling ratio of the array.
%   rejection[1]        Rejection fraction of the array.
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
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum = N;  
for j = 1:N
    for k = 1:N        
        if j ~= k
            Delta_x = positions(j,1) - positions(k,1);
            Delta_y = positions(j,2) - positions(k,2);
            dcphi = cos(phases(j) - phases(k));
            b_jk = sqrt(Delta_x^2 + Delta_y^2);
            arg = 2*pi * b_jk * theta_star / lambda;

            sum = sum + dcphi * (2 * besselj(1, arg) / arg);
        end
    end
end

ratio = sum / N;

rejection = 1/ratio;

end