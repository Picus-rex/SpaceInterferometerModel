function R = create_field(theta_x, theta_y, N, lambda, positions, A, ...
    combination, phase_shifts)
%CREATE_FIELD Computes the complex field distribution for a system of N
%apertures to see the fringes patterns.
%
% ARGUMENT INPUTS:
%   theta_x[matrix]     Meshgrid of ang. coordinates in x-direction [rad]
%   theta_y[matrix]     Meshgrid of ang. coordinates in y-direction [rad]
%   N[1]                Number of apertures [-]
%   lambda[1]           Wavelength [m]
%   positions[Nx2]      Matrix of (x, y) positions of the apertures [m]
%   A[Nx1]              Vector of aperture amplitudes [V/m]
%   combination[Nx1]    Vector of beam combiner coefficients [-]
%   phase_shifts[Nx1]   Vector of phase shifts for each aperture [rad]
%
% OUTPUTS:
%   R[matrix]           Complex field distribution matrix.
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%   Absil O. Astrophysical studies of extrasolar planetary systems using
%   infrared interferometric techniques [Internet]. 2006.
%
% VERSION HISTORY:
%   2025-02-27 -------- 1.0
%   2025-04-09 -------- 1.1
%                     - Moved to external file
%                     - Differentiate the type of provided input
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Differentiate between the case in which the independent variable is the
% angular position or the ray position. 
if isvector(phase_shifts)
    E = zeros(size(theta_y, 1), size(theta_x, 2));
else
    E = zeros(size(phase_shifts, 1), 1);
end

for k = 1:N

    xk = positions(k,1);
    yk = positions(k,2);
    
    % Phase term (differentiate into the two cases)
    if isvector(phase_shifts)
        phase_term = phase_shifts(k);
    else
        phase_term = phase_shifts(:, k);
    end
    phase_k = 2 * pi * (xk * theta_x + yk * theta_y) / lambda + phase_term;
    
    % Add contribution from each aperture
    E = E + combination(k) * A(k) * exp(1i * phase_k);

end

R = abs(E).^2;

end