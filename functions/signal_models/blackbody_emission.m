function Phi = blackbody_emission(T, lambda, delta_lambda, area, angular_dimension, albedo)
% BLACKBODY_EMISSION Computes the photon flux from a blackbody source.
% This function calculates the total photon flux (ph/m^2/s) emitted by a 
% blackbody at temperature T, integrating over a given spectral bandwidth 
% centered around lambda.
%
% INPUTS:
%   T[1]                Temperature of the blackbody [K].
%   lambda[1]           Central wavelength of the band [m].
%   delta_lambda[1]     Bandwidth extension (full width) [m].
%   area[1]             Area of a single collector [m^2].
%   angular_dimension[1]Angular length of the object to study [rad].
%   albedo[1]           (Optional) Albedo factor (0 to 1). Default is 1.
%
% OUTPUTS:
%   Phi[1]              Total photon flux [ph/m^2/s].
%
% NOTES:
%   - This version has been adapted from the original provided by Colin
%     Dandumont and modified to adapt to the structure here used.
%   - MUST INTEGRATE MORE DATA FOR THE STAR AND PLANETS!
%
% VERSION HISTORY:
%   2025-03-06 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
h = 6.62607015e-34;   % Planck's constant [J*s]
c = 2.99792458e8;     % Speed of light [m/s]
k = 1.380649e-23;     % Boltzmann's constant [J/K]

% Default albedo
if nargin < 6
    albedo = 1;
end

% Wavelength range for integration
lambda_min = lambda - delta_lambda / 2;
lambda_max = lambda + delta_lambda / 2;
n_lambda = 100;
lambda_vec = linspace(lambda_min, lambda_max, n_lambda);

% Spectral photon flux density [ph/s/m^2/m/sr]
F_lambda = (2 .* c ./ lambda_vec.^4) ./ (exp((h .* c ./ (lambda_vec .* k .* T))) - 1);

% [ph/s/m^2/m/sr] to [ph/s/m]
F_lambda = F_lambda * area * pi * (angular_dimension)^2;

% Integrate over the bandwidth to get the total photon flux
Phi = trapz(lambda_vec, F_lambda) * albedo;

end
