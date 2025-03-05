function [B_LZ, B_EZ] = ...
    compute_background_emissions(lambda, D, B, st_elat, st_dist, ...
                                 st_rad, st_teff, st_lum, TQB)
%COMPUTE_BACKGROUND_EMISSIONS Computes the local and exozodiacal background 
% emissions in ph/s/m.
% 
% INPUTS:
%   lambda              Wavelength [m]
%   nu                  Frequency [Hz]
%   D                   Telescope diameter [m]
%   B                   Telescope baseline [m]
%   st_elat             Ecliptic latitude of viewing direction [rad]
%   st_dist             Distance to the star [m]
%   st_rad              Star radius [m]
%   st_teff             Star temperature [K]
%   st_lum              Star luminosity [L_SUN]
%   TQB                 Throughput [-]
%
% OUTPUTS:
%   B_LZ                Local zodiacal disk emission [ph/s/m]
%   B_EZ                Exozodiacal disk emission [ph/s/m]
%
% NOTES:
%   - This version has been adapted from the original provided by Colin
%     Dandumont and modified to adapt to the structure here used.
%
% VERSION HISTORY:
%   2025-03-05 -------- 1.0b
%
% Authors: 
%  - Colin Dandumont, CSL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physical constants
h = 6.626e-34;     % Planck's constant [J*s]

% Compute local zodiacal contribution (assuming anti-solar direction)
sol_long = pi; %[rad]
B_LZ = compute_local_zodiac_contribution(st_elat, sol_long, lambda); %[Jy/sr]
B_LZ = B_LZ .* 10^(-26) ./ (lambda .* h); %[ph/s/m2/m/sr]

% Compute exozodiacal contribution
B_EZ = compute_exozodiac_contribution(lambda, D, st_dist, st_lum, st_teff, B, st_rad); %[Jy/sr]
B_EZ = B_EZ .* 10^(-26) ./ (lambda .* h); %[ph/s/m2/m/sr]

% Convert to ph/s/m
FOV = pi .* (0.61 .* lambda ./ D).^2; %[sr]
area = (pi .* (D/2).^2); %[m^2]
B_LZ = B_LZ .* TQB .* FOV .* 2 .* area; %[ph/s/m]
B_EZ = B_EZ .* TQB .* FOV .* 2 .* area; %[ph/s/m]

end
