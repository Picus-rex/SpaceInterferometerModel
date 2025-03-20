function eta_mod = compute_modulation_efficiency(eta_opt, areas, C_i)
%COMPUTE_MODULATION_EFFICIENCY Compute the asymptotic modulation efficiency
%of an array following the reference by Lay.
%
% INPUTS:
%   eta_opt[Nx1]        Optical efficiency of every aperture.
%   C_i[N1x1]           Baseline amplitude factors. [-]
%   areas[Nx1]          Areas of every aperture. [m^2]
%
% OUTPUTS:
%   eta_mod[1]          Asymptotic rms modulation efficiency of array;
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   None
%
% VERSION HISTORY:
%   2025-03-20 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Using formula 15, 16 on page 8.
num = sqrt(0.5 * sum(C_i.^2));
den = sum(eta_opt * areas);
eta_mod = num / den;

end