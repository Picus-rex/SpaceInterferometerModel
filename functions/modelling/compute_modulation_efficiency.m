function eta_mod = compute_modulation_efficiency(eta_opt, areas, C_i)
%COMPUTE_MODULATION_EFFICIENCY Compute the asymptotic modulation efficiency
%of an array following the reference by Lay.
%
% INPUTS:
%   eta_opt[Nx1]        Optical efficiency of every aperture.
%   C_i[N1xN2]          Baseline amplitude factors. [-]
%   areas[Nx1]          Areas of every aperture. [m^2]
%
% OUTPUTS:
%   eta_mod[N2]         Asymptotic rms modulation efficiency of array;
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%   Formula 15, 16 on page 8.
%
% NOTES:
%   None
%
% VERSION HISTORY:
%   2025-03-20 -------- 1.0
%   2025-04-30 -------- 1.1
%                     - Expanded the function to multiple points. 
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_mod = zeros(size(C_i, 2), 1);

for i = 1:length(eta_mod)
    num = sqrt(0.5 * sum(C_i(:, i).^2));
    den = sum(eta_opt * areas);
    eta_mod(i) = num / den;
end

end