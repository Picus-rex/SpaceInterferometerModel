function intensities = area2intensity(surfaces, eta_opt, eta_bc)
%AREA2INTENSITY Convert the surfaces of the apertures to the intensities
%using Lay's relation. 
%
% INPUTS:
%   surfaces[Nx1]       Vector of surfaces associated to each aperture.
%                       [m^2]
%   eta_opt[1]          Optical transmission efficiency. [-]
%   eta_bc[1]           Beam combiner efficiency. [-]
%
% OUTPUTS:
%   intensities[1xN]    Vector of intensity associated to each aperture.
%                       [m^2]
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005; 
%
% VERSION HISTORY:
%   2025-05-08 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diameters = 2 * sqrt(surfaces / pi);
intensities = sqrt(pi * diameters.^2 .* eta_opt * eta_bc);

end