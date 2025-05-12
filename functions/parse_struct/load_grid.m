function points = load_grid(direction)
%LOAD_GRID Automatic grid for computations when not specified. This
%function is part of the arguments default options.
%
% INPUTS:
%   direction           Either "x" or "y".
%
% OUTPUTS:
%   points              Meshgrid results in the respective direction.
%
% VERSION HISTORY:
%   2025-04-02 -------- 1.0
%   2025-05-12 -------- 1.1
%                     - Corrected interva in accordance with other
%                       simulations.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generation of the meshgrid
theta_range = mas2rad(linspace(-700, 700, 1000)); 
[points_x, points_y] = meshgrid(theta_range, theta_range);

switch direction
    case "x"
        points = points_x;
    case "y"
        points = points_y;
end

end