function x = default_arguments(label)
%DEFAULT_ARGUMENTS Return default arguments for simulations when not
%differently specified.
%
% INPUTS:
%   label[string]       Name of the argument to obtain. Can be:
%
%
% OUTPUTS:
%   x[vary]             Desired argument.
%
% VERSION HISTORY:
%   2025-05-12 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = mas2rad(linspace(-700, 700, 1000));
maps_to_compute = 1 : floor(length(theta) / 3) : length(theta);
stellar_angular_radius = 1.50444523662918e-09;

switch label
    case "theta"
        x = theta;
    case "maps_to_compute"
        x = maps_to_compute;
    case "stellar_angular_radius"
        x = stellar_angular_radius;
    otherwise
        error("Label %s does not exist!", label)
end

end