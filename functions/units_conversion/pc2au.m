function au = pc2au(pc)
%PC2AU Converts parsec distances to astronomical units.
%
% ARGUMENT INPUTS:
%   pc[any]             Distance measurement(s) in parsec.
%
% OUTPUTS:
%   au[any]             Distances measurement(s) in astronomical units.
%
% REFERENCES:
%   Conversion factor: 1 pc = 206265 au
%
% VERSION HISTORY:
%   2025-04-11 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

au = 206265 * pc;

end