function phase = opd2phase(optical_path, lambda)
%OPD2PHASE Convert optical path difference to phase and wrap to [-pi, pi].
%
% INPUTS:
%   optical_path[matrix]  Matrix where each row corresponds to a series and
%                         each column corresponds to a point in the series.
%   lambda[scalar]        Wavelength used for phase conversion.
%
% OUTPUTS:
%   phase[matrix]         Phase matrix in the same format as optical_path,
%                         wrapped to the range [-pi, pi].
%
% NOTES:
%   - The phase is calculated using the formula: phi = 2 * pi * OP / lambda.
%
% REFERENCES:
%   Lucas Viseur. Development of a performance modeling tool for nulling 
%   interferometry. Université de Liège; 2024. 
%
% VERSION HISTORY:
%   2025-03-28 -------- 1.0
%
% Author: Francesco De Bortolii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the phase
phase = wrapToPi(2 * pi * optical_path / lambda);

end