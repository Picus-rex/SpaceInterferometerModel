function optical_path = phase2opd(phase, lambda)
%PHASE2OPD Convert phase to optical path difference.
%
% INPUTS:
%   phase[matrix]         Phase matrix of any size.
%   lambda[scalar]        Wavelength used for phase conversion.
%
% OUTPUTS:
%   optical_path[matrix]  Optical path difference matrix in the same format
%                         as phase.
%
% NOTES:
%   - The optical path difference is calculated using the formula:
%     OPD = phi * lambda / (2 * pi).
%   - Any information on the actual optical path is however lost! This
%     function only return the OPD to the wrapped amount. 
%
% REFERENCES:
%   Lucas Viseur. Development of a performance modeling tool for nulling
%   interferometry. Université de Liège; 2024.
%
% VERSION HISTORY:
%   2025-04-02 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the optical path difference
optical_path = phase * lambda / (2 * pi);

end
