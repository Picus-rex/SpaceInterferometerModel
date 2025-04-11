function rad = arcsec2rad(arcsec)
%RAD2MAS Converts angular measurements from milliarcseconds to radians.
%
% ARGUMENT INPUTS:
%   arcsec[any]         Angular measurement(s) in arcseconds
%
% OUTPUTS:
%   rad[any]            Angular measurement(s) in radians
%
% REFERENCES:
%   Conversion factor: 1 rad = 206265 arcseconds
%
% VERSION HISTORY:
%   2025-04-11 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rad = ((pi/180)/3600) * arcsec;

end