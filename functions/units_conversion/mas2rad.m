function rad = mas2rad(mas)
%RAD2MAS Converts angular measurements from milliarcseconds to radians.
%
% ARGUMENT INPUTS:
%   mas[any]            Angular measurement(s) in milliarcseconds (mas)
%
% OUTPUTS:
%   rad[any]            Angular measurement(s) in radians
%
% REFERENCES:
%   Conversion factor: 1 radian = 206264806 milliarcseconds
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rad = ((pi/180)/3600 / 1e3) * mas;

end