function mas = rad2mas(rad)
%RAD2MAS Converts angular measurements from radians to milliarcseconds.
%
% ARGUMENT INPUTS:
%   rad[any]            Angular measurement(s) in radians
%
% OUTPUTS:
%   mas[any]            Angular measurement(s) in milliarcseconds (mas)
%
% REFERENCES:
%   Conversion factor: 1 radian = 206264806 milliarcseconds
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mas = ((pi/180)/3600 / 1e3)^(-1) * rad;

end