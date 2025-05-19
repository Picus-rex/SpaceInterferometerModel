function starTable = compute_habitable_zone(starTable)
%COMPUTE_HABITABLE_ZONE Calculates HZ limits (AU and rad) for main-sequence 
% stars.
%
% INPUTS:
%   startable[table]    Unique stars of the catalogue. With fields:
%       Nuniverse[Nex1] Universe associated to that planet.
%       Nstar[Ne x 1]   Univoque number of the star in the universe.
%       Rs[Ne x 1]      Radius of the star. [R_Sun]
%       Ts[Ne x 1]      Effective temperature of star. [K]
%       Ds[Ne x 1]      Distance to Sun of host star. [pc]
%       Ms[Ne x 1]      Mass of the star. [M_Sun]
%       Stype[Ne x 1]   Type of the host star (like 'F', 'G', 'K', 'M').
%       RA[Ne x 1]      Right ascension of the star [deg]
%       Dec[Ne x 1]     Declination of the star [deg]
%
% OUTPUTS:
%   startable[table]    Unique stars of the catalogue. With fields:
%       Nuniverse[Nsx1] Universe associated to that planet.
%       Nstar[Ns x 1]   Univoque number of the star in the universe.
%       Rs[Ns x 1]      Radius of the star. [R_Sun]
%       Ts[Ns x 1]      Effective temperature of star. [K]
%       Ds[Ns x 1]      Distance to Sun of host star. [pc]
%       Ms[Ns x 1]      Mass of the star. [M_Sun]
%       Stype[Ns x 1]   Type of the host star (like 'F', 'G', 'K', 'M').
%       RA[Ns x 1]      Right ascension of the star [deg]
%       Dec[Ns x 1]     Declination of the star [deg]
%       HZ_min[Ns x 1]  Minimum position of habitable zone [AU]
%       HZ_max[Ns x 1]  Maximum position of habitable zone [AU]
%       AngSep_HZ_min[Nsx1] Min pos of HZ in ang. sep. from Earth [rad]
%       AngSep_HZ_max[Nsx1] Max pos of HZ in ang. sep. from Earth [rad]
%
% NOTES:
%   This function requires that the coefficients file Kopparapu_HZ.csv
%   exists and is located in the others folder.
%
% REFERENCES:
%   Kopparapu RK, Ramirez R, Kasting JF, Eymet V, Robinson TD, Mahadevan S,
%   et al. Habitable zones around main-sequence stars: new estimates. ApJ.
%   2013 Feb 26;765(2):131.
%
% VERSION HISTORY:
%   2025-05-19 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
T_sun = 5778; % K

% Filter temperature range and keep only valid stars
validIdx = starTable.Ts >= 2600 & starTable.Ts <= 7200;
numScrapped = sum(~validIdx);
fprintf('Filtered out %d stars outside temperature range 2600K–7200K.\n', numScrapped);
starTable = starTable(validIdx, :);

% Load coefficients from external file
hzData = readmatrix('others/Kopparapu_HZ.csv'); 
coeffs_runaway = hzData(1,:); % Inner HZ
coeffs_maxgreen = hzData(2,:); % Outer HZ

% Preallocate result columns
n = height(starTable);
HZ_min = zeros(n,1);
HZ_max = zeros(n,1);
AngSep_HZ_min = zeros(n,1);
AngSep_HZ_max = zeros(n,1);

for i = 1:n
    T_eff = starTable.Ts(i);
    R_star = starTable.Rs(i);
    d_pc = starTable.Ds(i);

    % Stellar luminosity (L/L_sun)
    L_star = (R_star^2) * (T_eff / T_sun)^4;

    % Temperature offset
    T_star = T_eff - 5780;

    % S_eff for inner (runaway greenhouse)
    S_eff_in = coeffs_runaway(1) + ...
               coeffs_runaway(2)*T_star + ...
               coeffs_runaway(3)*T_star^2 + ...
               coeffs_runaway(4)*T_star^3 + ...
               coeffs_runaway(5)*T_star^4;

    % S_eff for outer (maximum greenhouse)
    S_eff_out = coeffs_maxgreen(1) + ...
                coeffs_maxgreen(2)*T_star + ...
                coeffs_maxgreen(3)*T_star^2 + ...
                coeffs_maxgreen(4)*T_star^3 + ...
                coeffs_maxgreen(5)*T_star^4;

    % Compute HZ distances in AU
    HZ_min(i) = sqrt(L_star / S_eff_in);
    HZ_max(i) = sqrt(L_star / S_eff_out);

    % Convert AU to angular separation in radians
    AngSep_HZ_min(i) = HZ_min(i) / (d_pc * 206265); % 1 pc = 206265 AU
    AngSep_HZ_max(i) = HZ_max(i) / (d_pc * 206265);
end

% Append to table
starTable.HZ_min = HZ_min;
starTable.HZ_max = HZ_max;
starTable.AngSep_HZ_min = AngSep_HZ_min;
starTable.AngSep_HZ_max = AngSep_HZ_max;
end
