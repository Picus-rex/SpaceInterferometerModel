function startable = extract_stars(exotable)
%EXTRACT_STARS Extracts unique host stars from an exoplanet table.
%
% INPUTS:
%   exotable[table]     Exoplanets table from P-Pop with fields:
%       Nuniverse[Nex1] Universe associated to that planet.
%       Rs[Ne x 1]      Radius of the star. [R_Sun]
%       Ts[Ne x 1]      Effective temperature of star. [K]
%       Ds[Ne x 1]      Distance to Sun of host star. [pc]
%       Ms[Ne x 1]      Mass of the star. [M_Sun]
%       Stype[Ne x 1]   Type of the host star (like 'F', 'G', 'K', 'M').
%       RA[Ne x 1]      Right ascension of the star [deg]
%       Dec[Ne x 1]     Declination of the star [deg]
%       Nstar[Ne x 1]   Univoque number of the star in the universe.
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
%
% VERSION HISTORY:
%   2025-05-19 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure required variables exist
requiredVars = {'Nuniverse','Nstar','Ds','Rs','Ts','Ms','Stype','RA','Dec'};
for var = requiredVars
    if ~ismember(var{1}, exotable.Properties.VariableNames)
        error('Missing required variable: %s', var{1});
    end
end

% Create a unique ID for each star in a given universe
starKeys = strcat(string(exotable.Nuniverse), '_', string(exotable.Nstar));

% Identify the first occurrence of each unique star (and keep order)
[~, uniqueIdx] = unique(starKeys, 'stable');

% Extract only the first occurrence of each star
startable = exotable(uniqueIdx, {'Nuniverse','Nstar','Ds','Rs','Ts','Ms','Stype','RA','Dec'});

end
