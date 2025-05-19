function summaryTable = summarise_HZ_by_type(starTable)
%SUMMARIZE_HZ_BY_TYPE Computes HZ statistics per spectral type.
%
% INPUTS:
%   starTable[table]    Unique stars of the catalogue. With fields:
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
% OUTPUTS:
%   summaryTable[table] Table that, for each stellar type, contains:
%       Stype[Nt x 1]   Type of the host star (like 'F', 'G', 'K', 'M').
%       HZ_min_10[Ntx1] Minimum HZ 10% percentile [AU]
%       HZ_min_90[Ntx1] Minimum HZ 90% percentile [AU]
%       HZ_min_mean[Ntx1] Minimum HZ mean value [AU]
%       ASHZ_min_10[Ntx1] Minimum HZ 10% percentile [rad]
%       ASHZ_min_90[Ntx1] Minimum HZ 90% percentile [rad]
%       ASHZ_min_mean[Ntx1] Minimum HZ mean value [rad]
%       HZ_max_10[Ntx1] Maximum HZ 10% percentile [AU]
%       HZ_max_90[Ntx1] Maximum HZ 90% percentile [AU]
%       HZ_max_mean[Ntx1] Maximum HZ mean value [AU]
%       ASHZ_max_10[Ntx1] Maximum HZ 10% percentile [rad]
%       ASHZ_max_90[Ntx1] Maximum HZ 90% percentile [rad]
%       ASHZ_max_mean[Ntx1] Maximum HZ mean value [rad]
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

% Unique spectral types
stypes = unique(starTable.Stype);
nTypes = numel(stypes);

% Preallocate arrays
Stype = cell(nTypes,1);
HZ_min_10 = zeros(nTypes,1);
HZ_min_90 = zeros(nTypes,1);
HZ_min_mean = zeros(nTypes,1);

ASHZ_min_10 = zeros(nTypes,1);
ASHZ_min_90 = zeros(nTypes,1);
ASHZ_min_mean = zeros(nTypes,1);

HZ_max_10 = zeros(nTypes,1);
HZ_max_90 = zeros(nTypes,1);
HZ_max_mean = zeros(nTypes,1);

ASHZ_max_10 = zeros(nTypes,1);
ASHZ_max_90 = zeros(nTypes,1);
ASHZ_max_mean = zeros(nTypes,1);

% Loop over each spectral type
for i = 1:nTypes
    typeMask = strcmp(starTable.Stype, stypes{i});
    HZmin = starTable.HZ_min(typeMask);
    HZmax = starTable.HZ_max(typeMask);
    ASHZmin = starTable.AngSep_HZ_min(typeMask);
    ASHZmax = starTable.AngSep_HZ_max(typeMask);

    Stype{i} = stypes{i};

    HZ_min_10(i) = prctile(HZmin, 10);
    HZ_min_90(i) = prctile(HZmin, 90);
    HZ_min_mean(i) = mean(HZmin);

    ASHZ_min_10(i) = prctile(ASHZmin, 10);
    ASHZ_min_90(i) = prctile(ASHZmin, 90);
    ASHZ_min_mean(i) =  mean(ASHZmin);

    HZ_max_10(i) = prctile(HZmax, 10);
    HZ_max_90(i) = prctile(HZmax, 90);
    HZ_max_mean(i) = mean(HZmax);

    ASHZ_max_10(i) = prctile(ASHZmax, 10);
    ASHZ_max_90(i) = prctile(ASHZmax, 90);
    ASHZ_max_mean(i) =  mean(ASHZmax);
end

% Construct summary table
summaryTable = table(Stype, ...
                     HZ_min_10, HZ_min_mean, HZ_min_90, ...
                     ASHZ_min_10, ASHZ_min_mean, ASHZ_min_90, ...
                     HZ_max_10, HZ_max_mean, HZ_max_90, ...
                     ASHZ_max_10, ASHZ_max_mean, ASHZ_max_90);
end
