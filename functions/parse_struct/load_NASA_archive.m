function T = load_NASA_archive()
%LOAD_NASA_ARCHIVE Loads a table of confirmed exoplanet data from the NASA
%archive.
%
% OUTPUTS:
%   T[table]            Table containing confirmed exoplanet data with the 
%                       following fields:
%       pl_name[string]              Exoplanet rame
%       st_name[string]              Star name
%       Tp[double]                   Equilibrium temp. of the planet [K]
%       Rp[double]                   Radius of the planet [Earth radii]
%       Stype[string]                Type of the star 
%       Ts[double]                   Effective temp. of the host star [K]
%       Rs[double]                   Radius of the host star [Solar radii]
%       AngSep[double]               Ang. sep. between pl. and star [rad]
%       ap[double]                   Sem-maj. ax. of planet's orbit [AU]
%
% REFERENCES:
%   NASA exoplanets archive
%
% VERSION HISTORY:
%   2025-04-10 -------- 1.0
%   2025-04-11 -------- 1.1
%                     - Updated new archive with more reliable content and
%                       less guesses required in the computation.
%                     - Removed references to external content for
%                       bolometric corrections. 
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load table from CSV, skipping header lines
opts = detectImportOptions('others/NASA_archive.csv', 'NumHeaderLines', 32);
data = readtable('others/NASA_archive.csv', opts);

% Remove all the data that is incomplete, leaving only valid exoplanets for
% each simulation
data = data(~isnan(data.pl_masse) & ~isnan(data.pl_dens) ...
    & ~isnan(data.sy_dist) & ~isnan(data.pl_eqt) ...
    & ~isnan(data.st_teff) & ~isnan(data.st_rad) ...
    & ~isnan(data.pl_rade) & ~isnan(data.pl_orbsmax), :);

angular_sep_rad = data.pl_orbsmax ./ pc2au(data.sy_dist); 

% Create new output table
T = table(data.pl_name, data.hostname, data.pl_eqt, data.pl_rade, ...
    data.st_spectype, data.st_teff, data.st_rad, angular_sep_rad, ...
    data.pl_orbsmax, 'VariableNames', {'pl_name', 'st_name', 'Tp', ...
    'Rp', 'Stype', 'Ts', 'Rs', 'AngSep', 'ap'});

end