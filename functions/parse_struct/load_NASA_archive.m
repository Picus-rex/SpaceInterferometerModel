function T = load_NASA_archive()
%LOAD_NASA_ARCHIVE Loads a table of confirmed exoplanet data from the NASA
%archive.
%
% OUTPUTS:
%   T[table]            Table containing confirmed exoplanet data with the 
%                       following fields:
%       kepid[double]                Kepler ID of the planet
%       kepoi_name[string]           Kepler Object of Interest (KOI) name
%       Tp[double]                   Equilibrium temp. of the planet [K]
%       Rp[double]                   Radius of the planet [Earth radii]
%       Ts[double]                   Effective temp. of the host star [K]
%       Rs[double]                   Radius of the host star [Solar radii]
%       AngSep[double]               Ang. sep. between pl. and star [rad]
%       ap[double]                   Sem-maj. ax. of planet's orbit [AU]
%
% REFERENCES:
%   Ecaut & Mamajek (2013, ApJS, 208, 9)
%   URL: http://adsabs.harvard.edu/abs/2013ApJS..208....9P
%   NASA exoplanets archive
%
% VERSION HISTORY:
%   2025-04-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load table from CSV, skipping header lines
opts = detectImportOptions('others/NASA_archive.csv', 'NumHeaderLines', 43);
data = readtable('others/NASA_archive.csv', opts);

% Load table from Mamajek
opts = detectImportOptions('others/Mamajek.csv');
mamajek = readtable('others/Mamajek.csv', opts);

% Clean data (remove rows with NaNs or invalid BCv)
valid = ~isnan(mamajek.Teff) & ~isnan(mamajek.BCv);
Teff_vals = mamajek.Teff(valid);
BCv_vals = mamajek.BCv(valid);

% Filter for CONFIRMED planets
confirmed = data(strcmp(data.koi_disposition, 'CONFIRMED'), :);

% Interpolate BCv for your stars' temperatures
T_star = confirmed.koi_steff; % or however you call it
BC = interp1(Teff_vals, BCv_vals, T_star, 'linear', 'extrap');

% Constants
T_sun = 5778;          % K
R_sun = 6.957e8;       % m
L_sun = 3.828e26;      % W
AU = 1.496e11;         % m
Mbol_sun = 4.74;       % -

% Preallocate arrays
n = height(confirmed);
distance_pc = zeros(n, 1);
angular_sep_rad = zeros(n, 1);

% For every confirmed planet...
for i = 1:n

    % Extract data based on header
    T_star = confirmed.koi_steff(i);
    R_star = confirmed.koi_srad(i); % in solar radii
    a = confirmed.koi_sma(i);       % in AU

    % Estimate luminosity ratio
    L_ratio = (R_star)^2 * (T_star / T_sun)^4;

    % Bolometric magnitude
    Mbol_star = Mbol_sun - 2.5 * log10(L_ratio);

    % Estimate bolometric correction
    M_J = Mbol_star - BC(i);

    % Apparent magnitude (use koi_jmag if available)
    m_J = confirmed.koi_jmag(i);
    if isnan(m_J)
        distance_pc(i) = NaN;
        angular_sep_rad(i) = NaN;
        continue;
    end

    % Compute distance in parsecs
    distance_pc(i) = 10^((m_J - M_J + 5) / 5);

    % Convert to AU
    distance_au = distance_pc(i) * 206265;

    % Angular separation θ = a / d [radians]
    angular_sep_rad(i) = a / distance_au;
end

% Create new output table
T = table(confirmed.kepid, confirmed.kepoi_name, ...
    confirmed.koi_teq, confirmed.koi_prad, ...
    confirmed.koi_steff, confirmed.koi_srad, ...
    angular_sep_rad, confirmed.koi_sma, ...
    'VariableNames', {'kepid', 'kepoi_name', 'Tp', ...
    'Rp', 'Ts', 'Rs', 'AngSep', 'ap'});

end