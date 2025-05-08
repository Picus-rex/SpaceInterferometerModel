function [ang_sep, stellar_ang_radius] = select_random_exoplanet(filename)


% Conversion factor arcsec to radians
arcsec_to_rad = pi / (180 * 3600);

data = readtable(filename, 'Delimiter', '\t');

% Select a random planet and extract values (see documentation)
planet = data(randi(height(data)), :);
ang_sep = planet.ang_sep * arcsec_to_rad;           % [rad]
Rs = planet.Rs * 6.955e8;                           % [m]
Ds = planet.dist * 3.086e16;                        % [m]

% Stellar angular radius using trigonometry
stellar_ang_radius = asin(Rs / Ds);

end
