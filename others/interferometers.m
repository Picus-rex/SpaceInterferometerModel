% Specify the file names
exoplanetFile = 'TestPlanetPopulation.txt';
observatoryFile = 'space_observatories.csv';

arcsec2rad = (pi/180)/3600;
h = 6.626e-34;                  % Planck constant [J s]
c = 3e8;                        % Speed of light [m/s]
k = 1.381e-23;                  % Boltzmann constant [J/K]
R_earth = 6.371e6;              % Earth radius [m]
R_sun   = 6.957e8;              % Solar radius [m]
Baseline = 10;                  % Baseline [m]
Diameter = 2.5;                 % Diameter [m]

B = @(T_val, lam) (2*h*c^2 ./ lam^5) ./ (exp(h*c./(lam*k*T_val)) - 1);


% Read the exoplanet catalogue
opts = detectImportOptions(exoplanetFile, 'Delimiter', {'\t',' '});
T = readtable(exoplanetFile, opts);

% Read the observatory data
opts_obs = detectImportOptions(observatoryFile, 'Delimiter', ';', 'VariableNamingRule','preserve');
observatories = readtable(observatoryFile, opts_obs);

% Initialize variables for results
numUniverses = 10;
detectablePlanets = zeros(numUniverses, height(observatories));

% Loop through each observatory
for obsIdx = 1:height(observatories)
    % Extract observatory data
    minWav = observatories{obsIdx, 'Min Wav [microm]'} * 1e-6; % Convert to meters
    maxWav = observatories{obsIdx, 'Max Wav [microm]'} * 1e-6; % Convert to meters
    baseline = observatories{obsIdx, 'Baseline [m]'};
    diameter = observatories{obsIdx, 'Aperture [m]'};
    flux_ratio = observatories{obsIdx, 'Nulling needed'} * 1e-4;
    
    % Calculate mean wavelength
    lambda = (minWav + maxWav) / 2;
    
    % Calculate IWA and OWA (in radians)
    IWA = lambda / baseline; % IWA in radians
    OWA = lambda / diameter;  % OWA in radians
    
    % Loop through each universe
    for universe = 0:numUniverses-1
        % Filter exoplanets for the current universe
        rows = T.nMC == universe;
        Tsub = T(rows, :);
        
        % Calculate the flux contrast
        contrast_1 = ((Tsub.rp * R_earth).^2 .* B(Tsub.Tp, lambda1)) ./ ((Tsub.Rs * R_sun).^2 .* B(Tsub.Ts, lambda1));
        
        % Check for detectable exoplanets
        detectableCount = sum(Tsub.ang_sep * arcsec2rad < OWA & Tsub.ang_sep * arcsec2rad > IWA & contrast_1 > flux_ratio);
        detectablePlanets(universe + 1, obsIdx) = detectableCount; % Store the count
    end
end

% Average the results over the universes
averageDetectablePlanets = mean(detectablePlanets, 1);

col = styling;

% Create a bar plot
figure;
b = bar(averageDetectablePlanets);

b.FaceColor = 'flat'; % Allow for individual bar colors
b.CData = col(1:3, :); % Assign the colors from the col matrix

set(gca, 'XTickLabel', observatories.Name, 'XTickLabelRotation', 45);
xlabel('Observatories');
ylabel('Detectable Exoplanets');
grid on;
colorbar('off')
styling(true, 10, 7, "exports/space_observatories", false);
