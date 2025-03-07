function data = convert_data(data)
%CONVERT_DATA Convert specific entries of configuration files from human
%readable to S.I. units. and compute some missing values.

% Conversions factors
pc2m = 3.0857e16;

% Conver fields
if isfield(data.environment, "target_distance")
    data.environment.target_distance = data.environment.target_distance * pc2m;
end

if isfield(data.instrument, "phase_shifts")
    data.instrument.phase_shifts = cell2mat(data.instrument.phase_shifts);
    data.instrument.phase_shifts = deg2rad(data.instrument.phase_shifts);
end

if isfield(data.environment, "disturbances") && ...
        isfield(data.environment.disturbances, "exozodiacal_inclination")
    data.environment.disturbances.exozodiacal_inclination = ...
        deg2rad(data.environment.disturbances.exozodiacal_inclination);
end

data.instrument.intensities = cell2mat(data.instrument.intensities);

% Compute elements
if isfield(data.environment, "target_distance") && isfield(data.environment, "star_radius")
    
    if isfield(data.environment, "stellar_angular_radius")
        warning("If target_distance and star_radius are specified, stellar_angular_radius is overwritten!")
    end

    data.environment.stellar_angular_radius = data.environment.star_radius ...
        / data.environment.target_distance; 

end

if isfield(data.environment, "target_distance") && isfield(data.environment, "exoplanet_radius")
    
    if isfield(data.environment, "exoplanet_angular_radius")
        warning("If target_distance and exoplanet_radius are specified, exoplanet_angular_radius is overwritten!")
    end

    data.environment.exoplanet_angular_radius = data.environment.exoplanet_radius ...
        / data.environment.target_distance; 

end

if isfield(data.instrument, "efficiencies")
    data.instrument.throughput = data.instrument.efficiencies.beam_combiner * data.instrument.efficiencies.optical_line;
end

% Creation of the grid for the response
[data.simulation.theta_range, data.simulation.theta_x, ...
    data.simulation.theta_y] = define_range(data.simulation.angular_extension);

% Positions (x, y) of the apertures 
data.instrument.positions = define_array(data.instrument.array, ...
    data.instrument.baseline, data.instrument.apertures_ratio);

% Plot array
if data.outputs.plot_array
    data.instrument.diameter = plot_apertures(data.instrument.positions, ...
        data.instrument.intensities, true);
else
    data.instrument.diameter = plot_apertures(data.instrument.positions, ...
        data.instrument.intensities, false);
end

end