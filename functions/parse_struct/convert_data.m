function data = convert_data(data)
%CONVERT_DATA Convert specific entries of configuration files from human
%readable to S.I. units.

% Conversions
pc2m = 3.0857e16;

if isfield(data.environment, "target_distance")
    data.environment.target_distance = data.environment.target_distance * pc2m;
end

if isfield(data.environment, "disturbances") && ...
        isfield(data.environment.disturbances, "exozodiacal_inclination")
    data.environment.disturbances.exozodiacal_inclination = ...
        deg2rad(data.environment.disturbances.exozodiacal_inclination);
end

end