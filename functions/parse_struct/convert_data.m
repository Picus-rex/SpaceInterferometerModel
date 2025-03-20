function data = convert_data(data)
%CONVERT_DATA Convert specific entries of configuration files from human
%readable to S.I. units. and compute some missing values in preparation to
%the analysis. 

fields_inst = fieldnames(data.instrument);
fields_envi = fieldnames(data.environment);
fields_simu = fieldnames(data.simulation);

% Conversion factors
pc2m = 3.0857e16;

for i = 1:length(fields_inst)
    switch fields_inst{i}
        case "target_distance"
            data.environment.target_distance = data.environment.target_distance * pc2m;
        case "phase_shifts"
            data.instrument.phase_shifts = cell2mat(data.instrument.phase_shifts);
            data.instrument.phase_shifts = deg2rad(data.instrument.phase_shifts);
            [data.instrument.baselines, data.instrument.unique_baselines] = classify_baselines(data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, false);
        case "efficiencies"
            data.instrument.throughput = data.instrument.efficiencies.beam_combiner * data.instrument.efficiencies.optical_line;
        case "intensities"
            data.instrument.intensities = cell2mat(data.instrument.intensities);
        case "array"
            data.instrument.positions = define_array(data.instrument.array, data.instrument.baseline, data.instrument.apertures_ratio);
    end

end


% Additional computations
[data.instrument.diameter, data.instrument.surfaces] = ...
    plot_apertures(data.instrument.positions, ...
    data.instrument.intensities, ...
    data.instrument.efficiencies.optical_line, ...
    data.instrument.efficiencies.beam_combiner, data.outputs.plot_array);



for i = 1:length(fields_envi)
    switch fields_envi{i}
        case "target_distance"
            data.environment.target_distance = data.environment.target_distance * pc2m;

            if isfield(data.environment, "star_radius")
                if isfield(data.environment, "stellar_angular_radius")
                    warning("If target_distance and star_radius are specified, stellar_angular_radius is overwritten!")
                end
                data.environment.stellar_angular_radius = data.environment.star_radius / data.environment.target_distance; 
            end

            if isfield(data.environment, "target_distance") && isfield(data.environment, "exoplanet_radius")
                if isfield(data.environment, "exoplanet_angular_radius")
                    warning("If target_distance and exoplanet_radius are specified, exoplanet_angular_radius is overwritten!")
                end
                data.environment.exoplanet_angular_radius = data.environment.exoplanet_radius / data.environment.target_distance; 
            end

        case "disturbances"
            if isfield(data.environment.disturbances, "exozodiacal_inclination")
                data.environment.disturbances.exozodiacal_inclination = deg2rad(data.environment.disturbances.exozodiacal_inclination);
            end
        case "exoplanet_position"
            data.environment.exoplanet_position = cell2mat(data.environment.exoplanet_position);
    end
end



for i = 1:length(fields_simu)
    switch fields_simu{i}
        case "angular_extension"
            [data.simulation.theta_range, data.simulation.theta_x, data.simulation.theta_y] = define_range(data.simulation.angular_extension);
    end
end


end