function data = convert_data(data)
%CONVERT_DATA Convert specific entries of configuration files from human
%readable to S.I. units. and compute some missing values in preparation to
%the analysis. 

% If the input is the path to the configuration file, convert it to a
% struct calling the appropriate function.
if ischar(data) || isstring(data)
    data = ReadYaml(data);
end

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
            if iscell(data.instrument.phase_shifts)
                data.instrument.phase_shifts = cell2mat(data.instrument.phase_shifts);
                data.instrument.phase_shifts = deg2rad(data.instrument.phase_shifts);
            end
            [data.instrument.baselines, data.instrument.unique_baselines] = classify_baselines(data.instrument.intensities, data.instrument.positions, data.instrument.phase_shifts, false);
        case "efficiencies"
            data.instrument.throughput = data.instrument.efficiencies.beam_combiner * data.instrument.efficiencies.optical_line;
        case "intensities" 
            if iscell(data.instrument.intensities)
                data.instrument.intensities = cell2mat(data.instrument.intensities);
            end
        case "combination"
            if iscell(data.instrument.combination)
                data.instrument.combination = cell2mat(data.instrument.combination);
            end
        case "array"
            if ~strcmp(data.instrument.array, "Custom")
                data.instrument.positions = define_array(data.instrument.array, data.instrument.baseline, data.instrument.apertures_ratio);
            else
                if ~isfield(data.instrument, "positions")
                    error("When the array is Custom, positions field must be manually defined in configuration!")
                end
                if iscell(data.instrument.positions)
                    data.instrument.positions = cell2mat(data.instrument.positions);
                end
            end
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
            if iscell(data.environment.exoplanet_position) 
                data.environment.exoplanet_position = cell2mat(data.environment.exoplanet_position);
            end
    end
end



for i = 1:length(fields_simu)
    switch fields_simu{i}
        case "phases"

            if strcmp(data.simulation.phases, "optimal")
                [data.simulation.U, data.instrument.combination, ...
                data.instrument.phase_shifts] = compute_optimal_splitting(...
                    data.instrument.apertures, data.instrument.baseline, ...
                    data.instrument.wavelength, data.instrument.positions(:, 1), ...
                    data.instrument.positions(:, 2), ...
                    data.environment.stellar_angular_radius, false);
            end

        case "angular_extension"
            [data.simulation.theta_range, data.simulation.theta_x, data.simulation.theta_y] = define_range(data.simulation.angular_extension);

        case "code_v_opd_file"
            
            if isfield(data.simulation.code_v_opd_file, "compensator")
                [op, opd, x, y] = load_opd(data.simulation.code_v_opd_file.compensator);

                if isfield(data.instrument, "wavelength")
                    pha = opd2phase(opd, data.instrument.wavelength);
                end
                
                data.simulation.code_v.nominal.op    =  op(:, 2); 
                data.simulation.code_v.nominal.opd   = opd(:, 2); 
                data.simulation.code_v.nominal.x     =   x(:, 2); 
                data.simulation.code_v.nominal.y     =   y(:, 2); 
                data.simulation.code_v.nominal.phase = pha(:, 2); 

                op(:, 1:2) = [];
                opd(:, 1:2) = [];
                x(:, 1:2) = [];
                y(:, 1:2) = [];
                pha(:, 1:2) = [];

                data.simulation.code_v.perturbed.op    =  op(:, 1:2:end); 
                data.simulation.code_v.perturbed.opd   = opd(:, 1:2:end); 
                data.simulation.code_v.perturbed.x     =   x(:, 1:2:end); 
                data.simulation.code_v.perturbed.y     =   y(:, 1:2:end); 
                data.simulation.code_v.perturbed.phase = pha(:, 1:2:end);

                data.simulation.code_v.corrected.op    =  op(:, 2:2:end); 
                data.simulation.code_v.corrected.opd   = opd(:, 2:2:end); 
                data.simulation.code_v.corrected.x     =   x(:, 2:2:end); 
                data.simulation.code_v.corrected.y     =   y(:, 2:2:end); 
                data.simulation.code_v.corrected.phase = pha(:, 2:2:end);

            else

                if isfield(data.simulation.code_v_opd_file, "nominal")
                    [data.simulation.code_v.nominal.op, ...
                     data.simulation.code_v.nominal.opd, ...
                     data.simulation.code_v.nominal.x, ...
                     data.simulation.code_v.nominal.y] = load_opd(data.simulation.code_v_opd_file.nominal);
    
                    if isfield(data.instrument, "wavelength")
                        data.simulation.code_v.nominal.phase = opd2phase(data.simulation.code_v.nominal.opd, data.instrument.wavelength);
                    end
                end
    
                if isfield(data.simulation.code_v_opd_file, "perturbed")
                    [data.simulation.code_v.perturbed.op, ...
                     data.simulation.code_v.perturbed.opd, ...
                     data.simulation.code_v.perturbed.x, ...
                     data.simulation.code_v.perturbed.y] = load_opd(data.simulation.code_v_opd_file.perturbed);
    
                    if isfield(data.instrument, "wavelength")
                        data.simulation.code_v.perturbed.phase = opd2phase(data.simulation.code_v.perturbed.opd, data.instrument.wavelength);
                    end
                end
            end
    end
end


end