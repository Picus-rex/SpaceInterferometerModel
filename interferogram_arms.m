%% INTERFEROGRAM ARMS
% See how the response depends on the phases on the x, y found rays from
% CODE V response.

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') 

% Loading of elements from configurations files
data = ReadYaml('config/linear_array.yml');
data = convert_data(data);

[optical_path, OPD, x_coords, y_coords] = load_opd("code_v/TestLens_21804_convex.txt");
[optical_path_n, OPD_n, x_coords_n, y_coords_n] = load_opd("code_v/TestLens_21804.txt");

delta_phi = opd2phase(OPD(:, 1), data.instrument.wavelength);
delta_phi_n = opd2phase(OPD_n(:, 1), data.instrument.wavelength);

% Computing of optimal splitting for analysis
[data.simulation.U, data.instrument.combination, ...
    data.instrument.phase_shifts] = compute_optimal_splitting(...
    data.instrument.apertures, data.instrument.baseline, ...
    data.instrument.wavelength, data.instrument.positions(:, 1), ...
    data.instrument.positions(:, 2), ...
    data.environment.stellar_angular_radius, false);

theta = linspace(-700, 700, 1000) * ((pi/180)/3600 / 1e3);
R = zeros(length(x_coords), length(theta));
ty = 0;
positions = data.instrument.positions;
lambda = data.instrument.wavelength;
nominal_phases = data.instrument.phase_shifts;
combination = data.instrument.combination;
A = data.instrument.intensities;
N = size(positions, 1);

for i = 1:length(theta)
    
    t = theta(i);
    E = zeros(length(x_coords), 1);

    for k = 1:N
    
        xk = positions(k,1);
        yk = positions(k,2);
        
        % Perturb first arm, leave others unchanged
        if k == 1
            phase_shifts = nominal_phases(k) + delta_phi;
        else
            phase_shifts = nominal_phases(k) + delta_phi_n;
        end
        
        % Phase term
        phase_k = 2 * pi * (xk * t + yk * ty) / lambda + phase_shifts;
        
        % Add contribution from each aperture
        E = E + combination(k) * A(k) * exp(1i * phase_k);
    
    end
    
    R(:, i) = abs(E).^2;

end

for i = [1, 10, 100, 200, 250]
    h = plot_value_on_image_plane(R(:, i), x_coords(:, 1), y_coords(:, 1), title="Intensity response", type="linear_1e0");
end

RF = sum(R, 1);

figure;
plot(theta / ((pi/180)/3600 / 1e3), RF, "LineWidth", 1.5);
xlabel("Observation angle [mas]")
ylabel("Response function")