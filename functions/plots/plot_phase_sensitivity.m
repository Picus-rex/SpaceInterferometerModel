function h = plot_phase_sensitivity(surfaces, intensities, positions, ...
    nominal_phases, eta_opt, lambda, theta_star, export_setup)
%PLOT_PHASE_SENSITIVITY Visualise how the modulation efficiency and the
%nulling ratio change with respect to the phase of the first element.
%
% INPUTS:
%   surfaces[Nx1]       Areas associated to each apertures [m^2].
%   intensities[Nx1]    Amplitudes associated to each aperture.
%   positions[Nx2]      Matrix of (x, y) positions of the apertures. [m]
%   nominal_phases[Nx1] Vector of phase shifts for each aperture 
%   eta_opt[1]          Optical line efficiency.
%   lambda[1]           Wavelength of observation. [m]
%   theta_star[1]       Angular dimension of the star. [rad]
%
% OPTIONAL INPUTS:
%   export_setup[struct]Struct with data to save exporting.
%
% OUTPUTS:
%   h[figure]           Handle of the figure. 
%
% VERSION HISTORY:
%   2025-05-06 -------- 1.0
%   2025-05-12 -------- 1.1
%                     - Added export setup options. 
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of apertures
N = length(nominal_phases);

% Initialise phases to be changed and variables...
phase = linspace(0.1, 360, 1000);
eta_mod = zeros(length(phase), 1);
nulling = zeros(length(phase), 1);

% For each phase...
for i = 1:length(phase)

    % Compute elements
    phase_shifts = [deg2rad(phase(i)), nominal_phases(2:end)];
    [~, unique_baselines] = classify_baselines(intensities, positions, phase_shifts, false);
    C_i = unique_baselines.C_i;
    
    % Store results
    eta_mod(i) = compute_modulation_efficiency(eta_opt, surfaces, C_i);
    [nulling(i), ~] = compute_nulling_ratio(N, intensities, ...
          phase_shifts, positions, theta_star, lambda);
end

% Plotting
figure;
plot(phase, eta_mod, "LineWidth", 1.5);
ylabel("Modulation efficiency")
xlim([0, 360])
xticks(0:30:360)
xlabel("Phase of the perturbed aperture [deg]")
grid minor

if isstruct(export_setup)
    glob_name = export_setup.name;
    export_setup.name = glob_name + "_modulation";
    export_figures("embedded", export_setup)
end

figure;
plot(phase, nulling, "LineWidth", 1.5);
ylabel("Nulling ratio")
set(gca, "YScale", "log")
xlim([0, 360])
xticks(0:30:360)
xlabel("Phase of the perturbed aperture [deg]")
grid minor

if isstruct(export_setup)
    export_setup.name = glob_name + "_nul_ratio";
    export_figures("embedded", export_setup)
end


end