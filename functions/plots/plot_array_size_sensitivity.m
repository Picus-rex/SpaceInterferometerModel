function plot_array_size_sensitivity(intensities, array, baseline, ...
    apertures_ratio, phase_shifts, theta_star, lambda)
%PLOT_ARRAY_SIZE_SENSITIVITY Visualise how dimensions affects the resulting
%metrics of nulling ratio and inner working angle in the detection of the
%system. 
%
% INPUTS:
%   intensities[1xN]    Amplitudes associated to each aperture.
%   array[string]       Type of array (see define_array()).
%   baseline[1]         First baseline of the array.
%   apertures_ratio[1]  Ratio between maximum and minimum size.
%   phase_shifts[1xN]   Phase associated to each aperture [rad]. 
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation (default: 1e-5). [m]
%
% VERSION HISTORY:
%   2025-05-08 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define and allocate
Bs = linspace(100, 600, 1000);
ARs = linspace(0.4, 6, 1000);
IWAs = zeros(length(Bs), 1);

% PART 1: EDIT THE BASELINE

% Plot nulling ratio
find_optimal_baseline_by_B(Bs(1), Bs(end), intensities, ...
    array, apertures_ratio, phase_shifts, theta_star, lambda);

% For every possible baseline...
for i = 1:length(Bs)
    
    % Update the array and compute the nulling ratio
    positions = define_array(array, Bs(i), apertures_ratio);
    [~, unique_baselines] = classify_baselines(intensities, ...
        positions, phase_shifts, false);
    IWAs(i) = compute_IWA(unique_baselines, lambda);

end

figure;
plot(Bs, rad2mas(IWAs), "LineWidth", 1.5);
xlabel("Maximum baseline")
ylabel("Inner working angle [mas]")
grid minor;

% PART 2: EDIT THE APERTURE RATIO

% Plot nulling ratio
find_optimal_baseline_by_AR(ARs(1), ARs(end), intensities, ...
    array, baseline, phase_shifts, theta_star, lambda);

% For every possible baseline...
for i = 1:length(ARs)
    
    % Update the array and compute the nulling ratio
    positions = define_array(array, baseline, ARs(i));
    [~, unique_baselines] = classify_baselines(intensities, ...
        positions, phase_shifts, false);
    IWAs(i) = compute_IWA(unique_baselines, lambda);
end

figure;
plot(ARs, rad2mas(IWAs), "LineWidth", 1.5);
xlabel("Aperture ratio")
ylabel("Inner working angle [mas]")
grid minor;

end