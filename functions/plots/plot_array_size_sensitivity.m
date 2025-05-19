function [optimal_B, optimal_AR, optimal_NR, optimal_IWA] = ...
    plot_array_size_sensitivity(intensities, array, baseline, ...
    apertures_ratio, phase_shifts, theta_star, lambda, export_setup, weight, autoplot)
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
% OPTIONAL INPUTS:
%   export_setup[struct]Struct with data to save exporting.
%   weight[1]           Value ∈ [0, 1] to prioritise the nulling ratio over 
%                       the IWA. When it's 0.5, both have the same
%                       importance. (default: 0.5)
%   autoplot[bool]      Create plots. (default: true)
%
% OUTPUTS
%   optimal_B[1]        Baseline that provides the best metrics in the 
%                       given range. [m]
%   optimal_AR[1]       Apertures ratio that provides the best metrics in
%                       the given range.
%   optimal_NR[1x2]     Nulling ratio resulting using the optimal baseline
%                       (first element) and the optimal AR (second elem).
%   optimal_IWA[1x2]    IWA resulting using the optimal baseline (first
%                       element) and the optimal AR (second elem). [rad]
%
% VERSION HISTORY:
%   2025-05-08 -------- 1.0
%   2025-05-12 -------- 1.1
%                     - Added export setup.
%   2025-05-15 -------- 1.2
%                     - Removed references to deprecated functions. 
%                     - Added optimiser and additional point plot.
%                     - Added outputs of optimisations.
%   2025-05-16 -------- 1.3
%                     - Added autoplot input argument. 
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    intensities 
    array 
    baseline 
    apertures_ratio 
    phase_shifts 
    theta_star 
    lambda = 1e-5
    export_setup = []
    weight = 0.5
    autoplot = true
end

% Define and allocate
Bs = linspace(100, 600, 1000);
ARs = linspace(0.4, 6, 1000);
IWAs = zeros(length(Bs), 1);
ratios = zeros(length(Bs), 1);
N = length(phase_shifts);

% PART 1: EDIT THE BASELINE

[optimal_B, optimal_NR(1), optimal_IWA(1)] = optimise_baseline(array, intensities, ...
     phase_shifts, theta_star, lambda, apertures_ratio, [Bs(1), Bs(end)], weight);

% For every possible baseline...
for i = 1:length(Bs)
    
    % Update the array and compute the nulling ratio
    positions = define_array(array, Bs(i), apertures_ratio);
    ratios(i) = compute_nulling_ratio(N, intensities, phase_shifts, ...
                            positions, theta_star, lambda);
    [~, unique_baselines] = classify_baselines(intensities, ...
        positions, phase_shifts, false);
    IWAs(i) = compute_IWA(unique_baselines, lambda);

end

if autoplot
    style_colors;
    
    figure; hold on;
    plot(Bs, ratios, "LineWidth", 1.5, "Color", ui_colours(6, :));
    plot(optimal_B, optimal_NR(1), '.', "MarkerSize", 10);
    set(gca, "YScale", "log");
    xlabel("Maximum baseline [m]")
    ylabel("Nulling ratio")
    grid minor;
    if isstruct(export_setup)
        export_setup.name = glob_name + "_nulling_baseline";
        export_figures("embedded", export_setup)
    end
    
    figure; hold on;
    plot(Bs, rad2mas(IWAs), "LineWidth", 1.5, "Color", ui_colours(2, :));
    plot(optimal_B, rad2mas(optimal_IWA(1)), '.', "MarkerSize", 10);
    xlabel("Maximum baseline [m]")
    ylabel("Inner working angle [mas]")
    grid minor;
    if isstruct(export_setup)
        export_setup.name = glob_name + "_iwa_baseline";
        export_figures("embedded", export_setup)
    end
end

% PART 2: EDIT THE APERTURE RATIO

[optimal_AR, optimal_NR(2), optimal_IWA(2)] = optimise_apertures_ratio(array, intensities, ...
     phase_shifts, theta_star, lambda, baseline, [ARs(1), ARs(end)], weight);

% For every possible baseline...
for i = 1:length(ARs)
    
    % Update the array and compute the nulling ratio
    positions = define_array(array, baseline, ARs(i));
    ratios(i) = compute_nulling_ratio(N, intensities, phase_shifts, ...
                            positions, theta_star, lambda);
    [~, unique_baselines] = classify_baselines(intensities, ...
        positions, phase_shifts, false);
    IWAs(i) = compute_IWA(unique_baselines, lambda);
end

if autoplot
    figure; hold on;
    plot(ARs, ratios, "LineWidth", 1.5, "Color", ui_colours(8, :));
    plot(optimal_AR, optimal_NR(2), '.', "MarkerSize", 10);
    set(gca, "YScale", "log");
    xlabel("Aperture ratio")
    ylabel("Nulling ratio")
    grid minor;
    if isstruct(export_setup)
        export_setup.name = glob_name + "_nulling_aspectratio";
        export_figures("embedded", export_setup)
    end
    
    figure; hold on;
    plot(ARs, rad2mas(IWAs), "LineWidth", 1.5, "Color", ui_colours(4, :));
    plot(optimal_AR, rad2mas(optimal_IWA(2)), '.', "MarkerSize", 10);
    xlabel("Aperture ratio")
    ylabel("Inner working angle [mas]")
    grid minor;
    if isstruct(export_setup)
        export_setup.name = glob_name + "_iwa_aspectratio";
        export_figures("embedded", export_setup)
    end
end

end