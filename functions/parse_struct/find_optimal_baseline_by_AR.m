function AR_opt = find_optimal_baseline_by_AR(minAR, maxAR, intensities, ...
    array, baseline, phase_shifts, theta_star, lambda, autoplot)
%FIND_OPTIMAL_BASELINE_BY_AR Find the best optimal aspect ratio to
%increment the nulling ratio.
%
% INPUTS:
%   minAR[1]            Minimum aperture ratio to consider. [-]
%   maxAR[1]            Maximum aperture ratio to consider. [-]
%   intensities[1xN]    Amplitudes associated to each aperture.
%   array[string]       Type of array (see define_array()).
%   baseline[1]         First baseline of the array.
%   phase_shifts[1xN]   Phase associated to each aperture [rad]. 
%   theta_star[1]       Angular dimension of the star. [rad]
%   lambda[1]           Wavelength of observation (default: 1e-5). [m]
%   autoplot[bool]      If true, visualise the plot (default: true).
%
% OUTPUTS:
%   AR_opt[1]            Optimal baseline to increment nulling ratio. [m]
%
% VERSION HISTORY:
%   2025-05-08 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    minAR
    maxAR
    intensities 
    array 
    baseline 
    phase_shifts 
    theta_star 
    lambda = 1e-5
    autoplot = true
end

% Define and allocate
N = length(phase_shifts);
ARs = linspace(minAR, maxAR, 1000);
ratios = zeros(length(ARs), 1);

% For every possible baseline...
for i = 1:length(ARs)
    
    % Update the array and compute the nulling ratio
    positions = define_array(array, baseline, ARs(i));
    ratios(i) = compute_nulling_ratio(N, intensities, phase_shifts, ...
                            positions, theta_star, lambda);
end

if autoplot
    figure;
    semilogy(ARs, ratios, "LineWidth", 1.5)
    xlabel("Apertures ratio")
    ylabel("Nulling ratio")
    grid minor;
end

% Extract the best point.
[~, i_opt] = min(ratios);
AR_opt = ARs(i_opt);

end