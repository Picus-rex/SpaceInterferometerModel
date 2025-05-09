function plot_interferogram_results(RFn, RFs, R, ratio, eta_mod, x, y, theta, maps_to_compute)
%PLOT_INTERFEROGRAM_RESULTS Create study plots for the interferogram based
%on perturbed simulations. 
%
% INPUTS:
%   RFn[1 x Nt]         Vector of nominal response function
%   RFs[Ns x Nt]        Matrix of perturbed response functions
%   R[Np x Nt x Ns]     3D matrix of the response map over apertures
%   ratio[Np x Ns]      Nulling ratio for every point for every simulation
%   eta_mod[Np x Ns]    Modulation eff. for every point for every sim. 
%   x[Np x 1]           X coordinates on the pupil screen. [m]
%   y[Np x 1]           Y coordinates on the pupil screen. [m]
%   theta[Nt x 1]       Vector of angular coordinates [rad] 
%   maps_to_compute[1xNm] Indices of the maps to include in the R 3D matrix
%
% VERSION HISTORY:
%   2025-05-09 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For all results requested with maps_to_compute, plots them.
for i = 1:size(R, 3)
    label = sprintf("Intensity response at observation angle of %.2f mas", rad2mas(theta(maps_to_compute(i))));
    plot_value_on_image_plane(R(:, i, 1), x(:, 1), y(:, 1), title=label, type="_1e0");
end

% Group results for plotting the reponse function
RFp = [RFn; RFs];
plot_response_function_theta(theta, RFp, "Normalize", true);

% Plot also the nulling ratio (which depends only on the simulation)
plot_value_on_image_plane(ratio(:, 1), x(:, 1), y(:, 1), title="Nulling ratio", type="log");

% Plot also the modulation efficiency
plot_value_on_image_plane(eta_mod(:, 1), x(:, 1), y(:, 1), title="Modulation efficiency", type="_1e0");

% Mean, Median, Max across Ns simulations
log_mean = mean(ratio, 2);
eff_mean = mean(eta_mod, 2);
log_median = median(ratio, 2);
log_max = max(ratio, [], 2); 

% eff_mean = mean(eff, 2);
plot_value_on_image_plane(eff_mean, x(:, 1), y(:, 1), ...
    title="Mean Modulation Efficiency", type="1e0");

plot_value_on_image_plane(log_mean, x(:, 1), y(:, 1), ...
    title="Mean Log Nulling Ratio", type="log");

plot_value_on_image_plane(log_median, x(:, 1), y(:, 1), ...
    title="Median Log Nulling Ratio", type="log");

plot_value_on_image_plane(log_max, x(:, 1), y(:, 1), ...
    title="Worst-case Log Nulling Ratio", type="log");

% Threshold for "good" nulling in log10 scale
threshold = 1e-4;
Ns = size(ratio, 2);
Np = size(ratio, 1);

% Count good cases for each pupil point
good_counts = sum(ratio < threshold, 2);
fraction_good = good_counts / Ns;

% Plot how frequently each position is "good"
plot_value_on_image_plane(fraction_good, x(:, 1), y(:, 1), ...
    title="Fraction of Good Nulling Simulations", type="linear");

% What % of the pupil is reliably good
good_pupil_fraction = sum(fraction_good > 0.9) / Np;  % e.g., 90% of sims are good
fprintf("Fraction of pupil with consistently good nulling: %.2f%%\n", 100 * good_pupil_fraction);

end
