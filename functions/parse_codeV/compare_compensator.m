function compare_compensator(opd_uncorr, null_uncorr, opd_corr, null_corr, x, y, export_setup)
%COMPARE_COMPENSATOR Provide some figures to compare the results between
%uncorrected and corrected simulations.
%
% INPUTS:
%   opd_uncorr[NpxNs]   OPD of uncorrected simulations from CODE V [m].
%   null_uncorr[NpxNs]  Nulling ratio of uncorrected simulations [-].
%   opd_corr[NpxNs]     OPD of corrected simulations from CODE V [m].
%   null_corr[NpxNs]    Nulling ratio of corrected simulations [-].
%   x[Npx1]             X coordinates at the pupil screen [m].
%   y[Npx1]             Y coordinates at the pupil screen [m].
%
% OPTIONAL INPUTS:
%   export_setup[struct]Struct with data to save exporting.
%
% VERSION HISTORY:
%   2025-05-07 -------- 1.0
%   2025-05-12 -------- 1.1
%                     - Added export setup options
%   2025-05-13 -------- 1.2
%                     - Boxplot replaced by boxchart for customisations.
%                     - Added colour scale at correlation plot and improved
%                       it with larger dots.
%                     - Normalised and inverted the map value.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7
    export_setup = [];
end

% Compute variables sizes
Ns = size(opd_corr, 2);
style_colors;

% Compute rms for uncorrected and corrected variables
rms_uncorr = sqrt( mean( opd_uncorr.^2, 1 ) );   % 1×Ns
rms_corr   = sqrt( mean( opd_corr.^2,   1 ) );   % 1×Ns

% PLOT 1 - Whether compensator reduces wavefront error
% Combine data into a single vector
ydata = [rms_corr, rms_uncorr] * 1e6;

% Create corresponding group labels
xgroupdata = [repmat("Corrected", 1, Ns), repmat("Uncorrected", 1, Ns)];

% Create the box chart
figure;
boxchart(categorical(xgroupdata), ydata, "BoxFaceColor", colours(5, :), "MarkerColor", colours(end, :));
ylabel('Root mean squares of OPD [µm]');
title('Reduction in wavefront error');
grid minor;

if iscell(export_setup) && isstruct(export_setup{1})
    glob_name = export_setup{1}.name;
    export_setup{1}.name = glob_name + "_boxplot";
    export_figures("embedded", export_setup{1})
end

% Convert nulling ratios to log values, then compute the mean for all
% values to flatten values
log_null_uncorr = log10(null_uncorr);    
log_null_corr   = log10(null_corr);      
dlog_null = log_null_corr - log_null_uncorr;     

log_null_uncorr_mean = mean(log_null_uncorr, 1);   % 1xNs
log_null_corr_mean   = mean(log_null_corr,   1); 

mean_log_uncorr = mean(log_null_uncorr(:, 1), 2);  % Np×1
mean_log_corr   = mean(log_null_corr(:, 1),   2);

% PLOT 2 - Overall trend: is compensation beneficial?
figure; 
histogram(dlog_null,30);
linear2power_thicks("x")
xlabel('Difference in nulling ratio');
ylabel('Count');
title('Improvement in log nulling ratio');
grid minor;

% PLOT 3
figure; hold on;
scatter(log_null_uncorr_mean, log_null_corr_mean, 30, rms_uncorr*1e6, 'filled');
colormap(darkBlue)
lims = [min(log_null_uncorr, [], "all") max(log_null_uncorr, [], "all")];
plot(lims, lims, 'k--','LineWidth',1);
xlabel('Uncorrected Nulling Ratio');
ylabel('Corrected Nulling Ratio');
colorbar; ylabel(colorbar,'Uncorrected OPD (RMS) [µm]');
title('Nulling Ratio: impact of compensation');
linear2power_thicks("xy")
axis equal;
grid minor;

if iscell(export_setup) && isstruct(export_setup{2})
    export_setup{2}.name = glob_name + "_correlation";
    export_figures("embedded", export_setup{2})
end

% Compute difference and plot it
dlog_map        = mean_log_corr - mean_log_uncorr; 
dlog_map = -dlog_map / max(abs(dlog_map), [], "all");

% PLOT 4
if iscell(export_setup) && isstruct(export_setup{3})
    export_setup{3}.name = glob_name + "_map";
    plot_value_on_image_plane(dlog_map, x, y, ...
        title="Mean difference of logarithm of Nulling Ratio", type="_1e0", embedded=export_setup{3});
else
    plot_value_on_image_plane(dlog_map, x, y, ...
    title="Mean difference of logarithm of Nulling Ratio", type="_1e0");
end


end