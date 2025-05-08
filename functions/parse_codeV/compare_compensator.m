function compare_compensator(opd_uncorr, null_uncorr, opd_corr, null_corr, x, y)
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
% VERSION HISTORY:
%   2025-05-07 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute variables sizes
Ns = size(opd_corr, 2);

% Compute rms for uncorrected and corrected variables
rms_uncorr = sqrt( mean( opd_uncorr.^2, 1 ) );   % 1×Ns
rms_corr   = sqrt( mean( opd_corr.^2,   1 ) );   % 1×Ns

% PLOT 1 - Whether compensator reduces wavefront error
figure; 
boxplot([rms_uncorr; rms_corr]',...
        [zeros(1,Ns) ones(1,Ns)], ...
        'Labels',{'Uncorrected','Corrected'});
ylabel('Root mean squares of OPD');
title('Reduction in wavefront error');
grid minor;

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
scatter(log_null_uncorr_mean, log_null_corr_mean, 15, rms_uncorr, 'filled');
lims = [min(log_null_uncorr, [], "all") max(log_null_uncorr, [], "all")];
plot(lims, lims, 'k--','LineWidth',1);
xlabel('Uncorrected Nulling Ratio');
ylabel('Corrected Nulling Ratio');
colorbar; ylabel(colorbar,'Initial RMS OPD');
title('Nulling Ratio: impact of compensation');
linear2power_thicks("xy")
axis equal;
grid minor;

% Compute difference and plot it
dlog_map        = mean_log_corr - mean_log_uncorr;  

% PLOT 4
plot_value_on_image_plane(dlog_map, x, y, ...
    title="Mean difference of logarithm of Nulling Ratio", type="_1e0");

end