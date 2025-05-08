function linear2power_thicks(ax)
%LINEAR2POWER_THICKS Convert thicks of the specified axes from linear (-1,
%1, 3, ...) to the same value in power (10^{-1}, 10^1, 10^3, ...) for plots
%that must be done in logarithmic format but can be read normally.
%
% INPUTS:
%   ax[string]          Either "x", "y" or "xy"
%
% VERSION HISTORY:
%   2025-05-07 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(ax, "x") || strcmp(ax, "xy")
    existingTicks = xticks;
    newTickLabels = arrayfun(@(t) sprintf('10^{%d}', t), existingTicks, 'UniformOutput', false);
    xticks(existingTicks); 
    xticklabels(newTickLabels); 
end

if strcmp(ax, "y") || strcmp(ax, "xy")
    existingTicks = yticks;
    newTickLabels = arrayfun(@(t) sprintf('10^{%d}', t), existingTicks, 'UniformOutput', false);
    yticks(existingTicks); 
    yticklabels(newTickLabels); 
end

end