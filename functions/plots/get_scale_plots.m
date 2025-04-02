function [scale, scale_tag] = get_scale_plots(type_label)
%GET_SCALE_PLOTS Extract the scale and corresponding label from a type 
% label for the use in integrated plots.
%
% INPUTS:
%   type_label[string]  Label containing the scale information in format
%                       'prefix_scale', e.g., 'm_1e-6' or 'N_1e9'.
%
% OUTPUTS:
%   scale[double]       Numeric inverted scale extracted from the type 
%                       label: if the scale is 1e-6, then scale is 1e6 so
%                       that it can be directly multiplied for scaling
%                       plots.
%   scale_tag[string]   Corresponding scale label, e.g., 'µm' for 1e-6.
%
% VERSION HISTORY:
%   2025-04-01 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the prefix and scale using regular expression
prefix_match = regexp(type_label, '^([a-zA-Z]+)_', 'tokens');
scale_match = regexp(type_label, '_(\d+e[+-]?\d+)', 'tokens');

% Assign prefix and scale
if ~isempty(prefix_match)
    prefix = prefix_match{1}{1};
else
    prefix = '';
end

if ~isempty(scale_match)
    scale = str2double(scale_match{1}{1});
else
    scale = 1;
end

% Define the scale to tag mapping for metric prefixes
metric_prefixes = containers.Map({1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e6, 1e9, 1e12, 1e15, 1e18, 1e21, 1e24}, ...
                                  {'y', 'z', 'a', 'f', 'p', 'n', 'µ', 'm', 'c', 'd', '', 'da', 'h', 'k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y'});

% Find the corresponding scale tag
if metric_prefixes.isKey(scale)
    scale_tag = [metric_prefixes(scale), prefix];
else
    % Handle cases where the scale is not in the predefined map
    scale_tag = sprintf('%s × 10^{%d}', prefix, log10(scale));
end

% Reverse the scale for direct multiplication
scale = 1 / scale;


end
