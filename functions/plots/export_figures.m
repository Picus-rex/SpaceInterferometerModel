function export_figures(export_settings)
%EXPORT_SETTINGS Save figures following specified settings in an automated
%way. 
%
% FIELDS:
%   Font_size
%   Width
%   Height
%   Name
%
% VERSION HISTORY:
%   2025-03-24 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    export_settings.embedded (1, 1) struct
    export_settings.font_size (1, 1) {mustBeNumeric} = 11
    export_settings.width (1, 1) {mustBeNumeric} = 10
    export_settings.height (1, 1) {mustBeNumeric} = 10
    export_settings.name (1, 1) string
end

if isfield(export_settings, "embedded")

    % Name can be passed externally in case of multiple subfigures
    if isfield(export_settings, "name")
        name = export_settings.name;
    else
        name = export_settings.embedded.name;
    end

    export_settings = export_settings.embedded;
    export_settings.name = name;
end

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', export_settings.font_size);

set(gcf, 'Units', 'centimeters', 'Position', [2, 2, export_settings.width, export_settings.height]);
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [export_settings.width, export_settings.height]);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, export_settings.name, '-dpdf', '-r300');


end