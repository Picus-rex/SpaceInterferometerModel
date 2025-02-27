function colours = styling(export, w, h, name, cbar)
    
if nargin == 0
    export = false;
end

if nargin <= 4
    cbar = true;
end

% Define the given colors as RGB values (normalized to [0,1])
colours = [0, 47, 97;   % #002f61
          54, 84, 125;  % #36547d
          94, 117, 148; % #5e7594
          130, 146, 170; % #8292aa
          163, 174, 190; % #a3aebe
          195, 200, 208; % #c3c8d0
          226, 226, 226] / 255; % #e2e2e2

% Define positions (normalized to [0,1])
x = linspace(0, 1, size(colours,1));

% Create a custom colormap with more interpolated points
numPoints = 256; % Define resolution of colormap
xq = linspace(0, 1, numPoints);
customColormap = interp1(x, colours, xq, 'linear');

% Display the colormap
if cbar
    colormap(customColormap);
    colorbar;
else
    colorbar("off");
end

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

if export
    set(gcf, 'Units', 'centimeters', 'Position', [2, 2, w, h]);
    set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [w, h]);
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, name, '-dpdf', '-r300');
end

end
