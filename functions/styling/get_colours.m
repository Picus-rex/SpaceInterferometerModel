function cols = get_colours(n, scale)
%GET_COLOURS Returns n colours following the defined colour scale in the
%project for specific plots.
arguments
    n 
    scale = "contrast"
end

% Define the original colours
style_colors;

if ~strcmp(scale, "contrast")
    colours_contrast = colours;
end

% Define original indices
numColours = size(colours_contrast, 1);
originalIndices = linspace(1, numColours, numColours);
targetIndices = linspace(1, numColours, n);

% Interpolate each RGB channel separately
cols = interp1(originalIndices, colours_contrast, targetIndices, 'linear');

end
