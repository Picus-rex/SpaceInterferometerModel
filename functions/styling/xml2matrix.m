function colours = xml2matrix(filename)
%XML2MATRIX Import a xml file and convert the colours to a matrix that can
%be used by colorbar.

% Read the XML file and get all color nodes
xmlDoc = xmlread(filename);
allColors = xmlDoc.getElementsByTagName('color');

% Initialize an empty matrix
colours = zeros(allColors.getLength(), 3);

% For each color node, convert hex to RGB
for k = 0:allColors.getLength()-1

    hexColor = char(allColors.item(k).getFirstChild.getData);
    rgb = hexToRgb(hexColor);
    colours(k+1, :) = rgb / 255; % Normalize to [0, 1]

end

end

function rgb = hexToRgb(hexColor)
    % Convert a hex color string to an RGB array
    hexColor = hexColor(2:end); % Remove the '#' character
    rgb = sscanf(hexColor, '%2x%2x%2x', [1 3])'; % Convert to RGB
end
