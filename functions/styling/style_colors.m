colours = [0, 47, 97;   % #002f61
          54, 84, 125;  % #36547d
          94, 117, 148; % #5e7594
          130, 146, 170; % #8292aa
          163, 174, 190; % #a3aebe
          195, 200, 208; % #c3c8d0
          226, 226, 226] / 255; % #e2e2e2

colours_contrast = [0, 47, 97;   % #002f61
          54, 84, 125;  % #36547d
          94, 117, 148; % #5e7594
          130, 146, 170; % #8292aa
          163, 174, 190; % #a3aebe
          195, 200, 208] / 255; % #c3c8d0
         

ui_colours = [
        43, 80, 170;   % Sapphire
        158, 43, 37;   % Auburn
        208, 207, 236; % Lavender (web)
        165, 230, 186; % Celadon
        255, 127, 17   % Orange (wheel)
    ] / 255;


colours = xml2matrix("others/autumn_night.xml");

% Create a custom colormap with more interpolated points
darkBlue = interp1(linspace(0, 1, size(colours,1)), colours, linspace(0, 1, 256), 'linear');

angle_darkBlue = [interp1(linspace(0, 1, size(colours, 1)), colours, linspace(0, 1, 128), 'linear'); ...
    interp1(linspace(0, 1, size(colours, 1)), colours, linspace(1, 0, 128), 'linear')];

% Create a new colormap with dark red at the zero position
% and interpolate the rest of the colors
darkBlueZero = zeros(256, 3);
darkBlueZero(1, :) = [139, 0, 0] / 255;

% Interpolate the remaining colors
darkBlueZero(2:256, :) = interp1(linspace(0, 1, size(colours, 1)), colours, linspace(1/255, 1, 255), 'linear');