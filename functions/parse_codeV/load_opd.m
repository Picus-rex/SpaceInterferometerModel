function [optical_path, OPD, x_coords, y_coords] = load_opd(filename)
%LOAD_OPD Load optical path difference data from a specified text file,
%exported following the CODE V implemented MACRO-PLUS.
%
% INPUTS:
%   filename[string]    The name of the text file containing the data.
%
% OUTPUTS:
%   optical_path[matrix]  Matrix where each row corresponds to a series and
%                         each column corresponds to a point in the series
%                         of the optical path.
%   OPD[matrix]           Matrix where each row corresponds to a series and
%                         each column corresponds to a point in the series
%                         of the OPD.
%   x_coords[matrix]      Matrix of x coordinates corresponding to the
%                         optical path.
%   y_coords[matrix]      Matrix of y coordinates corresponding to the
%                         optical path.
%
% NOTES:
%   - The file should contain series of data points, each introduced by
%     '==' followed by an iteration number.
%   - The first column is the optical path, and the remaining columns are
%     the x and y coordinates, respectively.
%   - The center value for each series should be 0, 0.
%
% VERSION HISTORY:
%   2025-03-28 -------- 1.0
%   2025-03-31 -------- 1.1
%                     - Conversion of the outputs to S.I. units.
%                     - Small changes in the structure of the txt files
%                       that reflect on the way the input is splitted into
%                       a matrix. 
%                     - Output of OPD as well as the optical path.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize cell arrays to store data for each series
optical_path = [];
x_coords = [];
y_coords = [];

% Open the file for reading
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open the file: %s', filename);
end

% Read the file line by line
current_series = [];
i = 0;
while ~feof(fid)
    line = fgetl(fid);
    if startsWith(line, '==')
        i = i + 1;
        % If a new series starts, store the previous series data
        if ~isempty(current_series)
            optical_path(:, i-1) = current_series(:, 1);
            x_coords(:, i-1) = current_series(:, 2);
            y_coords(:, i-1) = current_series(:, 3);
        end
        % Start a new series
        current_series = [];
    elseif ~isempty(line) && ~startsWith(line, '==')
        % Read the data points
        data = strsplit(line, ',');
        data = str2double(data); 
        if length(data) == 3
            current_series(end+1, :) = data;
        end
    end
end

% Store the last series data
if ~isempty(current_series)
    optical_path(:, i) = current_series(:, 1);
    x_coords(:, i) = current_series(:, 2);
    y_coords(:, i) = current_series(:, 3);
end

% Close the file
fclose(fid);

% Conversion to S.I. units
optical_path = optical_path * 1e-2;
x_coords = x_coords * 1e-2;
y_coords = y_coords * 1e-2;

% Compute the OPD with respect to the center value. Preallocatr and compute
OPD = optical_path;
for i = 1:size(optical_path, 2)
    OPD(:, i) = optical_path(:, i) - optical_path(ceil(size(optical_path, 1)/2), i);
end

end