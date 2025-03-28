function [optical_path, x_coords, y_coords] = load_opd(filename)
%LOAD_OPD Load optical path difference data from a specified text file,
%exported following the CODE V implemented MACRO-PLUS.
%
% INPUTS:
%   filename[string]    The name of the text file containing the data.
%
% OUTPUTS:
%   optical_path[matrix]  Matrix where each row corresponds to a series and
%                         each column corresponds to a point in the series.
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
        data = sscanf(line, '%f')';
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

% Compute the OPD with respect to the center value
for i = 1:size(optical_path, 2)
    optical_path(:, i) = optical_path(:, i) - optical_path(floor(size(optical_path, 1)/2), i);
end

end