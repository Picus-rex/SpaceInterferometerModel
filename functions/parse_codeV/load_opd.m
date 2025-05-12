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
%   2025-04-10 -------- 1.2
%                     - Completely rewritten to be significantly faster.
%   2025-05-12 -------- 1.3
%                     - Allow split importing for larger files.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read all lines from the file
lines = read_split_file(filename);

% Identify series start markers
markerIdx = find(startsWith(lines, '=='));
Ns_est = numel(markerIdx);

% Store all valid series temporarily
temp_series = cell(1, Ns_est);
Np_list = zeros(1, Ns_est);

% First pass: collect parsed data and count points
for s = 1:Ns_est
    if s < Ns_est
        dataBlock = lines(markerIdx(s)+1 : markerIdx(s+1)-1);
    else
        dataBlock = lines(markerIdx(s)+1 : end);
    end

    % Remove empty lines
    dataBlock = dataBlock(strlength(dataBlock) > 0);

    % Parse block
    combinedStr = join(dataBlock, newline);
    parsedData = sscanf(combinedStr, '%f,%f,%f');
    parsedData = reshape(parsedData, 3, []);
    
    temp_series{s} = parsedData';
    Np_list(s) = size(parsedData, 1);
end

% Find most common length (mode) across all series
Np = size(temp_series{1}, 1);
Ns = numel(temp_series);

if Ns < Ns_est
    warning('%d series discarded due to inconsistent point count.', Ns_est - Ns);
end

% Preallocate output matrices
optical_path = zeros(Np, Ns);
x_coords     = zeros(Np, Ns);
y_coords     = zeros(Np, Ns);

% Fill the matrices
for j = 1:Ns
    parsedData = temp_series{j};
    optical_path(:, j) = parsedData(:, 1) * 1e-2;
    x_coords(:, j)     = parsedData(:, 2) * 1e-2;
    y_coords(:, j)     = parsedData(:, 3) * 1e-2;
end

% Compute OPD
center_index = ceil(Np / 2);
OPD = optical_path - optical_path(center_index, :);
end




function lines = read_split_file(filename)
    % Extract folder, base name, and extension
    [folder, base, ext] = fileparts(filename);
    basefile = fullfile(folder, [base, ext]);

    % Initialize array for all lines
    lines = [];

    % Check if split parts exist (e.g., file.1.txt, file.2.txt, ...)
    part_index = 1;
    has_parts = false;
    while true
        part_file = fullfile(folder, sprintf('%s.%d%s', base, part_index, ext));
        if isfile(part_file)
            has_parts = true;
            part_lines = readlines(part_file);
            lines = [lines; part_lines];  %#ok<AGROW> <-- disable warning
            part_index = part_index + 1;
        else
            break;
        end
    end

    % If no parts were found, read the original file
    if ~has_parts
        if isfile(basefile)
            lines = readlines(basefile);
        else
            error('File "%s" not found.', basefile);
        end
    end
end
