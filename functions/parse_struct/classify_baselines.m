function [baselines, unique_baselines] = classify_baselines(amplitudes, ...
    positions, phases, export)
%CLASSIFY_BASELINES Classify interferometric baselines as imaging or 
% nulling as detailed in the refence.
%
% INPUTS:
%   amplitudes[Nx1]     Amplitudes of each apertures. [-]
%   positions[Nx2]      Matrix of (x, y) coordinates of the N apertures.[m]
%   phases[Nx1]         Vector of phase shifts associated with each 
%                       aperture. [rad]
%   export[bool]        If true, write to command window the results.
%
% OUTPUTS:
%   baselines[table]    Table containing the j, k indices, the positions
%                       shifts and the corresponding baselines, phases
%                       differences and the sine of the value along with
%                       the imaging flag:
%                        - imaging_flag = 1 if imaging baseline
%                        - imaging_flag = 0 if nulling baseline
%                        - imaging_flag = -1 otherwise (mixed situation)
%   unique_baselines[table]Table containing the the positions shifts and 
%                       the corresponding baselines, phases differences and 
%                       the sine of the value along with the imaging flag:
%                        - imaging_flag = 1 if imaging baseline
%                        - imaging_flag = 0 if nulling baseline
%                        - imaging_flag = -1 otherwise (mixed situation)
%                       The field contributing includes the couples
%                       contributing to that baselines as strings,
%                       separated by a dash "-" and the amplitudes factors
%                       C_i used to compute the PSF and the modulation
%                       efficiency from Lay.
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   - A baseline is unique if both baselines and phases differences are
%     different.
%
% VERSION HISTORY:
%   2025-02-14 -------- 1.0
%   2025-02-28 -------- 1.1
%                     - Corrected defintion based on the used convention
%                       (this does not correspond anymore to the reference)
%   2025-03-18 -------- 2.0
%                     - Now the function outputs a table both for all
%                       baselines and for the unique ones following the
%                       concept in Lay.
%                     - Revert to the correct definition, awaiting for
%                       confirmation.
%                     - The function now computes all the coefficients that
%                       are relative to the apertures.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(phases);

% Find all baselines except self term using combinations
baselines = zeros(N^2-N, 8);
unique_baselines = zeros(N^2-N, 6);
tol = 1e-3;

% Now consider all the possible combinations to extract all the baselines,
% computing data for baselines and wether the baseline is imaging (sin phi
% ~= 0) or nulling (cos phi ~= 0)
i = 1;
for j = 1:N
    for k = 1:N
        if j ~= k
            Delta_x = positions(j,1) - positions(k,1);
            Delta_y = positions(j,2) - positions(k,2);
            B = sqrt(Delta_x^2 + Delta_y^2);
            dphi = phases(j) - phases(k);
            sphi = sin(dphi);
            cphi = cos(dphi);

            if abs(sphi) > tol
                flag = 1;
            elseif abs(cphi) > tol
                flag = 0;
            else
                flag = -1;
            end

            baselines(i, :) = [j, k, Delta_x, Delta_y, B, dphi, sphi, flag];

            i = i + 1;
        end
    end
end

% Now uniform the description by extracting only the unique baselines by
% looking at the the entire array again.
k = 1;
for i = 1:size(baselines, 1)
    
    found = false;

    for j = 1:size(unique_baselines, 1)
        if abs(baselines(i, 5) - unique_baselines(j, 3)) < tol && ...
                abs(abs(baselines(i, 6)) - abs(unique_baselines(j, 4))) < tol 
            found = true;
            break 
        end
    end

    if ~found
        unique_baselines(k, :) = baselines(i, 3:end);
        k = k + 1;
    end
end

% Strip empty lines and make the baselines positive as for the procedure in
% Lay from eq. 8 (a negative value could never contribute to the expression)
unique_baselines = unique_baselines(1:k-1, :);
unique_baselines(:, 1:2) = abs(unique_baselines(:, 1:2));

% Add the contributing baselines to the array:
contributing = strings(size(unique_baselines, 1), size(baselines, 1));
maxks = zeros(size(baselines, 1), 1);
for i = 1:size(unique_baselines, 1)
    k = 1;
    for j = 1:size(baselines, 1)
        if abs(baselines(j, 5) - unique_baselines(i, 3)) < tol && ...
                 abs(abs(baselines(j, 6)) - abs(unique_baselines(i, 4))) < tol 
            contributing(i, k) = string(baselines(j, 1)) + "-" + string(baselines(j, 2));
            k = k + 1;
        end
    end
    maxks(i) = k - 1;
end

% Strip empty column
contributing = contributing(:, 1:max(maxks));

% Allocate space for baseline amplitude factors C_i; for every unique
% baseline follows equation 8 from Lay. At the end, avoiding double
% contributions, the sign is doubled.
C_i  = zeros(size(unique_baselines, 1), 1);
for i = 1:size(unique_baselines, 1)
    
    repeating_couples = strings(1, size(contributing, 2));
    n = 1;

    % Extract contributing baselines from the seen apertures before
    for l = 1:size(contributing, 2)
        if ~strcmp(contributing(i, l), "")
            
            % Transform back to numbers
            apertures = split(contributing(i, l), '-');
            j = str2double(apertures(1));
            k = str2double(apertures(2));
            
            found = false;
            
            for m = 1:length(repeating_couples)
                if strcmp(string(k)+"-"+string(j), repeating_couples(m)) || strcmp(string(j)+"-"+string(k), repeating_couples(m))
                    found = true;
                    break
                end
            end

            % Debug info
            if export && ~found
                fprintf("[%.0f] Contribution from apertures %.0f, %.0f\n", i, j, k)
            end
            
            if ~found
                % Compute the amplitude factor
                C_i(i) = C_i(i) + amplitudes(j) * amplitudes(k) * sin(phases(j) - phases(k));
                repeating_couples(n) = contributing(i, l);
                n = n + 1;
            end
   
        end
    end
end

% Convert to tables
baselines = array2table(baselines, "VariableNames", {'j', 'k', 'Delta_x', ...
    'Delta_y', 'B', 'Delta_phi', 'sin_Delta_phi', 'Imaging_Flag'});
unique_baselines = array2table(unique_baselines, 'VariableNames', {'Delta_x', ...
    'Delta_y', 'B', 'Delta_phi', 'sin_Delta_phi', 'Imaging_Flag'});
unique_baselines.contributing = contributing;
unique_baselines.C_i = C_i * 2;

if export
    disp('Baseline Classification (Imaging = 1, Nulling = 0):');
    disp(unique_baselines);
end

end
