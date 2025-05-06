function C_i = compute_baselines_amplitude_factor(num_unique_baselines, contributing, amplitudes, phases)
%COMPUTE_BASELINES_AMPLITUDE_FACTOR Compute the C_i parameter associated to
%each unique_baseline.
%
% INPUTS:
%   num_unique_baselines[1] Number of unique baselines. [-]
%   contributing[N1xN2] Matrix of strings that contains, for each unique
%                       baselines, the indices of the contributing
%                       apertures (e.g. "1-2", "3-4"). [strings]
%   amplitudes[Nx1]     Amplitudes of each apertures. [-]
%   phases[NxN3]        Vector of phase shifts associated with each 
%                       aperture. [rad]
%
% OUTPUTS:
%   C_i[N1x1]           Vector of amplitude factor for each baseline.
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   - A baseline is unique if both baselines and phases differences are
%     different.
%
% VERSION HISTORY:
%   2025-04-30 -------- 1.0
%   2025-05-06 -------- 1.1
%                     - Replaced strings with cells to increase speed of
%                       the executions.  
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate space for baseline amplitude factors C_i; for every unique
% baseline follows equation 8 from Lay. At the end, avoiding double
% contributions, the sign is doubled.
C_i  = zeros(num_unique_baselines, size(phases, 1));
for i = 1:num_unique_baselines
    
    repeating_couples = cell(1, size(contributing, 2));
    n = 1;

    % Extract contributing baselines from the seen apertures before
    for l = 1:size(contributing, 2)
        if ~isempty(contributing{i, l})
            
            % Transform back to numbers
            apertures = cell2mat(contributing{i, l});
            j = apertures(1);
            k = apertures(2);
            
            found = false;
            
            for m = 1:length(repeating_couples)
                if isequal({k, j}, repeating_couples{m}) || isequal({j, k}, repeating_couples{m})
                    found = true;
                    break
                end
            end
            
            if ~found
                % Compute the amplitude factor
                C_i(i, :) = C_i(i, :) + (amplitudes(j) * amplitudes(k) * sin(phases(:, j) - phases(:, k)))';
                repeating_couples{n} = contributing{i, l};
                n = n + 1;
            end
   
        end
    end
end

end