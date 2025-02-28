function imaging_baselines = classify_baselines(positions, phases, export)
%CLASSIFY_BASELINES Classify interferometric baselines as imaging or 
% nulling as detailed in the refence.
%
% INPUTS:
%   positions[Nx2]      Matrix of (x, y) coordinates of the N apertures.[m]
%   phases[Nx1]         Vector of phase shifts associated with each 
%                       aperture. [rad]
%   export[bool]        If true, write to command window the results.
%
% OUTPUTS:
%   imaging_baselines[Mx3]  Matrix where each row corresponds to a baseline 
%                       and contains:[aperture_i, aperture_j, imaging_flag]
%                       where:
%                        - aperture_i, aperture_j are indices of the 
%                          baseline pair
%                        - imaging_flag = 1 if imaging baseline
%                        - imaging_flag = 0 if nulling baseline
%                        - imaging_flag = NaN otherwise (mixed situation)
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   - A baseline may contribute to both imaging and nulling.
%
% VERSION HISTORY:
%   2025-02-14 -------- 1.0
%   2025-02-28 -------- 1.1
%                     - Corrected defintion based on the used convention
%                       (this does not correspond anymore to the reference)
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(positions, 1); 
num_baselines = nchoosek(N, 2);         % Number of unique baseline pairs

imaging_baselines = zeros(num_baselines, 3); 

index = 1;
for i = 1:N-1
    for j = i+1:N
        
        % Compute phase difference between aperture pairs
        delta_phi = abs(abs(phases(i)) - abs(phases(j)));
        
        % Determine if it is an imaging baseline
        if abs(delta_phi - pi) < 1e-6
            is_imaging = false;
        elseif delta_phi < 1e-6
            is_imaging = true;
        else
            is_imaging = NaN;
        end
        
        % Store result
        imaging_baselines(index, :) = [i, j, is_imaging];
        index = index + 1;
    end
end

if export
    disp('Baseline Classification (Imaging = 1, Nulling = 0):');
    disp(imaging_baselines);
end

end
