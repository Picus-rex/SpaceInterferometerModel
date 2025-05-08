function IWA = compute_IWA(unique_baselines, lambda)
%COMPUTE_IWA Compute the inner working angle using the simplified relation
%on nulling baselines.
%
% INPUTS:
%   unique_baselines[table] Table from the classification of baselines as
%                       given by classify_baselines.m.
%   lambda[1]           Reference wavelength. [m]
%
% OUTPUTS:
%   IWA[1]              Internal working angle. [rad]
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   - Uses the simplified relation based on nulling baselines.
%
% VERSION HISTORY:
%   2025-02-18 -------- 1.0
%   2025-03-18 -------- 2.0
%                     - FIX: computation of the modulation using the
%                       correct relation in the paper, using the new added
%                       input as coefficients C_i. 
%   2025-03-20 -------- 2.0.1
%                     - Moved efficiency computation to dedicated function.
%   2025-03-24 -------- 3.0
%                     - FIX: now uses the simplified relation which is more
%                       reliable than the old approach.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find B_null
nulling_baselines_idx = unique_baselines.Imaging_Flag == 0;
nulling_baselines = unique_baselines(nulling_baselines_idx, :);
B_null = max(nulling_baselines.B);

% Relation at page 9
IWA = lambda / (2 * B_null);

end