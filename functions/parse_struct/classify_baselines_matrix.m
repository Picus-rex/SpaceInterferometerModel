function [C_i, unique_baselines] = classify_baselines_matrix(amplitudes, positions, phases, tol)
%CLASSIFY_BASELINES_MATRIX This is a vectorised version of
%classify_baselines and compute_baselines_amplitude_factors for fast
%computations.
%
% INPUTS:
%   amplitudes[Nx1]     Amplitudes of each apertures. [-]
%   positions[Nx2]      Matrix of (x, y) coordinates of the N apertures.[m]
%   phases[MxN]         Matrix of phase shifts associated with each 
%                       aperture; each row is a point. [rad]
%   tol[1]              If given, tolerance for similarity (default: 1e-3).
%
% OUTPUTS:
%   C_i[N1xM]           Matrix of amplitude factor for each baseline (N1)
%                       and point (M).
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
%   - This is a combination of classify_baselines and 
%     compute_baselines_amplitude_factors that works with matrices for
%     speed, generated with LLMs and verified. 
%
% VERSION HISTORY:
%   2025-05-06 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    tol = 1e-3; 
end

[~, N] = size(phases);
amplitudes = amplitudes(:);
assert(length(amplitudes)==N, 'amp must be N×1 matching phases columns');

% reference phase‑row for grouping
phi0 = phases(1,:);

% all j<k pairs in same order as nested loops j=1:N, k=1:N
[J,K] = meshgrid(1:N,1:N);
mask   = J < K;
j_idx  = J(mask);  % P×1 ordering: (1,2),(1,3)...(1,N),(2,3)...
k_idx  = K(mask);
P      = numel(j_idx);

% geometry
dX = positions(j_idx,1) - positions(k_idx,1);
dY = positions(j_idx,2) - positions(k_idx,2);

% reference dphi0 for grouping, make column
dphi0 = (phi0(j_idx) - phi0(k_idx)).';   % P×1

% quantize by tol
qX   = round(abs(dX)/tol)*tol;           % P×1
qY   = round(abs(dY)/tol)*tol;           % P×1
qPhi = round(abs(dphi0)/tol)*tol;        % P×1

% grouping into G groups, preserve first-seen "stable" order
[uniqKeys, ~, ic] = unique([qX, qY, qPhi], 'rows', 'stable');
G = size(uniqKeys,1);

% output table of unique baselines
B = sqrt(uniqKeys(:,1).^2 + uniqKeys(:,2).^2);
unique_baselines = table(uniqKeys(:,1), uniqKeys(:,2), B, ...
                         'VariableNames',{'Delta_x','Delta_y','B'});

% weights w_p = a_j * a_k
w = amplitudes(j_idx) .* amplitudes(k_idx);

% M×P matrix of sin Δφ for each phase‑row
Dphi = phases(:,j_idx) - phases(:,k_idx);   % M×P
S    = sin(Dphi);                           % M×P

% weight and sum into G groups via sparse
W    = S .* (w.');                         % M×P
Sg   = sparse(1:P, ic, 1, P, G);            % P×G
Cmat = W * Sg;                              % M×G

% double per Lay eq. 8, and transpose to G×M
C_i = 2 * Cmat.';  

end
