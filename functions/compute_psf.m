function [PSF, theta_FWHM] = ...
    compute_psf(lambda, A, positions, phases, THETA, theta_range)
%COMPUTE_PSF Compute the Point Spread Function (PSF) using baselines
%
% INPUTS:
%   lambda[1]           Wavelength of observation. [m]
%   A[Nx1]              Vector of amplitudes for each collector.
%   positions[Nx2]      Matrix of (x, y) coordinates of the N. [m]
%   phases[Nx1]         Vector of phase shifts associated with each 
%                       aperture. [rad]
%   THETA[1x2]          (x, y) position of the point source [rad]
%   theta_range[Mx1]    Range of angular offsets to compute the PSF. [rad]
%
% OUTPUTS:
%   PSF[MxM]            Computed PSF over the given angular range.
%   theta_FWHM[1]       Angular resolution of the array [rad].
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   - Uses Bessel functionJ_0(2piB_1 theta / lambda) for convolution.
%
% VERSION HISTORY:
%   2025-02-14 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of apertures
N = length(A);

% Compute pairwise baselines and their contributions (C values) and store
% them in a vector where each row will be [|Delta_x|, |Delta_y|, C_jk]
baseline_list = []; 
for j = 1:N-1
    for k = j+1:N
        Delta_x = positions(j,1) - positions(k,1);
        Delta_y = positions(j,2) - positions(k,2);
        dphi = phases(j) - phases(k);
        C_pair = 2 * A(j) * A(k) * sin(dphi);

        % Use absolute values for grouping (for the Kronecker delta)
        baseline_list = [baseline_list; abs(Delta_x), abs(Delta_y), C_pair];
    end
end

% Group baselines with identical separations. [Delta_x_i, Delta_y_i, C_i]
tol = 1e-10;
unique_baselines = []; 

for i = 1:size(baseline_list,1)
    dx = baseline_list(i,1);
    dy = baseline_list(i,2);
    C_val = baseline_list(i,3);
    found = false;

    for j = 1:size(unique_baselines,1)
        if (abs(unique_baselines(j,1) - dx) < tol) && (abs(unique_baselines(j,2) - dy) < tol)
            unique_baselines(j,3) = unique_baselines(j,3) + C_val;
            found = true;
            break;
        end
    end

    if ~found
        unique_baselines = [unique_baselines; dx, dy, C_val];
    end
end

% Number of distinct baselines
num_baselines = size(unique_baselines,1);

% Compute the length of each baseline: B_i = sqrt((Delta_x_i)^2 + (Delta_y_i)^2)
B = sqrt(unique_baselines(:,1).^2 + unique_baselines(:,2).^2);

% Compute normalization factor: N_c and its rms as for the formula
THETA_x = THETA(1);
THETA_y = THETA(2);
N_c = 0;

for i = 1:num_baselines
    dx = unique_baselines(i,1);
    dy = unique_baselines(i,2);
    C_i = unique_baselines(i,3);
    N_c = N_c + C_i * sin( (2*pi/lambda) * (dx*THETA_x + dy*THETA_y) );
end

% For a single value, rms(N_c) is the absolute value; avoid diving by 0.
rms_Nc = abs(N_c);
if rms_Nc < 1e-6
    rms_Nc = 1;
end

% Compute angular resolution as the FWHM
num = 0;
den = 0;
for m = 1:num_baselines
    num = num + unique_baselines(m, 3)^2 * B(m);
    den = den + unique_baselines(m, 3)^2;
end
B_av = num / den;
theta_FWHM = 0.48 * lambda / B_av;

% Create the meshgrid over theta offsets
[theta_x_grid, theta_y_grid] = meshgrid(theta_range, theta_range);

% Compute the dirty map (PSF)
PSF = zeros(size(theta_x_grid));

for m = 1:num_baselines
    for n = 1:num_baselines

        C_m = unique_baselines(m,3);
        C_n = unique_baselines(n,3);
        B_n = B(n); 
        
        Delta_x_m = unique_baselines(m,1);
        Delta_y_m = unique_baselines(m,2);
        Delta_x_n = unique_baselines(n,1);
        Delta_y_n = unique_baselines(n,2);
        
        % Compute beta_x and beta_y for this pair (m,n)
        beta_x = (Delta_x_m*Delta_x_n*THETA_x + Delta_y_m*Delta_x_n*THETA_y + ...
                  Delta_y_m*Delta_y_n*THETA_x - Delta_x_m*Delta_y_n*THETA_y) / (B_n^2);
        beta_y = (Delta_x_m*Delta_y_n*THETA_x + Delta_y_m*Delta_y_n*THETA_y - ...
                  Delta_y_m*Delta_x_n*THETA_x + Delta_x_m*Delta_x_n*THETA_y) / (B_n^2);
              
        arg1 = (2*pi*B_n/lambda) * sqrt((theta_x_grid - beta_x).^2 + (theta_y_grid - beta_y).^2);
        arg2 = (2*pi*B_n/lambda) * sqrt((theta_x_grid + beta_x).^2 + (theta_y_grid + beta_y).^2);
        
        % Evaluate the Bessel functions J0 (besselj with order 0)
        J0_1 = besselj(0, arg1);
        J0_2 = besselj(0, arg2);
        
        PSF = PSF + C_m * C_n * (J0_1 - J0_2);
    end
end

% Normalize the PSF by the factor 2*rms(N_c)
PSF = PSF / (2 * rms_Nc);

end
