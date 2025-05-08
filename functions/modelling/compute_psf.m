function [PSF, theta_FWHM, C_i] = ...
    compute_psf(lambda, A, positions, phases, THETA, theta_range)
%COMPUTE_PSF Compute the Point Spread Function (PSF) using baselines
%
% INPUTS:
%   lambda[1]           Wavelength of observation. [m]
%   A[Nx1]              Vector of amplitudes for each collector.
%   positions[Nx2]      Matrix of (x, y) coordinates of the N. [m]
%   phases[Nx1]         Vector of phase shifts associated with each 
%                       aperture. [rad]
%   THETA[1x2]          (x, y) position of the point source. [rad]
%   theta_range[Mx1]    Range of angular offsets to compute the PSF. [rad]
%
% OUTPUTS:
%   PSF[MxM]            Computed PSF over the given angular range.
%   theta_FWHM[1]       Angular resolution of the array. [rad]
%   C_i[N1x1]           Baseline amplitude factors. [-]
%
% where N1 is the number of unique baselines.
%
% REFERENCES:
%   Lay OP. Imaging properties of rotating nulling interferometers. 2005;
%
% NOTES:
%   - Uses Bessel functionJ_0(2piB_1 theta / lambda) for convolution.
%
% VERSION HISTORY:
%   2025-02-14 -------- 1.0
%   2025-03-18 -------- 2.0
%                     - Moved baseline computations to different function.
%                     - Correction in the computation of amplitude factors.
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract baselines info from the respective function
[~, unique_baselines] = classify_baselines(A, positions, phases, false);
C_i = unique_baselines.C_i;
num_baselines = size(unique_baselines, 1);
B = unique_baselines.B;

% Compute normalization factor: N_c and its rms as for the formula
THETA_x = THETA(1);
THETA_y = THETA(2);
N_c = 0;

for i = 1:num_baselines
    dx = unique_baselines(i,:).Delta_x;
    dy = unique_baselines(i,:).Delta_y;
    C = C_i(i);
    N_c = N_c + C * sin( (2*pi/lambda) * (dx*THETA_x + dy*THETA_y) );
end

% For a single value, rms(N_c) is the absolute value; avoid diving by 0.
rms_Nc = abs(N_c);
if rms_Nc < 1e-6
    rms_Nc = 1;
end

% Compute angular resolution as the FWHM following relations 18, 19 from
% the tables.
num = 0;
den = 0;
for m = 1:num_baselines
    num = num + C_i(m)^2 * B(m);
    den = den + C_i(m)^2;
end
B_av = num / den;
theta_FWHM = 0.48 * lambda / B_av;

% Create the meshgrid over theta offsets
[theta_x_grid, theta_y_grid] = meshgrid(theta_range, theta_range);

% Compute the dirty map (PSF)
PSF = zeros(size(theta_x_grid));

for m = 1:num_baselines
    for n = 1:num_baselines

        C_m = C_i(m);
        C_n = C_i(n);
        B_n = B(n); 
        
        Delta_x_m = unique_baselines(m,:).Delta_x;
        Delta_y_m = unique_baselines(m,:).Delta_y;
        Delta_x_n = unique_baselines(n,:).Delta_x;
        Delta_y_n = unique_baselines(n,:).Delta_y;
        
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