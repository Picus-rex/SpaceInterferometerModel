function h = plot_psf_map(theta_range, PSF, export_settings)
%PLOT_PSF_MAP Create the point spread function map image given the map as
%computed by compute_psf with settings for exporting.
%
% INPUTS:
%   theta_range[Nx1]    Range of coordinates to consider. [rad]
%   PSF[NxN]            PSF map resulting from function. [-]
%  optional:
%   export_setting[struct] Setting for export. 
%
% OUTPUTS:
%   h                   Handle of the figure.
%
% VERSION HISTORY:
%   2025-03-24 -------- 1.0
%
% Author: Francesco De Bortoli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    export_settings = NaN;
end

% Call styling
style_colors;

% Conversion factor
conversion_rad2mas = 1e3 * (3600 * 180) / pi;

% Convert data and generate empty arrays
theta_range = theta_range * conversion_rad2mas;

h = figure; 
hold on;
imagesc(theta_range, theta_range, PSF);
xlabel('\theta_x [mas]');
ylabel('\theta_y [mas]');
colormap(darkBlue)
colorbar();
axis xy; 
axis equal;

if isstruct(export_settings)
    export_figures("embedded", export_settings)
end

end