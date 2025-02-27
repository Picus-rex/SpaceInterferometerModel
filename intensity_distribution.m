%% INTENSITY DISTRIBUTION
% Following pages 54-... of Viseur, L., the intensity of a 4 interferometer
% is derived, assuming a perfect signal). The procedure is generalisable to
% any N-apertures interferometer.

clc; clear; close;

%% Data

N = 4;                              % Aperture number

% For each aperture (3rd dimension), define the phase vector for the ideal 
% case (unperturbed). This should be optimised using the
% compute_optimal_splitting function. 
phases(:, :, 1) = pi * ones(1000);
phases(:, :, 2) = 0 * ones(1000);
phases(:, :, 3) = -pi * ones(1000);
phases(:, :, 4) = 0 * ones(1000);

%% Intesity

% The intensity is half the square of the amplitude
I = 0;

for i = 1:N
    for j = 1:N
        
        I = I + cos(phases(:, :, i) - phases(:, :, j));

    end
end

I = 0.5 * I;

%% Representation

figure;
imagesc(I);
xlabel('x');
ylabel('y');
title('Intensity Function for 4-Aperture Interferometer');
colorbar;
axis xy; axis equal;
