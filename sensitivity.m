%% SENSITIVITY ANALYSIS
% Based on
% 1. Lay OP. Systematic errors in nulling interferometers. 2004; 

clc; clear; close all;
addpath(genpath("."))
set(0, 'DefaultFigureWindowStyle', 'docked') % Change to NORMAL to export

%%

N = 4;                          % Number of apertures
lambda = 10e-6;                 % Wavelength [m]
delta_lambda = 0.5e-6;
B = 30;                         % Baseline [m]
theta_range = linspace(-0.5e-7, 0.5e-7, 2000); % Angular grid [rad]
[theta_x, theta_y] = meshgrid(theta_range, theta_range);

% Aperture positions (x,y) [m]
positions = [-B, 0; -B/3, 0; B/3, 0; B, 0];
diameter = 4;

% Ideal amplitudes and phases
A = 0.56 * ones(N, 1);
phi = [0, pi/2, pi, 3/2*pi];

Fstar = 7.5e5 / diameter;       % Stellar flux [photons/s/m^2]
theta_star = 3e-9;              % Stellar angular radius [rad]

Fplanet = 0.0775;               % Planet flux [photons/s/m^2]
theta_planet_x = 2.27e-7;       % Example planet angular separation [rad]
theta_planet_y = 0;             % Example planet angular separation [rad]

FEZ = 22.08;                            % EZ flux
theta_ez = 1.496e11 / (15 * 3.086e16);  % 

BLZ = 7.1e12;                   % LZ brightness [photons/s/m^2/sr]
DeltaOmega = 8.75e-12;          % Effective solid angle of each beam

% Perturbations
delta_a   = 0.001 * randn(N,1);      
delta_phi = 0.001 * randn(N,1);      
delta_x   = 0.001 * randn(N,1);     
delta_y   = 0.001 * randn(N,1);    

%% Star Contributions

B_star = zeros(N, N);
dB_dx = zeros(N, N);
dB_dy = zeros(N, N);

for j = 1:N
    for k = 1:N
        xjk = positions(j,1) - positions(k,1);
        yjk = positions(j,2) - positions(k,2);
        b_jk = sqrt(xjk^2 + yjk^2);
        arg = 2*pi * b_jk * theta_star / lambda;
        
        if j ~= k
            B_star(j,k) = 2 * Fstar * besselj(1, arg) / arg;
            dB_dx(j,k) = -2 * xjk * Fstar * besselj(2, arg) / b_jk^2;
            dB_dy(j,k) = -2 * yjk * Fstar * besselj(2, arg) / b_jk^2;
        else
            B_star(j,k) = 2 * Fstar * 0.5;
            dB_dx(j,k) = -2 * xjk * Fstar * 0;
            dB_dy(j,k) = -2 * yjk * Fstar * 0;
        end
        
    end
end

[C_A, C_phi, C_x, C_y, C_AA, C_Aphi, C_phiphi] = ...
    compute_coefficients(A, phi, B_star, dB_dx, dB_dy, diameter);

N_star_ideal = 0;
for j = 1:N
    for k = 1:N
        N_star_ideal = N_star_ideal + ...
            A(j) * A(k) * cos(phi(j)-phi(k)) * B_star(j,k);
    end
end

% Resulting perturbations
deltaN_star_1 = sum(C_A .* delta_a) + sum(C_phi .* delta_phi) + ...
    sum(C_x .* delta_x) + sum(C_y .* delta_y);
deltaN_star_2 = 0;
for j = 1:N
    for k = 1:N
        deltaN_star_2 = deltaN_star_2 + C_AA(j,k)*delta_a(j)*delta_a(k) ...
            + C_Aphi(j,k)*delta_a(j)*delta_phi(k) + ...
            C_phiphi(j,k)*delta_phi(j)*delta_phi(k);
    end
end

deltaN_star = deltaN_star_1 + deltaN_star_2;

%% Exozodiacal Contributions

B_EZ = zeros(N, N);
dB_dx = zeros(N, N);
dB_dy = zeros(N, N);

for j = 1:N
    for k = 1:N
        xjk = positions(j,1) - positions(k,1);
        yjk = positions(j,2) - positions(k,2);
        b_jk = sqrt(xjk^2 + yjk^2);
        arg = 2*pi * b_jk * theta_ez / lambda;
        
        if j ~= k
            B_EZ(j,k) = 2 * FEZ * besselj(1, arg) / arg;
            dB_dx(j,k) = -2 * xjk * FEZ * besselj(2, arg) / b_jk^2;
            dB_dy(j,k) = -2 * yjk * FEZ * besselj(2, arg) / b_jk^2;
        else
            B_EZ(j,k) = 2 * FEZ * 0.5;
            dB_dx(j,k) = -2 * xjk * FEZ * 0;
            dB_dy(j,k) = -2 * yjk * FEZ * 0;
        end
 
    end
end

[C_A_EZ, C_phi_EZ, C_x_EZ, C_y_EZ, C_AA_EZ, C_Aphi_EZ, C_phiphi_EZ] = ...
    compute_coefficients(A, phi, B_EZ, dB_dx, dB_dy, diameter);

N_EZ_ideal = 0;
for j = 1:N
    for k = 1:N
        N_EZ_ideal = N_EZ_ideal + A(j)*A(k)*cos(phi(j)-phi(k))*B_EZ(j,k);
    end
end 

deltaN_EZ_1 = sum(C_A_EZ .* delta_a) + sum(C_phi_EZ .* delta_phi) + ...
    sum(C_x_EZ .* delta_x) + sum(C_y_EZ .* delta_y);
deltaN_EZ_2 = 0;
for j = 1:N
    for k = 1:N
        deltaN_EZ_2 = deltaN_EZ_2 + C_AA_EZ(j,k)*delta_a(j)*delta_a(k) ...
            + C_Aphi_EZ(j,k)*delta_a(j)*delta_phi(k) + ...
            C_phiphi_EZ(j,k)*delta_phi(j)*delta_phi(k);
    end
end

deltaN_EZ = deltaN_EZ_1 + deltaN_EZ_2;

%% Local Zodiacal Contributions

N_LZ_ideal = BLZ * sum(A.^2 .* DeltaOmega);
C_A_LZ = 2 * BLZ * (A.^2 .* DeltaOmega);  
deltaN_LZ = sum(C_A_LZ .* delta_a);

%% Planet Contributions

N_planet_ideal = 0;
for j = 1:N
    for k = 1:N
        phase_term = (2*pi/lambda) * ...
            ((positions(j,1)-positions(k,1)) * theta_planet_x + ...
            (positions(j,2)-positions(k,2)) * theta_planet_y);

        N_planet_ideal = N_planet_ideal + A(j) * A(k) * ...
            (cos(phi(j)-phi(k)) * cos(phase_term) - ...
            sin(phi(j)-phi(k)) * sin(phase_term));
    end
end
N_planet_ideal = Fplanet * N_planet_ideal;

%% Total

N0 = N_star_ideal + N_EZ_ideal + N_LZ_ideal + N_planet_ideal;
deltaN_total = deltaN_star + deltaN_EZ + deltaN_LZ;

% Fractional (relative) perturbation
epsilon = deltaN_total / N0;

%% Comparison

[T_ideal, ~] = compute_response_function(lambda, N, positions, A, phi, [1, 1, 1, 1], theta_x, theta_y);

T_perturbed_norm = T_ideal * (1 + epsilon);

% Star region
x = theta_star * [-1, 1, 1, -1]; 
y = [-10, -10, 10, 10]; 

figure;
subplot(1,3,1);
imagesc(theta_range, theta_range, T_ideal);
colorbar;
xlabel('\theta_x (rad)'); ylabel('\theta_y (rad)');
title('Ideal Response Function');
axis xy; axis equal;

subplot(1,3,2);
imagesc(theta_range, theta_range, T_perturbed_norm);
colorbar;
xlabel('\theta_x (rad)'); ylabel('\theta_y (rad)');
title('Perturbed Response Function');
axis xy; axis equal;

subplot(1,3,3); hold on;
ylim([0, max(interp1(theta_range, T_perturbed_norm(1000, :), 1.5*theta_star*[-1, 1]))])
xlim(1.5*theta_star*[-1, 1]);
h = fill(x, y, 'b', 'FaceAlpha', 0.5);
plot(theta_range, T_ideal(1000, :), "LineWidth", 1.5);
plot(theta_range, T_perturbed_norm(1000, :), "LineWidth", 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C_A, C_phi, C_x, C_y, C_AA, C_Aphi, C_phiphi] = ...
    compute_coefficients(A, phi, B, dB_dx, dB_dy, diameter)

N = length(A);

% 1st order
C_A = zeros(N,1);
C_phi = zeros(N,1);
C_x = zeros(N,1);
C_y = zeros(N,1);

for j = 1:N

    sumA = 0; 
    sumPhi = 0; 
    sumX = 0; 
    sumY = 0;
    
    for k = 1:N
        sumA = sumA + A(k) * cos(phi(j)-phi(k)) * B(j,k);
        if j ~= k 
            sumPhi = sumPhi + A(k) * sin(phi(j)-phi(k)) * B(j,k);
        end
        sumX = sumX + A(j)*A(k)*cos(phi(j)-phi(k)) * dB_dx(j,k);
        sumY = sumY + A(j)*A(k)*cos(phi(j)-phi(k)) * dB_dy(j,k);
    end

    C_A(j) = 2 * A(j) * sumA;
    C_phi(j) = -2 * A(j) * sumPhi;
    C_x(j) = 2 * sumX;
    C_y(j) = 2 * sumY;
end

% 2nd order
C_AA = zeros(N,N);
C_Aphi = zeros(N,N);
C_phiphi = zeros(N,N);

for j = 1:N
    for k = 1:N
        
        C_AA(j,k) = A(j) * A(k) * cos(phi(j)-phi(k)) * B(j,k);

        if j ~= k
            
            C_Aphi(j,k) = 2 * A(j) * A(k) * sin(phi(j)-phi(k)) * B(j,k);
            C_phiphi(j,k) = A(j) * A (k) * cos(phi(j)-phi(k)) * B(j,k);

        else

            tempAphi = 0;
            tempPhiPhi = 0;

            for l = 1:N
                if l ~= j
                    tempAphi = tempAphi + A(l)*sin(phi(j)-phi(l)) * B(k,l);
                    tempPhiPhi = tempPhiPhi + A(l)*cos(phi(j)-phi(l))*B(j,l);
                end
            end

            C_Aphi(j,j) = -2 * A(j) * tempAphi;
            C_phiphi(j,j) = -A(j) * tempPhiPhi;
        end
    end
end

C_AA = C_AA * diameter;
C_Aphi = C_Aphi * diameter;
C_phiphi = C_phiphi * diameter;

end