clc; clear; close;

data = readtable('absorption.csv');

% Extract the columns
X = data.X;
Y = data.Y;

set(0, 'DefaultFigureWindowStyle', 'normal') % Change to NORMAL to export

col = styling;

% Create the plot
figure;
semilogx(X, 99*Y, '-', "LineWidth", 1.5, "Color", col(1, :)); 
xlabel('Wavelength [m]');
ylabel('Transmission [%]');
ylim([0, 100])
xlim([X(1), X(end)])
grid on;
colorbar('off')

styling(true, 10, 6, "exports/transmission_windows", false);