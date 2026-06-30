close all
load s170c_senstivity.mat
plot(s170c_sensitivty(:,1), s170c_sensitivty(:,2))
hold on;
xlabel("Wavelength (nm)")
ylabel("Responsitivity (A/W)")
% Example x-values where you want the vertical lines
x1 = 375;
x2 = 460;

% Get y-limits of current axes for shading range
ylims = ylim;

% Shade the area between x1 and x2
X_shade = [x1 x2 x2 x1];
Y_shade = [ylims(1) ylims(1) ylims(2) ylims(2)];
fill(X_shade, Y_shade, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); % gray shaded area

% Add vertical lines
xline(x1, 'r--', 'LineWidth', 1.5);
xline(x2, 'r--', 'LineWidth', 1.5);

% Optional labeling
text(x2+20, ylims(2), "Relevant Wavelengths", 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
hold off;

%%
clear; clc;

% Load the Excel file
filename = 'uv_vis.xlsx';

% Read the full spreadsheet into a table
T = readtable(filename, 'ReadVariableNames', false);

[data, txt, ~] = xlsread('uv_vis.xlsx');
wavelengths = data(1, 2:end);         % wavelengths [nm]
absorbance = data(2:end, 2:end);      % [wavelengths x resins]
resin_names = txt(2:end, 1);           % resin names

% Plot each spectrum
figure; hold on; grid on;
n_spectra = size(absorbance, 1);

for i = 1:n_spectra
    plot(wavelengths, absorbance(i, :), 'DisplayName', [resin_names{i}]);
end

xlabel('Wavelength (nm)', 'FontWeight', 'bold');
ylabel('Absorbance (a.u.)', 'FontWeight', 'bold');
title('UV-Vis Absorbance Spectra');
legend('Location', 'best');
xlim([min(wavelengths), max(wavelengths)]);
