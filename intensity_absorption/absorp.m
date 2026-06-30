clear; clc;

%% Parameters
l = 1.0; % path length in cm
ln10 = log(10);

%% Load calibrated spectrum [wavelength (nm), W/nm]
calSpec = readmatrix('calibrated_405f.txt');
lambda_cal = calSpec(:,1);
power_spec = calSpec(:,2);

% Total spectral power (for normalization)
P0 = trapz(lambda_cal, power_spec);

%% Load absorbance data from Excel
[data, txt, ~] = xlsread('uv_vis.xlsx');
wavelengths = data(1, 2:end)';         % wavelengths [nm]
absorbance = data(2:end, 2:end)';      % [wavelengths x resins]
resin_names = txt(2:end, 1);           % resin names

%% Interpolate spectral power to match absorbance grid
power_interp = interp1(lambda_cal, power_spec, wavelengths, 'linear', 0);

%% Compute power-weighted average absorbance, then alpha_eff
num_resins = size(absorbance, 2);
alpha_eff = zeros(num_resins, 1);
avg_abs = zeros(num_resins, 1);

for i = 1:num_resins
    A = absorbance(:,i);
    weighted_avg_A = trapz(wavelengths, A .* power_interp) / trapz(wavelengths, power_interp);
    avg_abs(i) = weighted_avg_A;
    alpha_eff(i) = ln10 * weighted_avg_A / l;  % units: cm^-1
end

%% Display results
fprintf('\n%-10s | %-15s | %-15s\n', 'Resin', 'Avg Absorbance', 'Alpha_eff (cm^-1)');
fprintf('%s\n', repmat('-',1,46));
for i = 1:num_resins
    fprintf('%-10s | %13.4f   | %13.5f\n', resin_names{i}, avg_abs(i), alpha_eff(i));
end

%% Optional: save to CSV
T = table(string(resin_names), avg_abs, alpha_eff, ...
    'VariableNames', {'Resin', 'AvgAbsorbance', 'AlphaEff_cm_inv'});
writetable(T, 'alpha_eff_from_absorbance.csv');
