clear; clc;

%% User input
P_meas = 0.000165353; % Measured power in W from PM100D (adjust this)

%% Load files
specData = readmatrix('ps_405f.txt');         % [nm, arbitrary units]
respData = readmatrix('s170c_responstivity.txt');     % [nm, A/W]

lambda_spec = specData(3:end,1);
S_rel = specData(3:end,2);

lambda_resp = respData(3:end,1);
R_lambda = respData(3:end,2);

%% Interpolate responsivity to match spectrometer wavelengths
R_interp = interp1(lambda_resp, R_lambda, lambda_spec, 'linear', 'extrap');

%% Normalize the spectrum (area under curve = 1)
S_rel_norm = S_rel / trapz(lambda_spec, S_rel);

%% Compute absolute spectral power distribution (W/nm)
S_abs = P_meas * S_rel_norm;

%% Plot result
figure;
plot(lambda_spec, S_abs, 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Spectral Power (W/nm)');
title('Calibrated Spectral Power Distribution');
grid on;

%% Save to file
output = [lambda_spec, S_abs];
writematrix(output, 'calibrated_spectrum.txt', 'Delimiter', 'tab');

disp('✅ Calibrated spectrum saved to "calibrated_spectrum.txt"');
