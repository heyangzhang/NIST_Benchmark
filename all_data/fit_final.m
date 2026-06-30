clear; clc;

% === Use script's directory ===
folderPath = pwd;   % current working directory

fileList = dir(fullfile(folderPath, '*_combined_trials.xlsx'));
resinResults = {};

% params = [k1, m, n]
lb = [0,   0.3, 0.3];
ub = [50,  3.0, 5.0];
x0 = [0.01, 1.2, 1.5];

% Weights for residual components
w1 = 1; w2 = 0;

for f = 1:length(fileList)
    filename = fileList(f).name;
    filepath = fullfile(folderPath, filename);
    resinID = regexp(filename, '(R\d+)', 'tokens', 'once');
    if isempty(resinID)
        warning('Could not parse resin ID from %s. Skipping.', filename);
        continue;
    end
    fprintf('Processing %s...\n', resinID{1});
    
    % Load data and split by trial column
    T = readtable(filepath);
    trials = split_trials_by_column(T);

    % --- Compute per-trial a_max and resin-level average a_max ---
    amax_list = nan(length(trials),1);
    for i = 1:length(trials)
        amax_list(i) = max(trials{i}.Alpha, [], 'omitnan');
    end
    amax_mean_resin = mean(amax_list, 'omitnan');   % <-- average a_max per resin

    % Define cost function (sum over trials)
    cost_fun = @(params) sum(arrayfun(@(i) ...
        trial_residual_combined(trials{i}, params, w1, w2), 1:length(trials)));

    % Setup fmincon problem
    opts = optimoptions('fmincon', 'Display','off', 'MaxFunctionEvaluations', 1e4);
    problem = createOptimProblem('fmincon', ...
        'x0', x0, 'objective', cost_fun, ...
        'lb', lb, 'ub', ub, 'options', opts);
    ms = MultiStart('UseParallel', true, 'Display', 'iter', ...
        'StartPointsToRun', 'all');

    % Run MultiStart
    if isempty(gcp('nocreate')), parpool; end
    [p_fit, err] = run(ms, problem, 50);

    % --- Read onset file and compute max onset (s) ---
    onsetFile = fullfile(folderPath, sprintf('%s_onset_times.xlsx', resinID{1}));
    maxOnset = NaN;
    if exist(onsetFile, 'file')
        try
            Ton = readtable(onsetFile);
            if any(strcmpi(Ton.Properties.VariableNames, 'OnsetTime_abs_s'))
                onsetVals = Ton{:, 'OnsetTime_abs_s'};
            elseif width(Ton) >= 2
                onsetVals = Ton{:, 2};
            else
                onsetVals = [];
            end
            if ~isempty(onsetVals)
                maxOnset = max(onsetVals, [], 'omitnan');
            end
        catch ME
            warning('Could not read onset file for %s: %s', resinID{1}, ME.message);
        end
    else
        warning('Onset file not found for %s: %s', resinID{1}, onsetFile);
    end

    % Store result  (now also saving Amax_mean)
    resinResults{end+1,1} = resinID{1};
    resinResults{end,2}   = p_fit;
    resinResults{end,3}   = err;
    resinResults{end,4}   = maxOnset;        % Max onset (s)
    resinResults{end,5}   = amax_mean_resin; % <-- Average a_max for this resin

    % --- Plot fits (per trial) ---
    fig = figure('Visible', 'off');  % Don't show plot
    tiledlayout(length(trials), 1);

    for i = 1:length(trials)
        nexttile
        t_fit    = trials{i}.Time(:);
        alpha_exp= trials{i}.Alpha(:);
        alpha0   = alpha_exp(1);
        a_max    = amax_list(i);  % per-trial a_max for modeling

        [t_model, alpha_model] = simulate_cure(p_fit, t_fit, alpha0, a_max);
        plot(t_fit, alpha_exp, 'o', t_model, alpha_model, '-');
        title(sprintf('Trial %d', i)); xlabel('Time (s)'); ylabel('\alpha');
        legend('Exp', 'Fit', 'Location', 'best');
        grid on;
    end

    % Save plot
    saveas(fig, fullfile(folderPath, [resinID{1}, '_fit_plot.png']));
    close(fig);
end

% Save results (now includes Amax_mean)
res = cell2table(resinResults, ...
    'VariableNames', {'ResinID', 'FittedParams', 'TotalError', 'MaxOnset_s', 'Amax_mean'});
writetable(res, fullfile(folderPath, 'global_fit_results.csv'));
disp('Done fitting all resins.');

%% === Helper Functions ===

function trials = split_trials_by_column(T)
    trialNums = unique(T{:,1});
    trials = cell(length(trialNums),1);
    for i = 1:length(trialNums)
        trials{i} = T(T{:,1} == trialNums(i), :);
        trials{i}.Properties.VariableNames = {'Trial', 'Time', 'Alpha', 'Rate'};
    end
end

function [t_out, alpha_out] = simulate_cure(params, tspan, alpha0, a_max)
    % params = [k1, m, n]
    k1 = params(1); m = params(2); n = params(3);
    ode_fun = @(t, a) reaction_rate(a, k1, m, n, a_max);
    try
        tspan = tspan(:);
        [t_out, alpha_out] = ode45(ode_fun, tspan, alpha0);
        if any(isnan(alpha_out)) || any(~isreal(alpha_out))
            error('Invalid ODE result');
        end
    catch
        t_out = tspan;
        alpha_out = nan(size(tspan));
    end
end

function r = reaction_rate(a, k1, m, n, a_max)
    % Clamp to physical ranges
    a     = max(a, 0);
    term1 = max(a, 0).^m;
    term2 = max(a_max - a, 0).^n;
    term3 = max(1 - a, 0);
    r = k1 .* term1 .* term2 .* term3;
end

function err = trial_residual_combined(trial, params, w1, w2)
    t = trial.Time(:);
    alpha_exp = trial.Alpha(:);
    rate_exp  = trial.Rate(:);
    alpha0    = alpha_exp(1);
    a_max     = max(alpha_exp, [], 'omitnan');   % per-trial a_max from data

    try
        % Rate residual
        rate_model = reaction_rate(alpha_exp, params(1), params(2), params(3), a_max);
        r1 = sum((rate_exp - rate_model).^2);

        % Alpha residual via ODE
        [t_model, alpha_model] = simulate_cure(params, t, alpha0, a_max);
        if length(alpha_model) ~= length(alpha_exp)
            alpha_model = interp1(t_model, alpha_model, t, 'linear', 'extrap');
        end
        r2 = sum((alpha_exp - alpha_model).^2);

        if any(isnan([r1 r2])) || any(~isreal([r1 r2]))
            err = 1e6;
        else
            err = w1 * r1 + w2 * r2;
        end
    catch
        err = 1e6;
    end
end
