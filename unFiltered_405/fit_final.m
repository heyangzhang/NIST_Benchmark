%% ============================================================
%  Global Fit with BOOTSTRAP-by-SETS (block bootstrap over trials)
%  - Solver unchanged: fmincon + MultiStart
%  - Each bootstrap replicate samples blocks (sets) of trials, not rows
%  - BOOT provides per-parameter Std & 95% CI and prediction bands
% ============================================================

clear; clc;

%% ------------ USER SETTINGS ------------
DOWNSAMPLE   = 1;        % keep every Nth row per trial (set to 1 to disable)
nStartsMS    = 15;       % MultiStart starts for the main fit (per resin)
nBoot        = 150;      % number of bootstrap resamples (per resin)
SET_SIZE     = 1;        % <-- size of a "set" (block) of trials to resample
                         % e.g., 2 or 3; must be >=1 and <= #trials

MaxFE_main   = 4000;     % fmincon MaxFunctionEvaluations for main fit
MaxFE_boot   = 1500;     % fmincon MaxFunctionEvaluations for each bootstrap refit
OptTol       = 1e-4;     % fmincon OptimalityTolerance
StepTol      = 1e-6;     % fmincon StepTolerance
rng(42);                 % reproducibility

% Parameter bounds and initial guess: params = [k1, m, n]
lb = [1,   0.3, 0.3];
ub = [500, 3.0, 5.0];
x0 = [10,  1.2, 2];

% Residual weights
w1 = 1;  % rate residual
w2 = 0;  % alpha residual (leave 0 for speed, set 1 to include)

% Use script's directory
folderPath = pwd;
fileList = dir(fullfile(folderPath, '*_combined_trials.xlsx'));
if isempty(fileList)
    error('No "*_combined_trials.xlsx" files found in %s', folderPath);
end

% Optimizer options (solver unchanged)
opts_main = optimoptions('fmincon', ...
    'Display','off', ...
    'MaxFunctionEvaluations', MaxFE_main, ...
    'OptimalityTolerance',    OptTol, ...
    'StepTolerance',          StepTol);

opts_boot = optimoptions('fmincon', ...
    'Display','off', ...
    'MaxFunctionEvaluations', MaxFE_boot, ...
    'OptimalityTolerance',    OptTol, ...
    'StepTolerance',          StepTol);

ms = MultiStart('UseParallel', true, 'Display', 'off', 'StartPointsToRun', 'all');

% Try to start a pool (optional)
if isempty(gcp('nocreate'))
    try, parpool; catch, warning('Could not start parallel pool.'); end
end

% Collect per-resin summary rows
summaryRows = {};

for f = 1:length(fileList)
    filename = fileList(f).name;
    filepath = fullfile(folderPath, filename);
    resinID = regexp(filename, '(R\d+)', 'tokens', 'once');
    if isempty(resinID)
        warning('Could not parse resin ID from %s. Skipping.', filename);
        continue;
    end
    rid = resinID{1};
    fprintf('\n=======================\nProcessing %s\n=======================\n', rid);

    % ---------- Load and split ----------
    T = readtable(filepath);
    trials = split_trials_by_column(T);  % builds Trial/Time/Alpha/Rate/a_max per trial

    nTrials = numel(trials);
    if SET_SIZE < 1 || SET_SIZE > nTrials
        warning('SET_SIZE=%d invalid for %s (nTrials=%d). Using SET_SIZE=1.', SET_SIZE, rid, nTrials);
        SET_SIZE = 1;
    end

    % ---------- Downsample (row-wise) ----------
    if DOWNSAMPLE > 1
        for i = 1:nTrials
            nrows = height(trials{i});
            if nrows == 0, continue; end
            keep = 1:DOWNSAMPLE:nrows;
            trials{i} = trials{i}(keep, :);  % subset whole table
        end
    end

    % --- Per-trial a_max and resin-level Amax stats ---
    amax_list = nan(nTrials,1);
    for i = 1:nTrials
        ai = trials{i}.a_max;
        amax_list(i) = ai(find(isfinite(ai),1,'first'));
        if isempty(amax_list(i)) || ~isfinite(amax_list(i))
            amax_list(i) = max(trials{i}.Alpha, [], 'omitnan');
        end
    end
    amax_mean_resin = mean(amax_list, 'omitnan');
    amax_std_resin  = std(amax_list,  'omitnan');

    % --- Optional onset file (mean minus outliers) ---
    onsetFile = fullfile(folderPath, sprintf('%s_onset_times.xlsx', rid));
    meanOnset = NaN;
    if exist(onsetFile, 'file')
        try
            Ton = readtable(onsetFile);
            vn = lower(string(Ton.Properties.VariableNames));
            cIdx = find(vn == "onsettime_abs_s" | vn == "onset_time_abs_s" | vn == "onsettime_s", 1, 'first');
            if isempty(cIdx) && width(Ton) >= 2, cIdx = 2; end
            if ~isempty(cIdx)
                onsetVals = Ton{:, cIdx};
                onsetVals = onsetVals(isfinite(onsetVals));
                if ~isempty(onsetVals)
                    cleanVals = rmoutliers(onsetVals);
                    meanOnset = mean(cleanVals, 'omitnan');
                end
            end
        catch ME
            warning('Could not read onset file for %s: %s', rid, ME.message);
        end
    end

    % ---------- Objective: sum of per-trial residuals ----------
    cost_fun = @(params) sum(arrayfun(@(i) ...
        trial_residual_combined(trials{i}, params, w1, w2), 1:nTrials));

    problem = createOptimProblem('fmincon', ...
        'x0', x0, 'objective', cost_fun, 'lb', lb, 'ub', ub, 'options', opts_main);

    % ---------- Main fit (fmincon + MultiStart, unchanged solver) ----------
    [p_fit, err] = run(ms, problem, nStartsMS);
    fprintf('Best fit [k1 m n] = [%.6g %.6g %.6g], error = %.6g\n', p_fit, err);

    % ---------- BOOTSTRAP-by-SETS (block bootstrap over trials) ----------
    fprintf('Bootstrapping by sets: %d reps, set size = %d (cluster blocks of trials)...\n', nBoot, SET_SIZE);
    p_boot = nan(nBoot, 3);

    % Precompute trial blocks (partition indices into consecutive sets)
    trialIdxAll = 1:nTrials;

    parfor b = 1:nBoot
        try
            % Random permutation of trials, then cut into blocks of size SET_SIZE
            perm = trialIdxAll(randperm(nTrials));
            blocks = mat2cell(perm, 1, [repmat(SET_SIZE,1,floor(nTrials/SET_SIZE)), mod(nTrials,SET_SIZE)]);
            if isempty(blocks{end}), blocks(end) = []; end

            % Number of blocks needed to (roughly) reconstitute original size
            nBlocks = numel(blocks);
            nDraw   = max(1, ceil(nTrials / SET_SIZE));

            % Sample blocks WITH replacement
            drawIdx = randi(nBlocks, [1, nDraw]);
            selTrials = [blocks{drawIdx}];  % concatenate selected blocks
            % If we selected fewer trials than original, top up by random picks
            if numel(selTrials) < nTrials
                selTrials = [selTrials, randsample(trialIdxAll, nTrials - numel(selTrials), true)];
            end
            % If we selected more, trim down
            selTrials = selTrials(1:nTrials);

            trials_boot = trials(selTrials);

            % Refit (warm start at p_fit for speed, single fmincon call)
            cost_fun_boot = @(params) sum(arrayfun(@(i) ...
                trial_residual_combined(trials_boot{i}, params, w1, w2), 1:numel(trials_boot)));
            [p_b, ~] = fmincon(cost_fun_boot, p_fit, [],[],[],[], lb, ub, [], opts_boot);
            p_boot(b,:) = min(max(p_b, lb), ub);
        catch
            p_boot(b,:) = [NaN NaN NaN];
        end
    end

    good = all(isfinite(p_boot),2);
    p_boot = p_boot(good,:);
    if isempty(p_boot)
        warning('No valid bootstrap samples for %s. Writing NaNs.', rid);
    end

    % ---- BOOT statistics (used everywhere) ----
    Std_BOOT    = std(p_boot, 0, 1, 'omitnan');
    CI_BOOT     = prctile(p_boot, [2.5 97.5], 1);  % 2.5/97.5 percentiles
    CI95_BOOT_L = CI_BOOT(1,:).';
    CI95_BOOT_U = CI_BOOT(2,:).';

    % ---------- Prediction band from bootstrap ensemble ----------
    % Common time base across trials
    t_all = unique(sort(cell2mat(cellfun(@(Ti) Ti.Time(:), trials, 'UniformOutput', false))));
    t_all = t_all(:);

    % Use median initial alpha across trials as representative initial condition
    alpha0_vec = cellfun(@(Ti) Ti.Alpha(1), trials);
    alpha0_bar = median(alpha0_vec, 'omitnan');

    % Simulate alpha(t) for a subset of bootstrap parameter vectors (cap to save time)
    nEnsemble = min(size(p_boot,1), 300);  % cap for speed
    sel = randsample(size(p_boot,1), nEnsemble, false);
    alpha_ENS = nan(numel(t_all), nEnsemble);
    for k = 1:nEnsemble
        [~, alpha_model] = simulate_cure(p_boot(sel(k),:), t_all, alpha0_bar, amax_mean_resin);
        alpha_ENS(:,k) = alpha_model;
    end
    alpha_lo = prctile(alpha_ENS, 2.5, 2);
    alpha_hi = prctile(alpha_ENS,97.5, 2);

    % ---------- Plots ----------
    % (1) Alpha(t) band (bootstrap-by-sets) + trial data + best-fit curve
    fig2 = figure('Visible','off');
    hold on; grid on;
    fill([t_all; flipud(t_all)], [alpha_lo; flipud(alpha_hi)], [0.8 0.9 1.0], ...
        'EdgeColor','none', 'FaceAlpha',0.5);
    for i = 1:nTrials
        plot(trials{i}.Time, trials{i}.Alpha, '.', 'MarkerSize',6);
    end
    [t_best, alpha_best] = simulate_cure(p_fit, t_all, alpha0_bar, amax_mean_resin);
    plot(t_best, alpha_best, 'k-', 'LineWidth',1.5);
    title(sprintf('%s: \\alpha(t) 95%% band (BOOT by sets, size=%d)', rid, SET_SIZE));
    xlabel('Time (s)'); ylabel('\alpha');
    legend({'95% band (BOOT-sets)','Exp (trials)','Best-fit'}, 'Location','best');
    saveas(fig2, fullfile(folderPath, [rid, sprintf('_alpha_bands_BOOT_sets%d.png', SET_SIZE)]));
    close(fig2);

    % (2) Per-trial fit panels (best fit)
    fig = figure('Visible', 'off');
    tiledlayout(ceil(nTrials/2), 2);
    for i = 1:nTrials
        nexttile
        t_fit     = trials{i}.Time(:);
        alpha_exp = trials{i}.Alpha(:);
        alpha0    = alpha_exp(1);

        amax_vec = trials{i}.a_max(:);
        idx = find(isfinite(amax_vec),1,'first');
        if isempty(idx), a_max = max(alpha_exp, [], 'omitnan');
        else, a_max = amax_vec(idx); end

        [~, alpha_model] = simulate_cure(p_fit, t_fit, alpha0, a_max);
        plot(t_fit, alpha_exp, 'o', t_fit, alpha_model, '-');
        title(sprintf('Trial %d', i)); xlabel('Time (s)'); ylabel('\alpha');
        legend('Exp', 'Fit', 'Location', 'best'); grid on;
    end
    saveas(fig, fullfile(folderPath, [rid, '_fit_plot.png']));
    close(fig);

    % ---------- Per-resin BOOT stats CSV ----------
    Param = {'k1';'m';'n'};
    Mean  = p_fit(:);
    Std   = Std_BOOT(:);
    CI95_L = CI95_BOOT_L(:);
    CI95_U = CI95_BOOT_U(:);
    stats_tbl = table(Param, Mean, Std, CI95_L, CI95_U, ...
        'VariableNames', {'Param','Mean','Std_BOOT','CI95_BOOT_L','CI95_BOOT_U'});
    writetable(stats_tbl, fullfile(folderPath, [rid, sprintf('_param_stats_BOOT_sets%d.csv', SET_SIZE)]));

    % Save bootstrap raw samples too
    if ~isempty(p_boot)
        Tboot = array2table(p_boot, 'VariableNames', {'k1','m','n'});
        writetable(Tboot, fullfile(folderPath, [rid, sprintf('_param_bootstrap_samples_sets%d.csv', SET_SIZE)]));
    end

    % ---------- Append global summary row (BOOT-only) ----------
    summaryRows(end+1, :) = { ...
        rid, ...
        p_fit, ...                   % FittedParams [k1 m n]
        err, ...
        meanOnset, ...
        amax_mean_resin, ...
        Std_BOOT(:).', ...           % [Std_k1 Std_m Std_n]
        CI95_BOOT_L(:).', ...        % lower CI per param
        CI95_BOOT_U(:).', ...        % upper CI per param
        amax_std_resin, ...          % Std of a_max
        SET_SIZE ...                 % record the block size used
    }; %#ok<SAGROW>
end

% ---------- Global summary table (BOOT-by-sets only) ----------
res = cell2table(summaryRows, ...
    'VariableNames', {'ResinID','FittedParams','TotalError', ...
                      'MeanOnset_s','Amax_mean', ...
                      'ParamStd_BOOT','CI95_BOOT_L','CI95_BOOT_U','Amax_std','BootSetSize'});
writetable(res, fullfile(folderPath, sprintf('global_fit_results_BOOT_sets%d.csv', SET_SIZE)));
disp('Done: BOOTSTRAP-by-SETS std/CI and prediction bands.');

%% ================= Helper Functions =================

function trials = split_trials_by_column(T)
    vn = lower(string(T.Properties.VariableNames));
    trialIdx = find(vn == "trial", 1);
    if isempty(trialIdx), error('Expected "Trial" column.'); end
    timeIdx  = find(vn == "time" | vn == "time_s_since_onset", 1);
    if isempty(timeIdx), error('Expected "Time" or "Time_s_since_onset".'); end
    alphaIdx = find(vn == "alpha", 1);
    if isempty(alphaIdx), error('Expected "Alpha" column.'); end
    rateIdx  = find(vn == "rate" | vn == "dalpha_dt", 1);
    if isempty(rateIdx), rateIdx = NaN; end
    amaxIdx  = find(vn == "a_max" | vn == "amax" | vn == "a__max", 1);
    if isempty(amaxIdx), amaxIdx = NaN; end

    trialNums = unique(T{:,trialIdx});
    trials = cell(numel(trialNums),1);
    for i = 1:numel(trialNums)
        rows = T{:,trialIdx} == trialNums(i);
        Ti = T(rows, :);
        Ti_out = table;
        Ti_out.Trial = Ti{:, trialIdx};
        Ti_out.Time  = Ti{:, timeIdx};
        Ti_out.Alpha = Ti{:, alphaIdx};
        if ~isnan(rateIdx), Ti_out.Rate = Ti{:, rateIdx};
        else, Ti_out.Rate = nan(height(Ti),1); end
        if ~isnan(amaxIdx), Ti_out.a_max = Ti{:, amaxIdx};
        else, Ti_out.a_max = nan(height(Ti),1); end
        trials{i} = Ti_out;
    end
end

function [t_out, alpha_out] = simulate_cure(params, tspan, alpha0, a_max)
    % params = [k1, m, n]
    k1 = params(1); m = params(2); n = params(3);
    ode_fun = @(t, a) reaction_rate(a, k1, m, n, a_max);
    try
        tspan = tspan(:);
        [t_out, alpha_out] = ode45(ode_fun, tspan, alpha0);
        if any(isnan(alpha_out)) || any(~isreal(alpha_out)), error('Invalid ODE'); end
    catch
        t_out = tspan; alpha_out = nan(size(tspan));
    end
end

function r = reaction_rate(a, k1, m, n, a_max)
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

    % per-trial a_max (first finite), fallback to max(alpha)
    amax_vec = trial.a_max(:);
    idx = find(isfinite(amax_vec),1,'first');
    if isempty(idx), a_max = max(alpha_exp, [], 'omitnan');
    else, a_max = amax_vec(idx); end

    alpha0 = alpha_exp(1);

    try
        % rate residual (model rate at measured alpha)
        rate_model = reaction_rate(alpha_exp, params(1), params(2), params(3), a_max);
        r1 = sum((rate_exp - rate_model).^2,'omitnan');

        % alpha residual via ODE (optional)
        if w2 ~= 0
            [t_model, alpha_model] = simulate_cure(params, t, alpha0, a_max);
            if numel(alpha_model) ~= numel(alpha_exp)
                alpha_model = interp1(t_model, alpha_model, t, 'linear', 'extrap');
            end
            r2 = sum((alpha_exp - alpha_model).^2,'omitnan');
        else
            r2 = 0;
        end

        if any(isnan([r1 r2])) || any(~isreal([r1 r2]))
            err = 1e6;
        else
            err = w1 * r1 + w2 * r2;
        end
    catch
        err = 1e6;
    end
end
