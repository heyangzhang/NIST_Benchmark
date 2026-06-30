clear; clc;

% === Use script's directory ===
folderPath = pwd;   % current working directory
fileList   = dir(fullfile(folderPath, '*.xlsx'));

% === Define Resin IDs ===
resinIDs = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6'};
groupedFiles = containers.Map;

% === Onset detection params ===
alphaThreshold = 0.01;    % primary threshold on alpha
minRise       = 1e-5;     % derivative fallback threshold
dt            = 0.1;      % interpolation step [s]
discardT      = 15;       % discard first 15 s BEFORE onset detection
postWindow    = 60;       % keep first 30 s AFTER onset

% === Group files by resin ID ===
names = {fileList.name};
for r = 1:numel(resinIDs)
    rid = resinIDs{r};
    groupedFiles(rid) = fileList(contains(names, rid) & ...
                                 contains(names, 'Kinetics', 'IgnoreCase', true));
end

% === Loop over each resin group ===
for g = 1:numel(resinIDs)
    resinID = resinIDs{g};
    groupFiles = groupedFiles(resinID);

    if isempty(groupFiles)
        warning('No files found for resin %s.', resinID);
        continue;
    end

    allData    = [];   % [Trial, Time_s_since_onset, Alpha, dAlpha_dt]
    onsetTimes = [];   % [Trial, OnsetTime_abs_s]

    figure('Name', ['Onset Detection - ' resinID]); hold on;

    for i = 1:numel(groupFiles)
        fileName = groupFiles(i).name;
        filePath = fullfile(folderPath, fileName);

        try
            % === Read data (time in col1, alpha in col2) ===
            try
                data = readmatrix(filePath);
            catch
                T = readtable(filePath);
                data = table2array(T(:,1:2));
            end

            if size(data,2) < 2
                warning('Expected at least 2 columns (time, alpha) in %s. Skipping.', fileName);
                continue;
            end

            time  = data(:,1);
            alpha = data(:,2);

            % Drop NaNs and sort by time
            valid = isfinite(time) & isfinite(alpha);
            time  = time(valid);   alpha = alpha(valid);
            [time, sortIdx] = sort(time);
            alpha = alpha(sortIdx);

            if numel(time) < 3
                warning('Too few points in %s (Trial %d). Skipping.', fileName, i);
                continue;
            end

            % === Interpolate to uniform grid ===
            tUniform    = (min(time):dt:max(time))';
            alphaInterp = interp1(time, alpha, tUniform, 'pchip');

            % === Discard the first 15 seconds BEFORE onset detection ===
            chopMask = tUniform >= discardT;
            if ~any(chopMask)
                warning('No data remains after discarding first %g s in %s (Trial %d).', discardT, fileName, i);
                continue;
            end
            tChop  = tUniform(chopMask);
            aChop  = alphaInterp(chopMask);

            % Derivative on chopped segment
            da_dt  = gradient(aChop, dt);

            % === Detect onset on chopped segment ===
            idxThr  = find(aChop > alphaThreshold, 1, 'first');
            onsetIdxLocal = idxThr;      % conservative earliest consistent point

            % Absolute onset time (relative to original start)
            onsetTimeAbs = tChop(onsetIdxLocal);
            onsetTimes   = [onsetTimes; i, onsetTimeAbs];

            % === Trim data from onset and keep first postWindow seconds ===
            tFromOnset    = tChop - onsetTimeAbs;
            keepFromOnset = tFromOnset >= 0 & tFromOnset <= postWindow;

            tTrim     = tFromOnset(keepFromOnset);
            aTrim     = aChop(keepFromOnset);
            da_dtTrim = da_dt(keepFromOnset);

            if isempty(tTrim)
                warning('No samples within %g s post-onset for %s (Trial %d).', postWindow, resinID, i);
                continue;
            end

            trialCol  = repmat(i, numel(tTrim), 1);
            allData   = [allData; [trialCol, tTrim, aTrim, da_dtTrim]]; %#ok<AGROW>

            % === Plot for verification ===
            plot(tUniform, alphaInterp, '-', 'DisplayName', sprintf('Trial %d', i));
            plot(onsetTimeAbs, aChop(onsetIdxLocal), 'ro', 'MarkerFaceColor', 'r');

        catch ME
            warning('Error reading file: %s — %s', fileName, ME.message);
        end
    end

    if isempty(allData)
        warning('No valid data found for resin %s.', resinID);
        continue;
    end

    % === Save combined trial data ===
    outputFile = fullfile(folderPath, sprintf('%s_combined_trials.xlsx', resinID));
    header     = {'Trial', 'Time_s_since_onset', 'Alpha', 'dAlpha_dt'};
    outData    = [header; num2cell(allData)];
    writecell(outData, outputFile);
    fprintf('✅ Combined data written for %s: %s\n', resinID, outputFile);

    % === Save onset times (absolute, from original start) ===
    onsetFile   = fullfile(folderPath, sprintf('%s_onset_times.xlsx', resinID));
    onsetHeader = {'Trial', 'OnsetTime_abs_s'};
    onsetData   = [onsetHeader; num2cell(onsetTimes-15)];
    writecell(onsetData, onsetFile);
    fprintf('📄 Onset times saved for %s: %s\n', resinID, onsetFile);

    % === Finalize plot ===
    title(['Onset Detection for ' resinID]);
    xlabel('Time (s)'); ylabel('Alpha');
    legend('show'); grid on;
end
