clear; clc;

% === Use script's directory ===
folderPath = pwd;   % current working directory
fileList   = dir(fullfile(folderPath, '*.xlsx'));

% === Define Resin IDs ===
resinIDs = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6'};
groupedFiles = containers.Map;

% === Onset detection params ===
alphaThreshold = 0.08;    % primary threshold on alpha
minRise       = 1e-5;    % derivative fallback threshold
dt            = 0.1;     % interpolation step [s]
discardT      = 15;      % discard first 15 s BEFORE onset detection
postWindow    = 15;      % keep first 15 s AFTER onset

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

    allData    = [];   % [Trial, Time_s_since_onset, Alpha, dAlpha_dt, a_max]
    onsetTimes = [];   % [Trial, OnsetTime_abs_s]

    for i = 1:numel(groupFiles)
        fileName = groupFiles(i).name;
        filePath = fullfile(folderPath, fileName);

        try
            for j = 4
                % === Read data ===
                try
                    data = readmatrix(filePath);
                catch
                    T = readtable(filePath);
                    data = table2array(T(:,1:4));
                end

                if size(data,2) < 2
                    warning('Expected at least 2 columns (time, alpha) in %s. Skipping.', fileName);
                    continue;
                end

                time  = data(:,1);
                alpha = data(:,j);

                % Drop NaNs and sort
                valid = isfinite(time) & isfinite(alpha);
                time  = time(valid);   
                alpha = alpha(valid);
                [time, sortIdx] = sort(time);
                alpha = alpha(sortIdx);

                if numel(time) < 3
                    warning('Too few points in %s (Trial %d). Skipping.', fileName, i);
                    continue;
                end

                % === Interpolate to uniform grid ===
                tUniform    = (min(time):dt:max(time))';
                alphaInterp = interp1(time, alpha, tUniform, 'pchip');

                % %%% Compute a_max BEFORE any chopping
                a_max_val = max(alphaInterp);

                % === Discard first 15 s ===
                chopMask = tUniform >= discardT;
                if ~any(chopMask)
                    warning('No data remains after discarding first %g s in %s (Trial %d).', discardT, fileName, i);
                    continue;
                end
                tChop  = tUniform(chopMask);
                aChop  = alphaInterp(chopMask);

                % Derivative
                da_dt  = gradient(aChop, dt);

                % === Detect onset ===
                idxThr  = find(aChop > alphaThreshold, 1, 'first');
                onsetIdxLocal = idxThr;

                onsetTimeAbs = tChop(onsetIdxLocal);
                onsetTimes   = [onsetTimes; i*10+j, onsetTimeAbs];

                % === Trim data from onset and keep postWindow ===
                tFromOnset    = tChop - onsetTimeAbs;
                keepFromOnset = tFromOnset >= 0 & tFromOnset <= postWindow;

                tTrim     = tFromOnset(keepFromOnset);
                aTrim     = aChop(keepFromOnset);
                da_dtTrim = da_dt(keepFromOnset);

                if isempty(tTrim)
                    warning('No samples within %g s post-onset for %s (Trial %d).', postWindow, resinID, i);
                    continue;
                end

                trialCol  = repmat(i*10 + j, numel(tTrim), 1);
                a_maxCol  = repmat(a_max_val, numel(tTrim), 1);   % %%% constant per trial
                allData   = [allData; [trialCol, tTrim, aTrim, da_dtTrim, a_maxCol]]; %#ok<AGROW>
            end
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
    header     = {'Trial', 'Time_s_since_onset', 'Alpha', 'dAlpha_dt', 'a_max'}; %%% added
    outData    = [header; num2cell(allData)];
    writecell(outData, outputFile);
    fprintf('✅ Combined data written for %s: %s\n', resinID, outputFile);

    % === Save onset times ===
    onsetFile   = fullfile(folderPath, sprintf('%s_onset_times.xlsx', resinID));
    onsetHeader = {'Trial', 'OnsetTime_abs_s'};
    onsetData   = [onsetHeader; num2cell(onsetTimes)];
    writecell(onsetData, onsetFile);
    fprintf('📄 Onset times saved for %s: %s\n', resinID, onsetFile);
end
