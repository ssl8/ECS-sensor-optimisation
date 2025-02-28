function multiFaultComparison_Refined()
    % MULTIFAULTCOMPARISON_REFINED
    %
    % Key Features (updated):
    % 1) Per-scenario NDCI figures (SP, S, U, and NDCI), each plotted once.
    % 2) Global (aggregated) mRMR ranking => single figure of overall results.
    % 3) Final figure comparing top-K sensors by NDCI vs. by mRMR.
    % 4) Coverage-based subset selection "per scenario" (union taken across scenarios).
    % 5) Additional figures displaying the selected sensor sets for each scenario 
    %    (by both NDCI and mRMR) and an overall summary figure.
    %
    % Assumptions:
    %   - Each sheet has row #1 => healthy (Deg=1), subsequent rows => fault states.
    %   - SP is defined by normalizing differences with the "healthy range",
    %     and S by normalizing with (1 - Deg) per fault row.
    %   - Coverage of a scenario is defined as the average NDCI of the chosen sensor
    %     subset divided by the maximum average NDCI (if using all sensors).
    %
    % Author: Burak SUSLU (modified)
    % Date:   2025-02-20

    clc;  clear all;

    %% ================ FILE & BASIC PARAMETERS ============================
    filename = 'SESAC1.xlsx';
    if ~isfile(filename)
        error('File "%s" not found.', filename);
    end

    % Tolerance to consider a sensor "constant" across all faults
    constThreshold    = 1e-3;
    coverageThreshold = 0.95;   % e.g., 95% coverage threshold

    % Attempt to read list of sheet names:
    try
        sheetList = sheetnames(filename);
    catch
        [~, sheetList] = xlsfinfo(filename);  % fallback for older MATLAB versions
    end
    numSheets = numel(sheetList);

    %% =========== ACCUMULATORS FOR GLOBAL (ALL SCENARIOS) ===============
    allScenarioNDCI    = [];     % rows = scenarios, columns = sensors
    allSensorNames     = {};     % will fix sensor order from first valid scenario
    scenarioLabels     = {};     % store sheet names or parsed "scenario" IDs
    aggregatedFaultData = [];    % for global mRMR
    aggregatedSheetIDs  = [];    % numeric scenario ID per row of aggregatedFaultData

    validCount = 0;  % how many sheets are processed successfully

    % For final coverage approach:
    scenarioNDICell  = {};     % store each scenario's [faultRows x sensors] NDCI
    % (Also, store per-scenario selected sensor sets for additional figures)
    scenarioSelectedSensors_NDCI = {};
    scenarioSelectedSensors_mRMR = {};

    %% ================= PROCESS EACH SCENARIO ============================
    for s = 1:numSheets
        sheetName = sheetList{s};
        [ndciVals, sensorNames, spVals, sVals, uVals, ...
         faultData, scenarioNDCI, scenarioValid] = processScenario(filename, sheetName);

        if ~scenarioValid
            warning('Skipping sheet "%s" due to invalid or insufficient data.', sheetName);
            continue;
        end
        validCount = validCount + 1;

        % If the sensor name list is empty so far, store it. Otherwise, check consistency:
        if isempty(allSensorNames)
            allSensorNames = sensorNames;
        elseif ~isequal(allSensorNames, sensorNames)
            warning('Sensor mismatch in sheet "%s" => using initial scenario sensor ordering.', sheetName);
        end

        scenarioLabels{validCount} = sheetName; %#ok<AGROW>
        allScenarioNDCI = [allScenarioNDCI; ndciVals]; %#ok<AGROW>

        % Accumulate fault data for global mRMR
        aggregatedFaultData = [aggregatedFaultData; faultData]; %#ok<AGROW>
        aggregatedSheetIDs  = [aggregatedSheetIDs; repmat(validCount, size(faultData,1),1)]; %#ok<AGROW>

        % Store scenario-level NDCI array for coverage approach:
        scenarioNDICell{validCount} = scenarioNDCI; %#ok<AGROW>

        % === Plot each scenario's NDCI figure exactly once here ===
        plotScenarioResults(sheetName, sensorNames, spVals, sVals, uVals, ndciVals);
    end

    if validCount == 0
        warning('No valid sheets processed. Exiting.');
        return;
    end

    %% ================ FILTER CONSTANT SENSORS ===========================
    sensorRange = max(aggregatedFaultData, [], 1) - min(aggregatedFaultData, [], 1);
    validSensors = (sensorRange > constThreshold);

    removedSensors = allSensorNames(~validSensors);
    if ~isempty(removedSensors)
        fprintf('\nExcluding constant sensors (range <= %.3g): %s\n', ...
            constThreshold, strjoin(removedSensors, ', '));
    end

    finalSensorNames = allSensorNames(validSensors);
    % Filter out invalid sensors from scenario NDCI as well
    allScenarioNDCI_valid = allScenarioNDCI(:, validSensors);

    %% ================ OVERALL NDCI RANKING =============================
    % We take the average across scenarios for each sensor
    avgNDCI = mean(allScenarioNDCI_valid, 1);
    [sortedNDCIvals, sortIdx] = sort(avgNDCI, 'descend');
    sortedSensorsNDCI = finalSensorNames(sortIdx);

    fprintf('\n=== Overall Sensor Ranking by NDCI ===\n');
    for i = 1:numel(sortedSensorsNDCI)
        fprintf('%2d) %s => Avg NDCI = %.4f\n', i, sortedSensorsNDCI{i}, sortedNDCIvals(i));
    end

    % Plot overall NDCI ranking
    figure('Name','Overall Sensor Ranking by NDCI','Color','white','Position',[150 150 800 400]);
    bar(sortedNDCIvals, 'FaceColor', [0.6 0.2 0.7]);
    title('Overall Sensor Ranking (NDCI)','FontSize',14,'FontWeight','bold');
    xlabel('Sensors (Ranked)'); ylabel('Average NDCI');
    set(gca, 'XTick', 1:numel(sortedSensorsNDCI), 'XTickLabel', sortedSensorsNDCI, 'XTickLabelRotation', 45);

    %% ================= GLOBAL mRMR RANKING =============================
    dataMatrix = aggregatedFaultData(:, validSensors);
    [sortedSensorsMRMR, mrmrScores] = doGlobalMRMR(dataMatrix, finalSensorNames, aggregatedSheetIDs);

    %% ============== COMPARE TOP-K NDCI VS. TOP-K MRMR ==================
    topK = min(5, numel(sortedSensorsNDCI));  % pick 5 or fewer if fewer sensors
    topNDCI = sortedSensorsNDCI(1:topK);
    topMRMR = sortedSensorsMRMR(1:topK);

    figure('Name','Comparison: Top Sensors by NDCI vs. mRMR',...
           'Color','white','Position',[200 200 1000 400]);
    subplot(1,2,1);
    bar(categorical(topNDCI), 1:topK, 'FaceColor', [0.3 0.3 0.9]);
    title('Top Sensors by NDCI','FontWeight','bold');
    ylabel('Rank'); set(gca, 'YTick', []);
    subplot(1,2,2);
    bar(categorical(topMRMR), 1:topK, 'FaceColor', [0.9 0.4 0.2]);
    title('Top Sensors by mRMR','FontWeight','bold');
    ylabel('Rank'); set(gca, 'YTick', []);

    %% =========== AUTOMATED COVERAGE (PER SCENARIO) =====================
    % For each scenario, we define coverage = (average NDCI using subset) /
    % (max average NDCI if using ALL sensors). Then we take the union
    % of the minimal subsets across all scenarios.
    minimalSetNDCI = {};
    minimalSetMRMR = {};

    for sc = 1:validCount
        scenarioName = scenarioLabels{sc};
        scenarioNDCI = scenarioNDICell{sc}(:, validSensors);  
        
        % Pick a minimal subset from the sortedSensorsNDCI
        subsetNDCI = pickMinimalScenarioSubset(scenarioNDCI, finalSensorNames, ...
            sortedSensorsNDCI, coverageThreshold);
        % Pick a minimal subset from the sortedSensorsMRMR
        subsetMRMR = pickMinimalScenarioSubset(scenarioNDCI, finalSensorNames, ...
            sortedSensorsMRMR, coverageThreshold);
        
        % Store per-scenario sensor selections (for additional figures)
        scenarioSelectedSensors_NDCI{sc} = subsetNDCI;
        scenarioSelectedSensors_mRMR{sc} = subsetMRMR;
        
        % Update overall union sets
        minimalSetNDCI = union(minimalSetNDCI, subsetNDCI);
        minimalSetMRMR = union(minimalSetMRMR, subsetMRMR);

        fprintf('\nScenario #%d (%s) minimal coverage by NDCI => ', sc, scenarioName);
        disp(subsetNDCI');
        fprintf('Scenario #%d (%s) minimal coverage by mRMR => ', sc, scenarioName);
        disp(subsetMRMR');
    end

    fprintf('\n=== FINAL Minimal Sensor Subset (NDCI-based) across all scenarios ===\n');
    disp(minimalSetNDCI');
    fprintf('\n=== FINAL Minimal Sensor Subset (mRMR-based) across all scenarios ===\n');
    disp(minimalSetMRMR');

    %% ===== NEW FIGURE: OVERALL SELECTED SENSOR SETS =====
    figure('Name','Overall Minimal Sensor Subset','Color','white','Position',[300 300 800 300]);
    subplot(1,2,1);
    title('Overall Selected Sensors by NDCI','FontSize',12,'FontWeight','bold');
    text(0.05, 0.5, strjoin(minimalSetNDCI, ', '), 'FontSize', 12);
    axis off;
    subplot(1,2,2);
    title('Overall Selected Sensors by mRMR','FontSize',12,'FontWeight','bold');
    text(0.05, 0.5, strjoin(minimalSetMRMR, ', '), 'FontSize', 12);
    axis off;

    %% ===== NEW FIGURES: PER-SCENARIO SELECTED SENSOR SETS =====
    for sc = 1:validCount
        figure('Name', sprintf('Scenario %d Selected Sensor Sets', sc),...
               'Color','white','Position',[350 350 800 300]);
        subplot(1,2,1);
        title(sprintf('Scenario %d: Sensors by NDCI', sc), 'FontSize',12,'FontWeight','bold');
        text(0.05, 0.5, strjoin(scenarioSelectedSensors_NDCI{sc}, ', '), 'FontSize', 12);
        axis off;
        subplot(1,2,2);
        title(sprintf('Scenario %d: Sensors by mRMR', sc), 'FontSize',12,'FontWeight','bold');
        text(0.05, 0.5, strjoin(scenarioSelectedSensors_mRMR{sc}, ', '), 'FontSize', 12);
        axis off;
    end
end

%% =========================== PROCESS SCENARIO ============================
function [ndciAvg, sensorNames, spVals, sVals, uVals, faultData, scenarioNDCI, valid] = processScenario(filename, sheetName)
    % PROCESSSCENARIO
    %
    % 1) Reads the sheet => row #1 is healthy, row #2..end => faults.
    % 2) Defines:
    %    SP => (absDiff) / (range(healthyVal) + eps)
    %    S  => (absDiff) / ((1 - deg(i)) + eps)
    %    U  => uniqueness among fault rows
    %    NDCI => (SP + S + U)/3
    % 3) Returns:
    %    ndciAvg   => [1 x numSensors], averaged across fault samples
    %    scenarioNDCI => [numFaultSamples x numSensors], row-wise NDCI 
    %    valid     => indicates whether scenario is valid

    valid = false;
    ndciAvg = [];
    spVals  = [];
    sVals   = [];
    uVals   = [];
    faultData  = [];
    sensorNames= {};
    scenarioNDCI = [];

    try
        T = readtable(filename, 'Sheet', sheetName);
    catch ME
        warning('Cannot read sheet "%s": %s', sheetName, ME.message);
        return;
    end

    if height(T) < 2 || width(T) < 2
        warning('Sheet "%s": insufficient data.', sheetName);
        return;
    end

    deg = T{:,1};           % first column => deg
    rawData = T{:,2:end};   % sensors
    sensorNames = T.Properties.VariableNames(2:end);
    nSensors    = numel(sensorNames);

    % row #1 => healthy
    healthyVals = rawData(1,:);
    % the rest => faults
    faultIdx = 2:height(T);
    if isempty(faultIdx)
        warning('No fault rows in "%s".', sheetName);
        return;
    end

    faultData = rawData(faultIdx,:);
    faultDeg  = deg(faultIdx);
    numFaults = numel(faultIdx);

    % Define sigma for SP as the range of the healthy row:
    sigmaSP = range(healthyVals) + 1e-12;

    spM = zeros(numFaults, nSensors);
    sM  = zeros(numFaults, nSensors);

    for i = 1:numFaults
        diffVal = abs(faultData(i,:) - healthyVals);
        spM(i,:) = diffVal / sigmaSP;   % SP: scaled by sigmaSP
        denomS = abs(1 - faultDeg(i)) + 1e-12;
        sM(i,:) = diffVal ./ denomS;      % S: scaled by (1 - deg)
    end

    spVals = mean(spM,1);
    sVals  = mean(sM,1);

    % Uniqueness among faultData
    uVals = computeUniqueness(faultData);

    % NDCI => average of (SP + S + U)/3 across fault samples
    ndciAvg = (spVals + sVals + uVals) / 3;

    % Compute row-wise NDCI (for coverage check)
    rowWiseNDCI = zeros(numFaults, nSensors);
    for i = 1:numFaults
        rowSP = spM(i,:);
        rowS  = sM(i,:);
        rowU  = uVals; % same for each row in this approach
        rowWiseNDCI(i,:) = (rowSP + rowS + rowU) / 3;
    end
    scenarioNDCI = rowWiseNDCI;
    valid = true;
end

%% ====================== UNIQUENESS CALCULATION ==========================
function U = computeUniqueness(faultData)
    [numFaults, nSensors] = size(faultData);
    if numFaults < 2
        U = ones(1, nSensors);
        return;
    end
    D = pdist(faultData', 'euclidean');
    Dmat = squareform(D);
    avgDist = zeros(1, nSensors);
    for j = 1:nSensors
        avgDist(j) = mean(Dmat(j, [1:j-1, j+1:end]));
    end
    maxDist = max(avgDist);
    if maxDist < 1e-12
        U = zeros(1, nSensors);
    else
        U = avgDist / maxDist;
    end
end

%% =================== SCENARIO FIGURE (SP, S, U, NDCI) ==================
function plotScenarioResults(sheetName, sensorNames, sp, sMetric, u, ndci)
    figure('Name',['Scenario: ', sheetName, ' (NDCI)'], 'Color','white', 'Position',[200 200 1100 600]);

    subplot(2,2,1);
    bar(sp, 'FaceColor', [0.2 0.6 0.8]);
    title('Separation Power (SP)','FontSize',12,'FontWeight','bold');
    set(gca, 'XTick', 1:numel(sensorNames), 'XTickLabel', sensorNames, 'XTickLabelRotation',45);
    ylabel('SP');

    subplot(2,2,2);
    bar(sMetric, 'FaceColor', [0.9 0.4 0.2]);
    title('Sensitivity (S)','FontSize',12,'FontWeight','bold');
    set(gca, 'XTick', 1:numel(sensorNames), 'XTickLabel', sensorNames, 'XTickLabelRotation',45);
    ylabel('S');

    subplot(2,2,3);
    bar(u, 'FaceColor', [0.4 0.8 0.3]);
    title('Uniqueness (U)','FontSize',12,'FontWeight','bold');
    set(gca, 'XTick', 1:numel(sensorNames), 'XTickLabel', sensorNames, 'XTickLabelRotation',45);
    ylabel('U');

    subplot(2,2,4);
    bar(ndci, 'FaceColor', [0.6 0.2 0.7]);
    title('Overall NDCI','FontSize',12,'FontWeight','bold');
    set(gca, 'XTick', 1:numel(sensorNames), 'XTickLabel', sensorNames, 'XTickLabelRotation',45);
    ylabel('NDCI');

    sgtitle(['Scenario: ', sheetName, '  (SP, S, U, NDCI)']);
    
    % Display a table of metrics in the command window
    T = table(sensorNames', sp', sMetric', u', ndci', 'VariableNames',{'Sensor','SP','S','U','NDCI'});
    disp(['--- Metrics for Scenario: ', sheetName, ' ---']);
    disp(T);
end

%% ========================= GLOBAL MRMR =================================
function [sortedSensors, mrmrScores] = doGlobalMRMR(dataMatrix, sensorNames, sheetIDs)
    if exist('fscmrmr','file') ~= 2
        warning('fscmrmr not found. Skipping global mRMR.');
        sortedSensors = sensorNames;
        mrmrScores = zeros(size(sensorNames));
        return;
    end

    % Build table for classification using scenario IDs as labels
    scenarioCats = categorical(sheetIDs);
    T = array2table(dataMatrix, 'VariableNames', sensorNames);
    [featIdx, featScores] = fscmrmr(T, scenarioCats);
    sortedSensors = sensorNames(featIdx);
    mrmrScores = featScores;

    fprintf('\n=== Sensor Ranking by mRMR ===\n');
    for i = 1:numel(sortedSensors)
        fprintf('%2d) %s => Score = %.4f\n', i, sortedSensors{i}, mrmrScores(i));
    end

    topK = min(5, numel(sortedSensors));
    topMRMR = sortedSensors(1:topK);
    fprintf('\n--- Top %d Sensors by mRMR ---\n', topK);
    disp(topMRMR');

    % Plot the global mRMR ranking
    figure('Name','Global mRMR Ranking','Color','white','Position',[250 250 900 400]);
    bar(mrmrScores, 'FaceColor', [0.2 0.6 0.8]);
    title('Global Sensor Ranking by mRMR','FontSize',14,'FontWeight','bold');
    set(gca, 'XTick', 1:numel(sortedSensors), 'XTickLabel', sortedSensors, 'XTickLabelRotation',45);
    xlabel('Sensors'); ylabel('mRMR Score');
end

%% ============ PICK MINIMAL SUBSET FOR A SINGLE SCENARIO ================
function minimalSet = pickMinimalScenarioSubset(scenarioNDCI, sensorNames, sensorRank, coverageThreshold)
    % PICKMINIMALSCENARIOSUBSET
    %   scenarioNDCI => [numFaultSamples x numSensors]
    %   sensorNames  => {1 x numSensors}
    %   sensorRank   => the sensor names in sorted order (by NDCI or mRMR)
    %   coverageThreshold => fraction (e.g., 0.95)
    %
    % We pick sensors (in sensorRank order) until the average NDCI of the selected
    % sensors reaches coverageThreshold of the maximum possible (i.e., using all sensors).

    if isempty(scenarioNDCI)
        minimalSet = {};
        return;
    end

    nSensorsAll = size(scenarioNDCI,2);
    allSensorsVal = mean(mean(scenarioNDCI,2));
    if allSensorsVal < 1e-12
        minimalSet = {};
        return;
    end

    usedIdx = [];
    coverageFrac = 0.0;
    for r = 1:numel(sensorRank)
        sName = sensorRank{r};
        idx = find(strcmp(sensorNames, sName), 1);
        if isempty(idx), continue; end
        usedIdx = [usedIdx, idx]; %#ok<AGROW>
        subVal = mean(mean(scenarioNDCI(:, usedIdx), 2));
        coverageFrac = subVal / allSensorsVal;
        if coverageFrac >= coverageThreshold
            break;
        end
    end

    minimalSet = sensorNames(usedIdx);
end
