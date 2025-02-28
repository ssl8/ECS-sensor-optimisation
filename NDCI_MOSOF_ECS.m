%% Multi-Objective Sensor Optimization using MOSOF Principles (NDCI-Centric)
% This script implements a MOSOF-inspired sensor selection framework where
% NDCI is central. All cases use the OEM sensor names and OEM NDCI values.
% For Airlines and MRO the additional attributes (cost, coverage, reliability/efficiency)
% are assigned randomly (within specified ranges) so that the objectives are not fully correlated.
%
% Four objectives are defined in each case.
% A multiobjective genetic algorithm (MOGA) is used to select exactly 2 sensors.
% Then, six pairwise 2D scatter plots are produced showing every solution (each solution represents a sensor pair)
% with sensor pair annotations. Feasible solutions are marked with blue circles,
% infeasible ones with red crosses, and the selected optimal pair is highlighted with a green star.
%
% Change 'scenario' to 'OEM', 'Airlines', or 'MRO' as needed.
%
% Author: Burak SUSLU (modified)
% Date:   2025-02-22

clc; clear; close all;
rng(42);  % set random seed for reproducibility

%% Select Scenario: 'OEM', 'Airlines', or 'MRO'
scenario = 'MRO';  % <-- Change this value as needed

switch scenario
    case 'OEM'
        % --- OEM Sensor Pool (9 sensors) ---
        sensors = [ ...
            struct('name','ThiSHX','cost',640, 'accuracy',0.95, 'ndci',34.66, 'mtbf',7100, 'compatibility',8.2);
            struct('name','ThoPHX','cost',570, 'accuracy',0.91, 'ndci',28.02, 'mtbf',6800, 'compatibility',7.8);
            struct('name','ThiRHX','cost',520, 'accuracy',0.96, 'ndci',22.26, 'mtbf',6400, 'compatibility',6.1);
            struct('name','TiT',   'cost',490, 'accuracy',0.93, 'ndci',15.77, 'mtbf',6700, 'compatibility',5.5);
            struct('name','ThiCHX','cost',460, 'accuracy',0.90, 'ndci',12.95, 'mtbf',6300, 'compatibility',4.8);
            struct('name','ToT',   'cost',430, 'accuracy',0.89, 'ndci',10.39, 'mtbf',6000, 'compatibility',3.9);
            struct('name','P_O',   'cost',410, 'accuracy',0.94, 'ndci',7.16,  'mtbf',5800, 'compatibility',2.7);
            struct('name','TciRHX','cost',390, 'accuracy',0.87, 'ndci',6.48,  'mtbf',5600, 'compatibility',1.3);
            struct('name','TciCHX','cost',370, 'accuracy',0.92, 'ndci',0.71,  'mtbf',5400, 'compatibility',0.9)  ...
        ];
        names    = {sensors.name};
        cost     = [sensors.cost];
        accuracy = [sensors.accuracy];
        ndci     = [sensors.ndci];
        mtbf     = [sensors.mtbf];
        compatibility = [sensors.compatibility];
        
        % OEM constraints (relaxed so several pairs are feasible)
        budget = 1150;          % Total cost <= 1150
        minNDCI = 40;           % Sum of ndci >= 40
        minAvgAccuracy = 0.92;  % Average accuracy >= 0.92
        minMTBF = 12000;        % Sum of mtbf >= 12000
        
        % Normalize attributes for objectives
        norm_ndci = ndci / max(ndci);
        norm_mtbf = mtbf / max(mtbf);
        norm_compat = compatibility / max(compatibility);
        
        % Define OEM objectives (4 objectives)
        % f1: Performance = - (0.4*norm_ndci + 0.6*accuracy)
        % f2: Cost (sum of cost)
        % f3: Reliability = - normalized mtbf
        % f4: Compatibility = - normalized compatibility
        objFcn = @(x) objectiveOEM(x, norm_ndci, accuracy, cost, norm_mtbf, norm_compat);
        consFcn = @(x) constraintOEM(x, cost, budget, ndci, minNDCI, accuracy, minAvgAccuracy, mtbf, minMTBF);
        
        objNames = {'Performance','Cost','Reliability','Compatibility'};
        
    case 'Airlines'
        % --- Airlines Sensor Pool using OEM names & ndci values ---
        % Define OEM names and OEM ndci values
        names = {'ThiSHX','ThoPHX','ThiRHX','TiT','ThiCHX','ToT','P_O','TciRHX','TciCHX'};
        ndciOEM = [34.66, 28.02, 22.26, 15.77, 12.95, 10.39, 7.16, 6.48, 0.71];
        % For each sensor, assign random values (within specified ranges) for cost, coverage, reliability
        nS = numel(names);
        cost = zeros(1,nS); 
        coverage = zeros(1,nS);
        reliability = zeros(1,nS);
        for i = 1:nS
            cost(i) = 260 + randi([-20,20]);  % cost in [240,280]
            coverage(i) = 0.90 + rand()*0.08;   % coverage in [0.90,0.98]
            reliability(i) = 5400 + randi([-300,300]); % reliability in [5100,5700]
        end
        % Build structure array
        sensors = struct('name',{},'cost',{},'coverage',{},'reliability',{},'ndci',{});
        for i = 1:nS
            sensors(i) = struct('name',names{i}, 'cost', cost(i), 'coverage', coverage(i), ...
                'reliability', reliability(i), 'ndci', ndciOEM(i));
        end
        
        % Airlines constraints (relaxed)
        budget = 700;           % increased budget
        minCoverage = 1.5;      % combined coverage >= 1.5
        minReliability = 10000; % total reliability >= 10000
        minNDCI = 30;           % sum of ndci >= 30
        
        % Normalize measures
        norm_ndci = ndciOEM / max(ndciOEM);
        norm_reliability = reliability / max(reliability);
        
        % Define Airlines objectives (4 objectives)
        % f1: Performance = - (0.5*norm_ndci + 0.5*coverage)
        % f2: Cost (sum of cost)
        % f3: Reliability = - normalized reliability
        % f4: Benefit-to-Cost = -((norm_ndci+coverage) ./ cost)
        objFcn = @(x) objectiveAirlines(x, cost, coverage, norm_ndci, norm_reliability);
        consFcn = @(x) constraintAirlines(x, cost, budget, coverage, minCoverage, reliability, minReliability, ndciOEM, minNDCI);
        
        objNames = {'Performance','Cost','Reliability','Benefit-to-Cost'};
        
    case 'MRO'
        % --- MRO Sensor Pool using OEM names & ndci values ---
        names = {'ThiSHX','ThoPHX','ThiRHX','TiT','ThiCHX','ToT','P_O','TciRHX','TciCHX'};
        ndciOEM = [34.66, 28.02, 22.26, 15.77, 12.95, 10.39, 7.16, 6.48, 0.71];
        nS = numel(names);
        cost = zeros(1,nS);
        coverage = zeros(1,nS);
        efficiency = zeros(1,nS);
        for i = 1:nS
            cost(i) = 260 + randi([-20,20]);  % cost in [240,280]
            coverage(i) = 0.90 + rand()*0.08;   % coverage in [0.90,0.98]
            efficiency(i) = 0.80 + rand()*0.15; % efficiency in [0.80,0.95]
        end
        sensors = struct('name',{},'cost',{},'coverage',{},'efficiency',{},'ndci',{});
        for i = 1:nS
            sensors(i) = struct('name',names{i}, 'cost', cost(i), 'coverage', coverage(i), ...
                'efficiency', efficiency(i), 'ndci', ndciOEM(i));
        end
        
        % MRO constraints (relaxed)
        budget = 700;           % increased budget
        minCoverage = 1.5;      % combined coverage >= 1.5
        minEfficiency = 1.5;    % total efficiency >= 1.5
        minNDCI = 30;           % sum of ndci >= 30
        
        norm_ndci = ndciOEM / max(ndciOEM);
        
        % Define MRO objectives (4 objectives)
        % f1: Performance = - (0.5*norm_ndci + 0.5*coverage)
        % f2: Cost (sum of cost)
        % f3: Efficiency = - efficiency
        % f4: Benefit-to-Cost = -((norm_ndci+coverage+efficiency) ./ cost)
        objFcn = @(x) objectiveMRO(x, cost, coverage, norm_ndci, efficiency);
        consFcn = @(x) constraintMRO(x, cost, budget, coverage, minCoverage, efficiency, minEfficiency, ndciOEM, minNDCI);
        
        objNames = {'Performance','Cost','Efficiency','Benefit-to-Cost'};
        
    otherwise
        error('Unknown scenario. Please choose OEM, Airlines, or MRO.');
end

%% GA Setup: The decision variable is a binary vector (exactly 2 sensors must be selected)
nvars = length(sensors);       
lb = zeros(1, nvars);
ub = ones(1, nvars);

options = optimoptions('gamultiobj', ...
    'PopulationSize', 200, ...
    'MaxGenerations', 500, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotpareto, ... % optional GA plot (can be commented out)
    'FunctionTolerance', 1e-6, ...
    'MutationFcn', @mutationadaptfeasible, ...
    'CreationFcn', @customCreationFcn, ...
    'CrossoverFcn', @crossoverintermediate, ...
    'SelectionFcn', @selectiontournament);

[x, fval, exitflag, output, ~, ~] = gamultiobj(objFcn, nvars, [], [], [], [], lb, ub, consFcn, options);

%% Validate and classify Pareto solutions
selectedPairIndex = [];
validSolutionFound = false;
feasibleIdx = []; % indices of solutions that satisfy constraints
infeasibleIdx = [];
for i = 1:size(x,1)
    selSensors = find(round(x(i,:)) > 0.5);
    [c, ceq] = consFcn(x(i,:));
    if all(c <= 0) && all(abs(ceq) <= 1e-6)
        feasibleIdx(end+1) = i;
        if length(selSensors) == 2 && ~validSolutionFound
            validSolutionFound = true;
            selectedPairIndex = i;
            fprintf('Optimal Sensor Configuration (%s):\n', scenario);
            disp(x(i,:));
            fprintf('Objective Values:\n');
            disp(fval(i,:));
            fprintf('Selected Sensors: %s and %s\n', names{selSensors(1)}, names{selSensors(2)});
        end
    else
        infeasibleIdx(end+1) = i;
    end
end

if ~validSolutionFound
    disp('No valid sensor pair found that meets all constraints.');
end

%% Figure: Pairwise Scatter Plots for All 6 Objective Combinations
% Create subplots for each 2-objective combination from the 4 objectives.
objPair = nchoosek(1:4,2);  % will create 6 subplots
figure('Name',sprintf('Pairwise Objective Scatter Plots (%s Scenario)',scenario),...
    'Color','w','Position',[50 50 1400 800]);

for k = 1:size(objPair,1)
    subplot(2,3,k);
    i1 = objPair(k,1); i2 = objPair(k,2);
    % Plot all Pareto solutions (each solution represents a sensor pair)
    scatter(fval(:,i1), fval(:,i2), 50, 'k', 'filled', 'MarkerFaceAlpha',0.4);
    hold on;
    % Overlay feasible solutions with blue circles
    if ~isempty(feasibleIdx)
        scatter(fval(feasibleIdx,i1), fval(feasibleIdx,i2), 70, 'b', 'o', 'LineWidth',1.5);
    end
    % Overlay infeasible solutions with red crosses
    if ~isempty(infeasibleIdx)
        scatter(fval(infeasibleIdx,i1), fval(infeasibleIdx,i2), 70, 'r', 'x', 'LineWidth',1.5);
    end
    % Highlight the selected optimal pair with a green star
    if validSolutionFound
        scatter(fval(selectedPairIndex,i1), fval(selectedPairIndex,i2), 120, 'g', 'p', 'filled');
    end
    xlabel(objNames{i1}, 'FontSize',12);
    ylabel(objNames{i2}, 'FontSize',12);
    title(sprintf('%s vs %s', objNames{i1}, objNames{i2}), 'FontSize',12);
    grid on;
    % Annotate each point with its sensor pair names
    for j = 1:size(x,1)
        sel = find(round(x(j,:)) > 0.5);
        if length(sel)==2
            pairStr = sprintf('%s & %s', names{sel(1)}, names{sel(2)});
            text(fval(j,i1), fval(j,i2), pairStr, 'FontSize',8, 'Color','m',...
                'VerticalAlignment','bottom','HorizontalAlignment','right');
        end
    end
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local Function Definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Custom Creation Function: Generate individuals with exactly 2 sensors selected.
function Population = customCreationFcn(GenomeLength, ~, options)
    Population = zeros(options.PopulationSize, GenomeLength);
    for i = 1:options.PopulationSize
        sel = randperm(GenomeLength, 2);
        Population(i, sel) = 1;
    end
end

%% OEM Objective Function (4 objectives: Performance, Cost, Reliability, Compatibility)
function f = objectiveOEM(x, norm_ndci, accuracy, cost, norm_mtbf, norm_compat)
    x = round(x);
    f1 = - (0.4 * sum(x .* norm_ndci) + 0.6 * sum(x .* accuracy));
    f2 = sum(x .* cost);
    f3 = - sum(x .* norm_mtbf);
    f4 = - sum(x .* norm_compat);
    f = [f1, f2, f3, f4];
end

%% OEM Constraint Function
function [c, ceq] = constraintOEM(x, cost, budget, ndci, minNDCI, accuracy, minAvgAccuracy, mtbf, minMTBF)
    x = round(x);
    c = [ sum(x .* cost) - budget, ...                
          minNDCI - sum(x .* ndci), ...                 
          minAvgAccuracy - (sum(x .* accuracy)/sum(x)), ... 
          minMTBF - sum(x .* mtbf) ];
    ceq = sum(x) - 2;  % Exactly 2 sensors selected
end

%% Airlines Objective Function (4 objectives: Performance, Cost, Reliability, Benefit-to-Cost)
function f = objectiveAirlines(x, cost, coverage, norm_ndci, norm_reliability)
    x = round(x);
    f1 = - (0.5 * sum(x .* norm_ndci) + 0.5 * sum(x .* coverage));
    f2 = sum(x .* cost);
    f3 = - sum(x .* norm_reliability);
    benefit = (norm_ndci + coverage);
    f4 = - sum(x .* (benefit ./ cost));
    f = [f1, f2, f3, f4];
end

%% Airlines Constraint Function
function [c, ceq] = constraintAirlines(x, cost, budget, coverage, minCoverage, reliability, minReliability, ndci, minNDCI)
    x = round(x);
    c = [ sum(x .* cost) - budget, ...
          minCoverage - sum(x .* coverage), ...
          minReliability - sum(x .* reliability), ...
          minNDCI - sum(x .* ndci) ];
    ceq = sum(x) - 2;
end

%% MRO Objective Function (4 objectives: Performance, Cost, Efficiency, Benefit-to-Cost)
function f = objectiveMRO(x, cost, coverage, norm_ndci, efficiency)
    x = round(x);
    f1 = - (0.5 * sum(x .* norm_ndci) + 0.5 * sum(x .* coverage));
    f2 = sum(x .* cost);
    f3 = - sum(x .* efficiency);
    benefit = (norm_ndci + coverage + efficiency);
    f4 = - sum(x .* (benefit ./ cost));
    f = [f1, f2, f3, f4];
end

%% MRO Constraint Function
function [c, ceq] = constraintMRO(x, cost, budget, coverage, minCoverage, efficiency, minEfficiency, ndci, minNDCI)
    x = round(x);
    c = [ sum(x .* cost) - budget, ...
          minCoverage - sum(x .* coverage), ...
          minEfficiency - sum(x .* efficiency), ...
          minNDCI - sum(x .* ndci) ];
    ceq = sum(x) - 2;
end
