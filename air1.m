clc; clear; close all;

rng(42);

% Define the sensors with explicit qualifications
sensors = [
    struct('name', 'Sensor 1', 'cost', 250, 'coverage', 0.95, 'reliability', 5500);
    struct('name', 'Sensor 2', 'cost', 300, 'coverage', 0.92, 'reliability', 5400);
    struct('name', 'Sensor 3', 'cost', 220, 'coverage', 0.93, 'reliability', 5600);
    struct('name', 'Sensor 4', 'cost', 280, 'coverage', 0.96, 'reliability', 5700);
    struct('name', 'Sensor 5', 'cost', 260, 'coverage', 0.91, 'reliability', 5300);
    struct('name', 'Sensor 6', 'cost', 240, 'coverage', 0.94, 'reliability', 5800);
    struct('name', 'Sensor 7', 'cost', 270, 'coverage', 0.97, 'reliability', 5900);
    struct('name', 'Sensor 8', 'cost', 230, 'coverage', 0.90, 'reliability', 5200);
    struct('name', 'Sensor 9', 'cost', 290, 'coverage', 0.89, 'reliability', 5100);
    struct('name', 'Sensor 10', 'cost', 250, 'coverage', 0.98, 'reliability', 6000);
    struct('name', 'Sensor 11', 'cost', 276, 'coverage', 0.84, 'reliability', 5550);
    struct('name', 'Sensor 12', 'cost', 306, 'coverage', 0.85, 'reliability', 5450);
    struct('name', 'Sensor 13', 'cost', 226, 'coverage', 0.82, 'reliability', 5650);
    struct('name', 'Sensor 14', 'cost', 286, 'coverage', 0.83, 'reliability', 5750);
    struct('name', 'Sensor 15', 'cost', 266, 'coverage', 0.86, 'reliability', 5350);
    struct('name', 'Sensor 16', 'cost', 245, 'coverage', 0.87, 'reliability', 5850);
    struct('name', 'Sensor 17', 'cost', 275, 'coverage', 0.80, 'reliability', 5950);
    struct('name', 'Sensor 18', 'cost', 235, 'coverage', 0.81, 'reliability', 5250);
    struct('name', 'Sensor 19', 'cost', 295, 'coverage', 0.88, 'reliability', 5150);
    struct('name', 'Sensor 20', 'cost', 255, 'coverage', 0.89, 'reliability', 4800)
];

numSensors = length(sensors);

% Extract sensor properties
cost = [sensors.cost];
coverage = [sensors.coverage];
reliability = [sensors.reliability];

% Define the constraints
budget = 800; % Total budget in dollars
minCoverage = 2.7; % Minimum sensor coverage
minReliability = 15000; % Minimum total MTBF
numSelectedSensors = 3; % Number of sensors to select

% Define the options for the genetic algorithm
options = optimoptions('gamultiobj', ...
    'PopulationSize', 500, ...
    'MaxGenerations', 1500, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotpareto, ...
    'FunctionTolerance', 1e-6, ...
    'MutationFcn', @mutationadaptfeasible, ...
    'CreationFcn', @(GenomeLength, FitnessFcn, options) customCreationFcn(GenomeLength, FitnessFcn, options, numSelectedSensors), ...
    'CrossoverFcn', @crossoverintermediate, ...
    'SelectionFcn', @selectiontournament);

% Number of variables is the number of sensors
nvars = numSensors;
lb = zeros(1, nvars);  % Lower bounds of variables
ub = ones(1, nvars);  % Upper bounds of variables

% Run the multi-objective genetic algorithm
[x, fval, exitflag, output, population, scores] = gamultiobj(@(x)objectiveFcn(x, coverage, cost, reliability), ...
    nvars, [], [], [], [], lb, ub, @(x)constraintFcn(x, cost, budget, coverage, minCoverage, reliability, minReliability, numSelectedSensors), options);

% Validate the selection and print results
validSolutionFound = false;
selectedPairIndex = [];
for i = 1:size(x, 1)
    selectedPair = find(x(i, :) > 0.5);
    if length(selectedPair) == numSelectedSensors
        validSolutionFound = true;
        selectedPairIndex = i;
        disp('Optimal Sensor Configuration:');
        disp(x(i, :));
        disp('Objective Function Values:');
        disp(fval(i, :));
        disp('Selected Sensors:');
        disp(selectedPair);
        break;
    end
end

if ~validSolutionFound
    disp('No valid set of sensors found that meet all constraints.');
end

% Diagnostic code to analyze constraint violations
for i = 1:size(population, 1)
    [c, ceq] = constraintFcn(population(i, :), cost, budget, coverage, minCoverage, reliability, minReliability, numSelectedSensors);
    disp(['Individual ', num2str(i), ' Constraint Violations:']);
    disp(c');
end

% Plot the Pareto front with labels for Cost vs Coverage and Reliability
figure;
plot3(fval(:, 2), -fval(:, 1), -fval(:, 3), 'bo');
xlabel('Total Cost');
ylabel('Negative Total Coverage');
zlabel('Negative Reliability');
title('Pareto Front for Sensor Selection');
grid on;
hold on;

% Annotate each point with its corresponding sensors
for i = 1:size(fval, 1)
    sensorsSelected = find(x(i, :) > 0.5); % Adjust the threshold to ensure only binary selections
    if length(sensorsSelected) == numSelectedSensors % Ensure only solutions with exactly 3 sensors are labeled
        text(fval(i, 2), -fval(i, 1), -fval(i, 3), ['Sensors: ', num2str(sensorsSelected)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end

% Highlight the selected pair
if validSolutionFound
    plot3(fval(selectedPairIndex, 2), -fval(selectedPairIndex, 1), -fval(selectedPairIndex, 3), 'r*', 'MarkerSize', 10);
    legend('Pareto Front', 'Selected Pair');
end

% Nested custom creation function to generate an initial population
function Population = customCreationFcn(GenomeLength, ~, options, numSelectedSensors)
    Population = randi([0, 1], options.PopulationSize, GenomeLength);
    % Ensure initial population contains feasible individuals
    for i = 1:options.PopulationSize
        while sum(Population(i, :)) ~= numSelectedSensors
            Population(i, :) = randi([0, 1], 1, GenomeLength);
        end
    end
end

% Define the multi-objective function
function f = objectiveFcn(x, coverage, cost, reliability)
    % Ensure x is binary
    x = round(x);

    % Initialize the output vector f
    f = zeros(1, 3); % Each individual returns a 3-element row vector

    % Objective 1: Maximize sensor coverage (negative for minimization)
    f(1) = -sum(x .* coverage);

    % Objective 2: Minimize total cost of ownership
    f(2) = sum(x .* cost);

    % Objective 3: Maximize system reliability (negative for minimization)
    f(3) = -sum(x .* reliability);
end

% Define the constraints function for the genetic algorithm
function [c, ceq] = constraintFcn(x, cost, budget, coverage, minCoverage, reliability, minReliability, numSelectedSensors)
    % Ensure x is binary
    x = round(x);

    % Inequality constraints
    c = [
        sum(x .* cost) - budget; % Total cost should not exceed budget
        -sum(x .* coverage) + minCoverage; % Total sensor coverage should be at least minCoverage
        -sum(x .* reliability) + minReliability % Total reliability should be at least minReliability
    ];
    
    % Equality constraint: Exactly numSelectedSensors sensors should be selected
    ceq = sum(x) - numSelectedSensors;
end
