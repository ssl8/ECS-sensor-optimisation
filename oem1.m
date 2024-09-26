clc; clear; close all;

% Define the sensors with explicit qualifications
sensors = [
    struct('name', 'Sensor 1', 'cost', 250, 'accuracy', 0.95, 'infoGain', 80, 'mtbf', 5500);
    struct('name', 'Sensor 2', 'cost', 300, 'accuracy', 0.92, 'infoGain', 75, 'mtbf', 5400);
    struct('name', 'Sensor 3', 'cost', 220, 'accuracy', 0.93, 'infoGain', 85, 'mtbf', 5600);
    struct('name', 'Sensor 4', 'cost', 280, 'accuracy', 0.96, 'infoGain', 90, 'mtbf', 5700);
    struct('name', 'Sensor 5', 'cost', 260, 'accuracy', 0.91, 'infoGain', 70, 'mtbf', 5300);
    struct('name', 'Sensor 6', 'cost', 240, 'accuracy', 0.94, 'infoGain', 88, 'mtbf', 5800);
    struct('name', 'Sensor 7', 'cost', 270, 'accuracy', 0.97, 'infoGain', 95, 'mtbf', 5900);
    struct('name', 'Sensor 8', 'cost', 230, 'accuracy', 0.90, 'infoGain', 65, 'mtbf', 5200);
    struct('name', 'Sensor 9', 'cost', 290, 'accuracy', 0.89, 'infoGain', 60, 'mtbf', 5100);
    struct('name', 'Sensor 10', 'cost', 250, 'accuracy', 0.98, 'infoGain', 94, 'mtbf', 6000);
    struct('name', 'Sensor 11', 'cost', 276, 'accuracy', 0.94, 'infoGain', 85, 'mtbf', 5550);
    struct('name', 'Sensor 12', 'cost', 236, 'accuracy', 0.85, 'infoGain', 91, 'mtbf', 5450);
    struct('name', 'Sensor 13', 'cost', 226, 'accuracy', 0.92, 'infoGain', 74, 'mtbf', 5650);
    struct('name', 'Sensor 14', 'cost', 246, 'accuracy', 0.93, 'infoGain', 84, 'mtbf', 5750);
    struct('name', 'Sensor 15', 'cost', 266, 'accuracy', 0.86, 'infoGain', 93, 'mtbf', 5350);
    struct('name', 'Sensor 16', 'cost', 242, 'accuracy', 0.87, 'infoGain', 87, 'mtbf', 5850);
    struct('name', 'Sensor 17', 'cost', 252, 'accuracy', 0.80, 'infoGain', 77, 'mtbf', 5950);
    struct('name', 'Sensor 18', 'cost', 235, 'accuracy', 0.91, 'infoGain', 87, 'mtbf', 5250);
    struct('name', 'Sensor 19', 'cost', 225, 'accuracy', 0.88, 'infoGain', 90, 'mtbf', 5150);
    struct('name', 'Sensor 20', 'cost', 255, 'accuracy', 0.89, 'infoGain', 88, 'mtbf', 5800)
];

numSensors = length(sensors);

% Extract sensor properties
cost = [sensors.cost];
accuracy = [sensors.accuracy];
infoGain = [sensors.infoGain];
mtbf = [sensors.mtbf];

% Define the constraints
budget =480; % Total budget in dollars
minInfoGain = 135; % Minimum information gain
minAvgAccuracy = 0.9; % Minimum average accuracy
minTotalMTBF = 11000; % Minimum total MTBF

% Define the options for the genetic algorithm
options = optimoptions('gamultiobj', ...
    'PopulationSize', 200, ...
    'MaxGenerations', 500, ...
    'ParetoFraction', 0.35, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotpareto, ...
    'FunctionTolerance', 1e-6, ...
    'MutationFcn', @mutationadaptfeasible, ... % Use mutationadaptfeasible
    'CreationFcn', @customCreationFcn, ... % Use custom creation function
    'CrossoverFcn', @crossoverintermediate, ... % Use intermediate crossover
    'SelectionFcn', @selectiontournament); % Use tournament selection

% Number of variables is the number of sensors
nvars = numSensors;
lb = zeros(1, nvars);  % Lower bounds of variables
ub = ones(1, nvars);  % Upper bounds of variables

% Run the multi-objective genetic algorithm
[x, fval, exitflag, output, population, scores] = gamultiobj(@(x)objectiveFcn(x, accuracy, cost, infoGain, mtbf), ...
    nvars, [], [], [], [], lb, ub, @(x)constraintFcn(x, cost, budget, infoGain, minInfoGain, accuracy, minAvgAccuracy, mtbf, minTotalMTBF), options);

% Validate the selection and print results

% Validate the selection and print results
validSolutionFound = false;
selectedPairIndex = [];
for i = 1:size(x, 1)
    selectedPair = find(x(i, :) > 0.5);
    if length(selectedPair) == 2
        validSolutionFound = true;
        selectedPairIndex = i;
        disp('Optimal Sensor Configuration:');
        disp(x(i, :));
        disp('Objective Function Values:');
        disp(fval(i, :));
        disp('Selected Pair of Sensors:');
        disp(selectedPair);
        break;
    end
end

if ~validSolutionFound
    disp('No valid pair of sensors found that meet all constraints.');
end



% Diagnostic code to analyze constraint violations
for i = 1:size(population, 1)
    [c, ceq] = constraintFcn(population(i, :), cost, budget, infoGain, minInfoGain, accuracy, minAvgAccuracy, mtbf, minTotalMTBF);
    disp(['Individual ', num2str(i), ' Constraint Violations:']);
    disp(c');
end

% Plot the Pareto front with labels for Cost vs Accuracy/Precision and Information Gain
figure;
plot3(fval(:, 2), -fval(:, 1), -fval(:, 3), 'bo');
xlabel('Total Cost');
ylabel('Negative Total Accuracy');
zlabel('Negative Information Gain');
title('Pareto Front for Sensor Selection');
grid on;
hold on;

% Annotate each point with its corresponding sensors
for i = 1:size(fval, 1)
    sensorsSelected = find(x(i, :) > 0.5); % Adjust the threshold to ensure only binary selections
    if length(sensorsSelected) == 2 % Ensure only solutions with exactly 2 sensors are labeled
        text(fval(i, 2), -fval(i, 1), -fval(i, 3), ['Sensors: ', num2str(sensorsSelected)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
end

% Highlight the selected pair
if validSolutionFound
    plot3(fval(selectedPairIndex, 2), -fval(selectedPairIndex, 1), -fval(selectedPairIndex, 3), 'r*', 'MarkerSize', 10);
    legend('Pareto Front', 'Selected Pair');
end

% Custom creation function to generate an initial population
function Population = customCreationFcn(GenomeLength, ~, options)
    Population = randi([0, 1], options.PopulationSize, GenomeLength);
    % Ensure initial population contains feasible individuals
    for i = 1:options.PopulationSize
        while sum(Population(i, :)) ~= 2
            Population(i, :) = randi([0, 1], 1, GenomeLength);
        end
    end
end

% Define the multi-objective function
function f = objectiveFcn(x, accuracy, cost, infoGain, ~)
    % Ensure x is binary
    x = round(x);

    % Initialize the output vector f
    f = zeros(1, 3); % Each individual returns a 3-element row vector

    % Objective 1: Maximize sensor performance (Accuracy)
    f(1) = -sum(x .* accuracy);

    % Objective 2: Minimize total cost of ownership
    f(2) = sum(x .* cost);

    % Objective 3: Maximize information gain
    f(3) = -sum(x .* infoGain);
end

% Define the constraints function for the genetic algorithm
function [c, ceq] = constraintFcn(x, cost, budget, infoGain, minInfoGain, accuracy, minAvgAccuracy, mtbf, minTotalMTBF)
    % Ensure x is binary
    x = round(x);

    % Inequality constraints
    c = [
        sum(x .* cost) - budget; % Total cost should not exceed budget
        -sum(x .* infoGain) + minInfoGain; % Total information gain should be at least minInfoGain
        -sum(x .* accuracy) / sum(x) + minAvgAccuracy; % Average accuracy should be at least minAvgAccuracy
        -sum(x .* mtbf) + minTotalMTBF % Total MTBF should be at least minTotalMTBF
    ];
    
    % Equality constraint: Exactly 2 sensors should be selected
    ceq = sum(x) - 2;
end
