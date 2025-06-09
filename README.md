# ECS Sensor Optimisation

This repository contains a small collection of MATLAB scripts for multi-objective
sensor selection using genetic algorithms. Each script configures a set of
sensors with different performance metrics and runs `gamultiobj` from MATLAB's
Global Optimization Toolbox to find an optimal subset.

## Prerequisites

You can clone the repository from GitHub using:

```bash
git clone https://github.com/username/ecs-sensor-optimisation.git
```

- MATLAB with the Global Optimization Toolbox installed.

## Usage

Run any of the MATLAB scripts from the MATLAB command window or using the
`run` command. For example:

```matlab
run('air1.m');
```

The scripts will display optimization progress in the command window and plot
the Pareto front of solutions. After completion, the selected sensors and
their metrics are printed.

## Scripts

- **`air1.m`** – Optimizes for coverage, cost, and reliability.
- **`mro1.m`** – Optimizes for coverage, cost, and data processing efficiency.
- **`oem1.m`** – Optimizes for accuracy, cost, information gain, and MTBF.

Each script defines its own list of sensors and constraints in the header
section. Adjust these parameters to match your specific scenario.

## License

This project is released under [CC0 1.0](LICENSE).
