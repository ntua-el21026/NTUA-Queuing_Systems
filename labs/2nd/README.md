# Lab 2 - M/M/1 and Finite Birth-Death Queues

## Goal

This assignment connects the theoretical M/M/1 formulas with numerical CTMC calculations. It also studies a finite-capacity birth-death queue where the arrival and service rates depend on the state.

## Implementation

The Octave programs cover the computational parts of the assignment:

- `programs/lab2_part2.m`: scans the service rate of an M/M/1 queue with arrival rate `lambda = 10`, computes utilization, response time, mean system population, throughput, and empty-system probability, and plots each metric against `mu`.
- `programs/lab2_part3a.m`: builds the generator matrix for an M/M/1/3 birth-death process with state-dependent arrival rates, computes stationary probabilities, blocking probability, mean number of customers, equilibrium throughput, and transient convergence from the empty state.
- `programs/lab2_part3b.m`: repeats the served-customers and transient analysis for several service rates to show how faster service changes the stationary distribution and convergence behavior.

The scripts use the Octave `queueing` package, mainly `ctmcbd`, `ctmc`, and `qsmm1`.

## Inputs and Outputs

- Input material is in `assignment/`, including the lab statement and starter/help code.
- Generated steady-state and transient plots are committed in `programs/output/`.
- Exported code PDFs are committed in `programs/pdf/`.

## Organization

- `assignment/`: lab statement, introduction material, and starter/help code.
- `programs/`: implemented Octave scripts.
- `programs/output/`: generated plots and terminal captures.
- `programs/pdf/`: exported PDFs of the code.
- `report/`: LaTeX source and rendered report.
- `submission/`: final submitted report and code files.

## Run

Load the `queueing` package before running the scripts if your Octave installation does not load it automatically.

```bash
cd programs
octave lab2_part2.m
octave lab2_part3a.m
octave lab2_part3b.m
```
