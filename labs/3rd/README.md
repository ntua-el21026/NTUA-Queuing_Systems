# Lab 3 - M/M/1/10 Simulation

## Goal

This assignment focuses on simulation accuracy for a finite-capacity queue. The model is an M/M/1/10 system with arrival rates `lambda = 1, 5, 10`, service rate `mu = 5`, and finite capacity `K = 10`. The goal is to verify the simulator, estimate steady-state quantities, and examine the effect of removing an initial transient period.

## Implementation

The main implementation is `programs/lab3.m`. It:

- Prints a debug trace for the first 30 transitions to validate the event logic.
- Runs simulations for several `(T, B)` pairs, where `T` is the number of transitions and `B` is the transient cutoff.
- Tracks time spent in each state, rejection probability, mean number in the system, sojourn times, and convergence of the estimated mean population.
- Produces grouped bar charts of state probabilities with and without transient removal.
- Produces convergence plots for the estimated mean number in the system.

The simulator is implemented directly in Octave using random event selection and exponential event times.

## Inputs and Outputs

- Input material is in `assignment/`, including the lab statement and starter/help code.
- Generated state-probability and convergence plots are committed in `programs/output/`.
- Exported code and output PDFs are committed in `programs/pdf/`.

## Organization

- `assignment/`: lab statement, introduction material, and starter/help code.
- `programs/lab3.m`: simulation implementation.
- `programs/output/`: generated figures and output logs for the configured run lengths.
- `programs/pdf/`: exported PDFs of the code and output.
- `report/`: LaTeX source and rendered report.
- `submission/`: final submitted report and code file.

## Run

The implementation does not require an Octave package beyond the base runtime.

```bash
cd programs
octave --eval "lab3"
```
