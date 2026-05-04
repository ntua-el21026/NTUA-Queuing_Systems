# Lab 1 - Probability Distributions and Poisson Processes

## Goal

This assignment introduces the probability distributions used throughout queueing theory. It studies the Poisson distribution, the binomial-to-Poisson approximation, exponential service/interarrival times, the memoryless property, and the Poisson counting process.

## Implementation

The Octave programs are split into three parts:

- `programs/lab1_part1.m`: plots Poisson PMFs for several rates, compares binomial distributions with the limiting Poisson distribution, computes mean and variance values, and studies Poisson splitting/superposition for external and internal calls.
- `programs/lab1_part2.m`: plots exponential PDFs and CDFs, computes the distribution of the minimum of independent exponential variables, and verifies the memoryless property numerically.
- `programs/lab1_part3.m`: simulates a Poisson counting process from exponential interarrival times, estimates the event rate for increasing sample sizes, and discusses event counts in non-overlapping time windows.

The scripts use the Octave `statistics` package.

## Organization

- `assignment/`: lab statement, introduction material, Octave tutorial, and starter/help code.
- `programs/`: implemented scripts.
- `programs/output/`: saved plots and terminal output for the three parts.
- `programs/pdf/`: exported PDFs of the code.
- `report/`: LaTeX source and rendered report.
- `submission/`: final submitted report and code files.

## Run

```bash
cd programs
octave lab1_part1.m
octave lab1_part2.m
octave lab1_part3.m
```
