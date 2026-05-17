# Lab 5 - Web-Server DES and Jackson Networks

## Goal

This assignment has two themes. Part A builds a discrete-event simulation of a single FIFO web server where every request passes through two service phases. Part B analyzes an open Jackson queueing network with five M/M/1 queues and probabilistic routing.

## Implementation

The implementation is split into two Octave entry points:

- `programs/lab5_part1.m`: simulates the web-server system. It logs arrivals, departures, queue sizes, response times, and throughput; tests stationarity for several arrival rates; computes confidence intervals from independent replications; removes an initial transient period guided by MQS(t); and verifies Little's Law.
- `programs/lab5_part2.m`: solves the Jackson network. It computes effective arrival rates, traffic intensities, ergodicity, mean number of clients per queue, end-to-end delay by Little's Law, the current bottleneck, the maximum admissible `lambda1`, and a delay curve as `lambda1` approaches the stability boundary.

Part A uses local random samplers and fallback quantile approximations, with optional use of the Octave `statistics` package when available.

## Inputs and Outputs

- Input material is in `assignment/`, including the lab statement and starter/help code.
- Generated web-server and Jackson-network plots are committed in `programs/output/`.
- Exported code and output PDFs are committed in `programs/pdf/`.

## Organization

- `assignment/`: lab statement, introduction material, and starter/help code.
- `programs/lab5_part1.m`: web-server discrete-event simulation.
- `programs/lab5_part2.m`: Jackson network analysis.
- `programs/output/`: generated figures, output logs, and output PDFs.
- `programs/pdf/`: exported PDFs of the code and output.
- `report/`: LaTeX source and rendered report.
- `submission/`: final submitted report and code files.

## Run

The scripts can use the Octave `statistics` package when available, but include fallback approximations for the required quantiles.

```bash
cd programs
octave --eval "lab5_part1"
octave --eval "lab5_part2"
```
