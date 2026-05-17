# NTUA Queuing Systems

Coursework repository for the NTUA ECE Queuing Systems course. It contains the assignment statements, Octave/MATLAB implementations, generated figures and logs, LaTeX reports, and final submission material for Labs 1-6.

## Repository Layout

```text
docs/
  lectures/             Course lecture PDFs.
  project_structure/    Helper script and text snapshot of the repository tree.

labs/
  README.md             Index of the lab folders and their common structure.
  1st/                  Probability distributions and Poisson processes.
  2nd/                  M/M/1 and finite birth-death queue analysis.
  3rd/                  M/M/1/K simulation with transient removal.
  4th/                  Multi-server telephone center and Erlang C sizing.
  5th/                  Discrete-event web-server simulation and Jackson networks.
  6th/                  Time-series workload forecasting.
```

Each lab folder follows the same organization:

- `assignment/`: original assignment PDF, introduction PDF, and any starter/help code.
- `programs/`: working Octave/MATLAB source files, generated output text, figures, and exported code PDFs.
- `report/`: LaTeX report source and rendered PDF.
- `submission/`: final submitted report and code files.

Some labs include extra inputs, such as `dataset.mat` in Lab 6.

## Lab Summary

| Lab | Main topic | Implementation focus |
| --- | --- | --- |
| 1 | Poisson, binomial, exponential distributions, and Poisson counting processes | Uses Octave statistics functions to plot PMFs/PDFs/CDFs, verify moments and memorylessness, and simulate event arrivals. |
| 2 | Analytical and numerical M/M/1 and M/M/1/K models | Uses the Octave `queueing` package for steady-state and transient CTMC calculations. |
| 3 | Simulation of an M/M/1/10 queue | Implements an event-driven simulator, compares results before/after transient removal, and studies convergence for multiple run lengths. |
| 4 | Telephone center dimensioning | Models an M/M/50 telephone system and computes Erlang C metrics such as offered load, waiting probability, service level, and ASA. |
| 5 | Web-server DES and open Jackson queueing network | Implements a two-phase FIFO web-server simulator, confidence intervals, Little's Law checks, and analytical Jackson network calculations. |
| 6 | Time-series workload forecasting | Compares linear-regression and seasonal-difference filter models for GPU request forecasting using prediction intervals and median CRPS. |

## Requirements

- GNU Octave.
- Octave `statistics` package for Lab 1, parts of Lab 5, and Lab 6.
- Octave `queueing` package for Labs 2 and 4.
- XeLaTeX for rebuilding the reports, because the LaTeX files use `fontspec`, `unicode-math`, and Greek text.

Example package loading is already included in the scripts where needed:

```octave
pkg load statistics
pkg load queueing
```

## Running the Code

Run scripts from their own `programs/` directory so relative output paths and generated figures are easy to manage. Examples:

```bash
cd labs/1st/programs
octave lab1_part1.m
octave lab1_part2.m
octave lab1_part3.m

cd ../../2nd/programs
octave lab2_part2.m
octave lab2_part3a.m
octave lab2_part3b.m

cd ../../3rd/programs
octave --eval "lab3"

cd ../../6th/programs
octave --eval "lab6"
```

Generated figures, terminal captures, CSV files, and exported PDFs are already stored under each lab's `programs/output/` and `programs/pdf/` directories where applicable. Lab 6 writes regenerated figures and CSV metrics to the current directory, so run it from `labs/6th/programs` and move regenerated artifacts into `programs/output/` when refreshing the committed outputs.

## Reports

Each lab report lives in `labs/<lab>/report/`. The rendered PDF is committed together with `main.tex`. To rebuild a report manually:

```bash
cd labs/1st/report
xelatex main.tex
```

Some reports include generated figures from the corresponding `programs/output/` directory, so copy or symlink figures into the report build directory if rebuilding from scratch.
