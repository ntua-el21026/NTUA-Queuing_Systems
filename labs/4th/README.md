# Lab 4 - Telephone Center and Erlang C

## Goal

This assignment studies the sizing of telephone-center resources. The first part models an IP PBX with 200 employees, 50 outgoing trunks, Poisson call arrivals, exponential call durations, and an infinite waiting queue. The second part computes Erlang C service metrics for a multi-agent call center.

## Implementation

The main implementation is `programs/lab4.m`. It:

- Computes offered load from call arrival rate and mean holding time.
- Checks the Erlang C stability condition.
- Computes `P0`, probability of waiting, service level at a selected waiting threshold, average speed of answer, and occupancy.
- Cross-checks the waiting probability from the Octave `erlangc` function against the closed-form expression.

The script uses the Octave `queueing` package.

## Organization

- `assignment/`: lab statement, extra material, introduction material, and starter/help code.
- `programs/lab4.m`: Erlang C calculation script.
- `programs/output/`: terminal output and exported output PDF.
- `programs/pdf/`: exported PDFs of the code and output.
- `report/`: LaTeX source and rendered report.
- `submission/`: final submitted report and code file.

## Run

```bash
cd programs
octave lab4.m
```
