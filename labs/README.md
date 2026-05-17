# Labs

This folder contains the six NTUA Queuing Systems lab assignments, their Octave/MATLAB implementations, generated artifacts, reports, and submission files.

## Lab Index

| Folder | Topic | Main entry points |
| --- | --- | --- |
| `1st/` | Probability distributions and Poisson processes | `programs/lab1_part1.m`, `programs/lab1_part2.m`, `programs/lab1_part3.m` |
| `2nd/` | M/M/1 and finite birth-death queues | `programs/lab2_part2.m`, `programs/lab2_part3a.m`, `programs/lab2_part3b.m` |
| `3rd/` | M/M/1/10 queue simulation | `programs/lab3.m` |
| `4th/` | Telephone-center sizing and Erlang C | `programs/lab4.m` |
| `5th/` | Web-server discrete-event simulation and Jackson networks | `programs/lab5_part1.m`, `programs/lab5_part2.m` |
| `6th/` | Time-series workload forecasting | `programs/lab6.m` |

## Common Structure

Most lab folders use this layout:

- `assignment/`: assignment PDFs, introduction material, starter code, and input data when provided.
- `programs/`: working Octave/MATLAB implementations.
- `programs/output/`: generated figures, terminal output, CSV files, and output PDFs.
- `programs/pdf/`: exported PDFs of code and selected output.
- `report/`: LaTeX report source and rendered report PDF.
- `submission/`: final submitted code and report material.

## Running Scripts

Run each script from its own `programs/` directory unless the lab README says otherwise:

```bash
cd 1st/programs
octave lab1_part1.m

cd ../../6th/programs
octave --eval "lab6"
```

The lab-level README files contain the complete run commands and package requirements for each assignment.
