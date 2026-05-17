# Lab 6 - Time-Series Workload Forecasting

## Goal

This assignment studies time-series forecasting for GPU resource-request data. The first 1200 hourly observations are used as training data, and the remaining observations are forecast with 95% prediction intervals. The implementation compares linear-regression forecasting models against filter-based seasonal-difference models.

## Implementation

The main implementation is `programs/lab6.m`. It:

- Loads `dataset.mat` and validates that it contains the `values` vector.
- Uses the first 1200 samples as training data and forecasts the remaining 47 samples in the committed dataset.
- Evaluates five linear-regression models with trend and seasonal terms for 12-hour, 24-hour, and 168-hour periods.
- Evaluates five filter models based on seasonal differencing with lags 1, 24, and 168.
- Selects the best model in each family using median CRPS.
- Reports MAE, RMSE, prediction-interval coverage, and median CRPS.
- Writes model-comparison metrics to `lab6_model_comparison_results.csv`.
- Generates forecast, residual, error, and interval-comparison figures.

For the committed run, the selected regression model is `LR3: linear trend + daily 24 + weekly 168`, and the selected filter model is `F3: Delta_168`. The regression model has the lower median CRPS in the final comparison.

The script uses the Octave `statistics` package for normal-distribution functions.

## Inputs and Outputs

- `assignment/dataset.mat`: original dataset distributed with the assignment.
- `programs/dataset.mat`: copy used by the executable script.
- `programs/output/`: generated figures, terminal output, CSV metrics, and output PDF from the committed run.
- `programs/pdf/`: exported PDFs of the code and output.

## Organization

- `assignment/`: lab statement, introduction material, forecasting starter code, and dataset.
- `programs/lab6.m`: forecasting and model-selection implementation.
- `programs/output/`: generated figures, comparison CSV, and output logs.
- `programs/pdf/`: exported PDFs of the code and output.
- `report/`: LaTeX source and rendered report.
- `submission/`: final submitted code file.

## Run

Run from `programs/` so `dataset.mat` is in the current directory:

```bash
cd programs
octave --eval "lab6"
```

The script writes regenerated PNG and CSV artifacts to the current directory. Move refreshed artifacts into `programs/output/` if you want to replace the committed run outputs.
