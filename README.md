# Marginal Distributionally Robust Kalman Filter (MDRKF)

This repository contains the implementation of the Marginal Distributionally Robust Kalman Filter (MDRKF), a robust estimation technique for systems with uncertain parameters. The MDRKF algorithm combines the principles of Kalman Filtering with Distributionally Robust Optimization to provide enhanced estimation performance under uncertainty.

## Installation

To get started with the MDRKF project, clone this repository to your local machine:

```bash
git clone https://github.com/lylechan42/MDRKF.git
cd MDRKF
```

## Usage

The primary scripts for the MDRKF algorithm are:

- `dynamic_estimation_MC.m`: Implements the Monte Carlo Simulation to evaluate our methods for dynamic estimation.
- `dynamic_estimation_test.m`: Test script to run our MDRKF algorithm briefly.

## Utilities

The `utils` directory contains various utility functions used in the MDRKF algorithm:

- `DRO.m`: Distributionally Robust Optimization implementation.
- `Frank_Wolfe.m`: Frank-Wolfe algorithm for convex optimization.
- `KF.m`: Standard Kalman Filter implementation.
- `MDRO.m`: Marginal Distributionally Robust Optimization implementation.
- `genSigma.m`: Function to generate covariance matrices.
- `getCorrData.m`: Function to retrieve trajectory with correlated measurement.
- `getData.m`: Function to retrieve trajectory with uncorrelated measurement.
- `getRanData.m`: Function to retrieve trajectory with random parameters.
- `getSys.m`: Function to retrieve system parameters.
- `isPSD.m`: Function to check if a matrix is positive semi-definite.
- `my_dist.m`: Custom distance function.
- `reconstruct.m`: Function to reconstruct estimation matrices.
- `save_estimations.m`: Function to save estimation results.
- `tau_update.m`: Function to update estimation with tau divergence based DRO.

## License

This project is licensed under the MIT License - see the [LICENSE](https://chat.openai.com/c/LICENSE) file for details.
