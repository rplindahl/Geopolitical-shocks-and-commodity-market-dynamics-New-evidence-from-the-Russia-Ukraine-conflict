# Geopolitical-shocks-and-commodity-market-dynamics-New-evidence-from-the-Russia-Ukraine-conflict

**This project utilizes code adjusted from the baseline of [Boer et al. (2023)](https://onlinelibrary.wiley.com/doi/epdf/10.1002/jae.2952?getft_integrator=sciencedirect_contenthosting&src=getftr&utm_source=sciencedirect_contenthosting)**

This repository contains the MATLAB codebase and data for estimating the impact of geopolitical shocks—specifically the Russian invasion of Ukraine—on agricultural, energy, and financial markets using an event-based Structural Vector Autoregression (SVAR) framework.

## Structure
* `main.m`: The primary execution script. It loads the data, specifies the SVAR model, performs partial identification, runs the moving block bootstrap (MBB) for confidence intervals, and outputs the Impulse Response Functions (IRFs) and summary statistics.
* `functions/`: A directory containing all necessary custom MATLAB functions for VAR estimation, bootstrap sampling, SVAR identification, and hypothesis testing.

*Note: The code uses `ann_base` as the default event-date input. To reproduce specific disaggregated analyses (e.g., Figure 7), change the event-date input to `ann_china` or `ann_attack`.*

## Methodology Highlights
* **Identification:** Employs heteroskedasticity-based partial identification to isolate the geopolitical shock.
* **Bootstrapping:** Uses Moving Block Bootstrap (MBB) with 50 draws by default to generate 68% and 95% confidence intervals for the IRFs.
* **Testing:** Includes automated descriptive statistics, Augmented Dickey-Fuller (ADF), and Phillips-Perron (PP) stationarity tests for all variables.

## How to Run
1. Clone or download this repository.
2. Open MATLAB and navigate to the repository's root directory.
3. Ensure you have the **Econometrics Toolbox** and **Optimization Toolbox** installed.
4. Ensure you're able to access the data as described in the paper. (The authors do not have permission to share the financial data from LSEG Workspace)
5. Open `main.m`.
6. *(Optional)* Adjust the `nrep` variable (default is `50`) if you want to speed up the bootstrap process for testing, or increase it for final paper-quality results. 
7. Run `main.m`. The script will automatically add the `functions/` folder to your path, execute the model, and generate the IRF plots and statistical tables.
