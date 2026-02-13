# Standard-Model-Fitting

This repository contains MATLAB implementations for fitting the **Standard Model of diffusion in white matter**. Exemplified here on Monte Carlo simulation data, specifically addressing intra-axonal water components. Fitting code also include **Standard Model of diffusion, phase and relaxation in white matter** accounting for susceptibility-induced relaxation (R2 and R2*) and orientation-dependent phase shifts.

## 🔬 Project Overview
The project uses a **Separable Least Squares** approach. It leverages a Lebedev grid for spherical integration and is optimized for **GPU acceleration** using MATLAB's `gpuArray` and parallelized kernel operations.

### Key Features:
* **SMPR Fitting:** Standard Model Parameter Regression.
* **ISO Fitting:** Isotropic-constrained fitting.
* **GPU Optimized:** High-performance Jacobian calculations using persistent variables.
* **Dynamic Pathing:** Location-independent script execution (handles the `Run/` folder context automatically).

