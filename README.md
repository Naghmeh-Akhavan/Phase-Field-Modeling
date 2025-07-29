# Phase Field Model (PFM) Simulation

This repository contains the official implementation of the paper:

**"Phase-Field Modeling of Border Cell Cluster Migration in Drosophila"**  
Naghmeh Akhavan, et al., 2025

We provide a MATLAB-based simulation framework for modeling collective cell migration in the Drosophila egg chamber using phase-field methods. The model captures interactions among nurse cells, the oocyte, and the migrating border cell cluster, incorporating chemoattractant signaling, adhesion, and interfacial tension. A key feature is the Tangential Interface Migration (TIM) force, which enables contact-mediated, tangential movement along nurse cell surfaces. The framework supports dynamic simulations with visualization tools for analyzing migration behavior and chemo-mechanical interactions.

## Quick Start

1. **Configure Parameters**: Modify `config.m` to set up your simulation parameters
2. **Run Simulation**: Execute `go_solver.m` to run the phase field model
3. **Visualize Results**: Use `go_show_results.m` to generate visualizations

## Project Structure

```
pfm_github/
├── README.md                    # This file
├── go_solver.m                  # Main solver script
├── go_show_results.m            # Visualization script
├── config.m                     # Configuration parameters
├── solve_pfm.m                  # Core phase field solver (steady state)
├── solve_pfm_concen.m           # Concentration dynamics solver
├── results/                     # Output data files (.mat)
└── visuals/                     # Generated visualizations (videos, plots)
```

## Configuration (`config.m`)

The `config.m` file contains all simulation parameters organized into several categories:

### Biological Parameters

- **Nurse Cells**: 6 cells with configurable radii, positions, and target volumes (red cells)
- **Oocyte**: Single large cell (blue cell)
- **Cluster Cell**: Mobile cell that can respond to $F_\text{chem}$ and $F_\text{TIM}$ (green cell)

```matlab
% Cell type I (Nurse cells)
params.cell_radi_nurse = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005];
params.x_centers_nurse = [1.5, 2.3, 3.1, 1.4, 2.3, 3.0];
params.y_centers_nurse = [1.6, 1.9, 2., 3.3, 3.1, 3];
params.target_volumes_nurse = [0.5, 0.6, 1., 0.6, 0.8, 1.1];
% Cell type II (Oocyte)
params.cell_radi_oocyte = 0.0005;
params.x_centers_oocyte = 3.4;
params.y_centers_oocyte = 2.5;
params.target_volumes_oocyte = 1.95;
% Cell type III (Cluster)
params.cell_radi_cluster = 0.0005;
params.x_centers_cluster = .95;
params.y_centers_cluster = 2.5;
params.target_volumes_cluster = 0.2;
```

#### Physical Parameters
- **Diffusion coefficients** (`params.dval`): Controls thickness of interface
- **Adhesion parameters** (`params.eta`): Cell-cell adhesion strength
- **Repulsion parameters** (`params.beta`): Cell-cell repulsion strength
- **Epithelial interactions** (`params.beta_s`, `params.eta_s`): wall/epithelial effects

### Discretization Parameters

```matlab
params.x_length = 5;           % Domain length
params.y_length = 5;           % Domain width  
params.h_space = 0.05;         % Spatial mesh size
params.h_time = 0.05;          % Time step size
params.t_final = 500;          % Final simulation time
```

### Model Setup

These boolean flags control which components of the simulation are executed:

```matlab
params.run_steadyQ = 1;        % Run steady-state solver
params.run_concenQ = 1;        % Run with concentration dynamics
params.run_tensionQ = 1;       % Include TIM forces
```

- **`1`**: Enable/run the specified component
- **`0`**: Disable/skip the specified component

**Component Descriptions:**
- **`run_steadyQ`**: Controls execution of the steady-state phase field solver that establishes initial cell configurations
- **`run_concenQ`**: Enables concentration dynamics (requires steady-state solution)
- **`run_tensionQ`**: Includes Tangential Interface Migration (TIM) forces for contact-mediated cell movement



## Main Solver (`go_solver.m`)

The main solver script orchestrates the simulation in two phases:

### Phase 1: Steady State Solution
The steady state solution only needs to be computed once for a given set of biological parameters. Once completed, the results are automatically saved in the `results/` folder with a unique filename based on the parameter configuration (`param_id_string`). For subsequent runs with the same parameters, you can set `params.run_steadyQ = 0` to skip this phase and load the existing steady state solution.

### Phase 2: Concentration Dynamicsd

This phase builds upon the steady state solution to incorporate concentration dynamics and chemotaxis. The solver uses the final cell configuration from Phase 1 as the initial condition. To enable or disable concentration dynamics, modify `params.run_concenQ` in `config.m`:

- **`params.run_concenQ = 1`**: Run simulation with concentration dynamics, chemical gradients, and chemotactic cell migration
- **`params.run_concenQ = 0`**: Skip concentration dynamics (only steady state cell configuration will be available)

#### Concentration Calculation Method
The `params.fixed_ConcenQ` parameter controls how the concentration is calculated:

- **`params.fixed_ConcenQ = 0`**: Uses a fixed radius for concentration calculations
- **`params.fixed_ConcenQ = 1`**: Calculates concentration based on cross-sectional area by computing the dynamic radius `r` 

## Visualization (`go_show_results.m`)

You can generate different types of visual outputs by modifying the parameters at the top of `go_show_results.m`. 

### Visualization Types


1. **3D Surface Videos** (`type_visual = 1`): Three-dimensional cell surfaces over time
2. **Contour Boundary Videos** (`type_visual = 2`): Cell boundaries as contour lines
3. **Filled Contour Videos** (`type_visual = 3`): Filled contour representations
   

### Data Types to Visualize
- **Cells** (`type_f = 1`): Phase field variables for all cell types
- **Concentration** (`type_f = 2`): Chemical concentration fields
- **TIM force** (`type_f = 3`): TIM force $F_\text{TIM}$
- **Chemoattractant** (`type_f = 4`): Chemical force $F_\text{chem}$

Figure 3 (type_visual =1, type_f=1), Figure 4 (type_visual =2, type_f=1), Figure 7 (type_visual =2, type_f=1), Figure 8 (type_visual =5). 

### Other plots
1. **Residuals Plot** (`type_visual = 4`): Convergence analysis / Relative change of the cluster phase
2. **Cluster Dynamics** (`type_visual = 5`): Position and velocity


## Author

**Naghmeh Akhavan** - July 2025

---

*For questions or issues, please refer to the comments within individual MATLAB files or contact the author Naghmeh Akhavan (nakhava1@umbc.edu) or Bradford E. Peercy (bpeercy@umbc.edu).*