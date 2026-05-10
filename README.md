# Single-Train Energy-Optimal Driving Strategy Simulator
### Guangzhou Metro Line 7 — Multi-Algorithm Benchmark Suite

> MATLAB-based optimization framework for eco-driving on urban rail transit.  
> Compares Dynamic Programming, Convex Optimization, MILP, and multiple MOEAs  
> on real route and rolling stock data from Guangzhou Line 7.

---

## Overview

This project implements and benchmarks multiple optimization algorithms for the **energy-optimal train control problem** on an urban metro system. Given a fixed journey time constraint (scheduled time), the goal is to find a speed profile that minimizes traction energy consumption.

The framework supports:
- **8 inter-station segments** (IS01–IS08) of Guangzhou Metro Line 7
- **Both directions** (up/down)
- **Multiple solvers** ranging from exact methods (DP, CO, MILP) to population-based heuristics (NSGA-II, MOPSO, SPEA2, MOEAD, DE)
- **Experimental campaigns** for ablation, benchmarking, and robustness analysis

---

## Project Structure

```
new_single_train_simulator_cakra/
│
├── main/                          # Entry points & orchestrators
│   ├── main.m                     # Master runner (single or benchmark mode)
│   ├── run_single.m               # Single solver run
│   └── run_benchmark_compare.m    # Multi-solver multi-seed benchmark
│
├── route/                         # Guangzhou Line 7 route data
│   └── Guangxhou_line_7/
│       ├── up/                    # 8 segment .mat files (IS01–IS08, up direction)
│       └── down/                  # 8 segment .mat files (down direction)
│
├── rollingstocks/                 # Rolling stock specifications
│   ├── rollingstock_Guangzhou_L7.m
│   └── rollingstock_Guangzhou_L7_res_upd.m
│
├── driving_strategy/              # Physics simulator & decision variable encoding
│   ├── simulation_fun_CC_CR.m     # Main simulation wrapper
│   └── simulation_fun_CC_CR_base.m
│
├── optimiser/                     # Optimization algorithm implementations
│   ├── nsga2/                     # NSGA-II (4 variants)
│   ├── mopso_main.m               # Multi-Objective PSO
│   ├── spea2_main.m               # SPEA2
│   ├── moead_main.m               # MOEA/D
│   ├── de_moea_main.m             # Differential Evolution MOEA
│   ├── dp_main.m                  # Dynamic Programming (wrapper)
│   └── dqn_main.m                 # Deep Q-Network
│
├── Dynamic_Programming/           # Standalone DP implementations
│   ├── Basic/
│   │   └── TEST_DP_ONE_SEGMENT.m  # DP with Lagrangian relaxation (dx=1m)
│   │       └── results/           # Per-segment summary & .mat outputs
│   ├── COPILOT/
│   │   └── TEST_DP_MULTI_SEGMENT_COPILOT.m
│   └── Cloud_AI/
│       └── DP_IMPROVED_CAI.m      # Improved DP (Illinois λ-search, Gaussian smoothing)
│           └── results/
│
├── Convex Optimisation/           # Convex relaxation solver (CVX + SDPT3/MOSEK)
│   ├── main_CO.m
│   └── results/                   # Per-segment CO results
│
├── MILP/                          # Mixed Integer LP solver (YALMIP + MOSEK/Gurobi)
│   └── main_MILP.m
│
├── experiments/                   # Experimental campaigns (E1–E5)
│   ├── main_experiments.m         # Master orchestrator
│   ├── run_E1_ablation.m
│   ├── run_E2_tight_sweep.m
│   ├── run_E3_benchmarking.m
│   ├── run_E4_robustness.m
│   └── utils/                     # Metrics, plotting, Wilcoxon test helpers
│
├── setup/                         # Project & parallel pool initialization
├── dp_target_time_results/        # Pre-computed DP results (all 8 segments)
└── third_party/YALMIP/            # YALMIP modeling library
```

---

## Route & Rolling Stock

### Guangzhou Metro Line 7 — Inter-Station Segments

| Segment | Distance | Scheduled Time | Speed Limit |
|---------|----------|---------------|-------------|
| IS01 | 1,120 m | 129 s | 80 km/h |
| IS02 | 1,908 m | 169 s | 80 km/h |
| IS03 | 2,172 m | 184 s | 80 km/h |
| IS04 | 1,642 m | 177 s | 80 km/h |
| IS05 | 2,116 m | 185 s | 80 km/h |
| IS06 | 2,365 m | 219 s | 80 km/h |
| IS07 | 2,406 m | 214 s | 80 km/h |
| IS08 | 3,778 m | 329 s | 80 km/h |

### Rolling Stock (AW0 — Guangzhou Line 7)

| Parameter | Value |
|-----------|-------|
| Vehicle mass | 204 t |
| Inertial mass (with rotary allowance λ=0.08) | 220.3 t |
| Max speed | 80 km/h |
| Max traction power | 3,716.8 kW |
| Max braking power | 3,911.2 kW |
| Davis resistance | A=27 kN, B=0, C=0.0042 kN/(km/h)² |
| Max traction acceleration | ~1.42 m/s² |
| Max braking deceleration | ~1.72 m/s² |

---

## Algorithms Implemented

| Algorithm | Type | File | Notes |
|-----------|------|------|-------|
| **DP** (Lagrangian relaxation) | Exact | `Dynamic_Programming/Basic/TEST_DP_ONE_SEGMENT.m` | dx=1m, Illinois λ-search in CAI version |
| **DP Improved (CAI)** | Exact | `Dynamic_Programming/Cloud_AI/DP_IMPROVED_CAI.m` | dx=2m, Illinois, Gaussian smoothing |
| **Convex Optimization** | Exact (relaxed) | `Convex Optimisation/main_CO.m` | CVX + SDPT3/MOSEK |
| **MILP** | Exact (PWL approx) | `MILP/main_MILP.m` | YALMIP + MOSEK/Gurobi |
| **NSGA-II** (4 variants) | MOEA | `optimiser/nsga2/` | Original, RL+SDE, Improved, Bayes-RL |
| **MOPSO** | MOEA | `optimiser/mopso_main.m` | Multi-Objective PSO |
| **SPEA2** | MOEA | `optimiser/spea2_main.m` | Strength Pareto EA2 |
| **MOEA/D** | MOEA | `optimiser/moead_main.m` | Decomposition-based |
| **DE-MOEA** | MOEA | `optimiser/de_moea_main.m` | Differential Evolution |
| **DQN** | RL | `optimiser/dqn_main.m` | Deep Q-Network |

### DP Result Summary (Baseline, T = T_sched)

| Segment | E_trac (kWh) | T_actual (s) | Status |
|---------|-------------|--------------|--------|
| IS01 | 8.927 | 128.560 | matched |
| IS02 | 19.999 | 168.695 | matched |
| IS03 | 25.484 | 184.405 | matched |
| IS04 | 16.637 | 176.759 | matched |
| IS05 | 22.970 | 185.453 | matched |
| IS06 | 17.472 | 218.795 | matched |
| IS07 | 22.847 | 214.391 | matched |
| IS08 | 41.844 | 329.535 | closest |

---

## Requirements

### MATLAB
- MATLAB R2019b or later recommended
- Parallel Computing Toolbox (optional, for `parfor` acceleration)
- Statistics and Machine Learning Toolbox (for Wilcoxon test in experiments)

### Optimization Solvers
| Solver | Required for | How to get |
|--------|-------------|------------|
| **CVX** | Convex Optimization | [cvxr.com/cvx](http://cvxr.com/cvx) — free academic |
| **YALMIP** | MILP | Bundled in `third_party/YALMIP/` |
| **MOSEK** | CO + MILP (recommended) | [mosek.com](https://www.mosek.com) — free academic |
| **SDPT3** | CO (alternative) | Bundled with CVX |
| **Gurobi** | MILP (optional, fastest) | [gurobi.com](https://www.gurobi.com) — free academic |

> **Minimum setup:** CVX + SDPT3 (bundled) is sufficient to run all scripts.  
> MOSEK is strongly recommended for MILP and large CO problems.

---

## Quick Start

### 1. Clone / Open Project

```matlab
cd 'path\to\new_single_train_simulator_cakra'
```

### 2. Setup CVX (first time only)

```matlab
cd '<cvx-install-folder>'
cvx_setup
```

### 3. Run Dynamic Programming (all segments)

```matlab
% Basic DP — results saved to Dynamic_Programming/Basic/results/
run('Dynamic_Programming/Basic/TEST_DP_ONE_SEGMENT.m')

% Improved DP (Illinois + smoothing) — results to Dynamic_Programming/Cloud_AI/results/
run('Dynamic_Programming/Cloud_AI/DP_IMPROVED_CAI.m')
```

### 4. Run Convex Optimization

```matlab
% Edit main_CO.m to select segments and solver, then:
run('Convex Optimisation/main_CO.m')
```

### 5. Run MILP

```matlab
% Edit main_MILP.m to select segment and solver, then:
run('MILP/main_MILP.m')
```

### 6. Run MOEA Benchmarking Experiments

```matlab
% Full benchmark (E1–E4) for one segment:
ACTIVE_SEG = struct('name','IS04','file','Guangzhou_Line7_IS04_5.200-6.842km.mat','T_sched',177);
ACTIVE_NRUNS = 30;
run('experiments/main_experiments.m')
```

---

## Experimental Campaigns

| Study | Script | Description |
|-------|--------|-------------|
| **E1 — Ablation** | `run_E1_ablation.m` | Component contribution analysis (4 configs × 30 runs) |
| **E2 — Time Sweep** | `run_E2_tight_sweep.m` | Energy vs. time slack (8 slack levels × 30 runs) |
| **E3 — Benchmarking** | `run_E3_benchmarking.m` | 6 methods vs. DP baseline (30 runs, tight & loose slack) |
| **E4 — Robustness** | `run_E4_robustness.m` | 30 runs × 18 uncertainty scenarios |
| **E5 — LLM Advisor** | `run_E5_llm_advisor.m` | LLM-guided parameter tuning study |

---

## Key Configuration (main.m)

```matlab
% Select segment
ACTIVE_SEG = 'IS04';          % or 'IS01'..'IS08', or 'all'

% Select direction
DIRECTION  = 'up';            % 'up' | 'down'

% Select driving strategy
driving_strategy = 'CC_CR';   % 'CC' | 'CR' | 'CC_CR'

% Select optimizer
solver = 'nsga2';             % 'nsga2' | 'mopso' | 'spea2' | 'moead' | 'de_moea' | 'dp' | 'dqn'

% NSGA-II variant
nsga2_variant = 'rl_sde';     % 'original' | 'rl_sde' | 'improved' | 'bayes_rl'

% Parallel execution
parallel_use = true;
```

---

## DP Implementation Notes

The standalone DP (`Dynamic_Programming/`) uses **Lagrangian relaxation**:

```
minimize  J = E_traction (kWh) + λ · T (s)
```

By sweeping λ (energy-time trade-off weight), the solver finds speed profiles for multiple target travel times. The `DP_IMPROVED_CAI.m` version adds:

- **Illinois algorithm** for λ-search (superlinear convergence vs. bisection)
- **Gaussian post-processing** for smooth speed profiles (display only)
- **dx=2m spatial grid** for ~5× faster runtime with negligible accuracy loss
- **Regenerative braking estimate** reported as informational output (not used in optimization)

---

## Author

**Cakra Adipura Wicaksana**  
Universitas Sultan Ageng Tirtayasa (Untirta)  
[cakraadipura@untirta.ac.id](mailto:cakraadipura@untirta.ac.id)

---

## License

This project is for academic research purposes.  
Route data and rolling stock parameters are based on publicly available data from Guangzhou Metro Line 7.

---

*Generated with assistance from Claude Code (Anthropic) — claude-sonnet-4-6*
