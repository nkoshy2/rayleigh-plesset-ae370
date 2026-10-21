# Experiments

One script per report figure. Each script imports the integrators from `src/rp/`,
writes its figure to `../figures/`, and can be run on its own:

```
python experiments/<script>.py
```

Planned (names are placeholders until the figures are settled):

| Script | Figure | Rubric item |
|---|---|---|
| `fig_caseA_rayleigh_check.py` | Case A numerical vs exact R(t) and collapse time | 3 |
| `fig_convergence_caseB.py` | Error vs step size (RK4) and vs tolerance/cost (adaptive) | 3 |
| `fig_step_size_history.py` | Adaptive step size through collapse and rebound | 2, 3 |
| `fig_collapse_scaling.py` | R_min/R0, peak wall speed, peak gas pressure vs pressure ratio and R0 | 4 |
| `fig_frequency_response.py` | Steady response amplitude vs forcing frequency for several amplitudes | 4 |
| `fig_regime_map.py` | Which term (gas, viscosity, surface tension) controls the collapse | 4 |
