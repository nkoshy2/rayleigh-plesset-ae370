# Experiments

One script per report figure. Each imports the model and integrators from `src/rp/`
and writes a PDF (used by the report) and a PNG preview to `../figures/`. Run from the
repository root after `pip install -e .`:

```
python experiments/run_all.py              # every figure, a few minutes
python experiments/fig_tolerance.py        # or one at a time
```

Each script also prints the numbers quoted in the report.

| Script | Figure(s) | What it shows | Rubric |
|---|---|---|---|
| `fig_verification_caseA.py` | `fig_verification_caseA` | Convergence to the exact Rayleigh solution; collapse time to R = 1e-3 R0 | 3(i) |
| `fig_convergence_caseB.py` | `fig_convergence_caseB` | Orders 4 and 5 on the full equation; energy conservation (inviscid) | 3(i) |
| `fig_verification_caseC.py` | `fig_verification_caseC` | Small-amplitude forced response vs linear theory | 3(i) |
| `fig_stability_stepsize.py` | `fig_stability`, `fig_stepsize` | h*lambda along real trajectories vs stability limits; adaptive step history | 2, 3(ii) |
| `fig_tolerance.py` | `fig_tolerance` | Error of R_min, peak wall speed, peak gas pressure vs rtol | 3(ii) |
| `q1_pressure_sweep.py` | `q1_pressure_sweep` | R_min and peak gas pressure vs pressure ratio for three R0 (Table 2) | 4 (Q1) |
| `q2_resonance.py` | `q2_resonance` | Steady forced response vs frequency for three forcing amplitudes | 4 (Q2) |
| `q3_physics.py` | `q3_physics` | R_min with viscosity or surface tension removed, vs R0 | 4 (Q3) |

`common.py` holds the shared plot style (reference categorical palette in fixed order;
a one-hue blue ramp for ordered sweeps such as the pressure ratio).
