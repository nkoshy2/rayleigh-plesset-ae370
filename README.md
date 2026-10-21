# Cavitation Bubble Collapse: Rayleigh-Plesset Dynamics with Adaptive Runge-Kutta

AE 370 (Aerospace Numerical Methods), Fall 2026, Group Project 1.

## What this is

Pumps that feed rocket engines spin fast enough that the pressure on their blades can drop
below the liquid's vapor pressure. Tiny bubbles form there, get carried into higher pressure,
and collapse hard. Those collapses erode blades and can throttle pump performance.

We model a single spherical bubble with the Rayleigh-Plesset equation, a nonlinear ODE for the
bubble radius R(t), and solve it with time integrators we write ourselves. A fixed-step RK4
scheme is the baseline. An adaptive embedded Runge-Kutta pair is the main method, since the
collapse squeezes the relevant time scale by orders of magnitude right when it matters most.
We use the solver to ask how hard a bubble collapses as a function of its size and the
pressure jump it sees, how it responds to oscillating pressure, and which physics controls the
collapse in which regime.

The full premise is in [`report/sections/premise.tex`](report/sections/premise.tex), with
working notes in [`docs/premise-notes.md`](docs/premise-notes.md).

## Repository layout

```
src/rp/
  model.py         Dimensionless Rayleigh-Plesset RHS, Jacobian, test cases A/B/C,
                   exact and reference solutions (Rayleigh time, energy, linear response)
  integrators.py   From-scratch RK4, fixed-step DP5, adaptive Dormand-Prince 5(4)
                   with FSAL, step-size control and event location
  diagnostics.py   Collapse quantities: R_min, peak wall speed, peak gas pressure, Mach
tests/             Unit tests (pytest)
experiments/       One script per report figure, plus run_all.py
figures/           Figure PDFs produced by experiments/ and used in the report
report/            LaTeX report: main.tex, references.bib, sections/
docs/              Project statement, premise notes, AI usage log, contributions
```

No time integrator from an outside library (for example `scipy.integrate.solve_ivp`) is used
to advance the solution. That is a project requirement.

## Install

```
python -m venv .venv
.venv\Scripts\activate          # Windows; use `source .venv/bin/activate` on macOS/Linux
pip install -r requirements.txt
pip install -e .
```

## Run

```
pytest                              # 24 tests, a few seconds
python experiments/run_all.py       # regenerate every figure in figures/ (about 2 minutes)
```

A minimal example:

```python
import numpy as np
from rp.model import case_B, rhs
from rp.integrators import dp5

prm = case_B(R0=100e-6, Pi=5.0)            # bubble hit by a 5x pressure step
sol = dp5(lambda t, y: rhs(t, y, prm), (0.0, 3.0), [1.0, 0.0],
          rtol=1e-8, atol=1e-10, positive=0)
print(sol.y[:, 0].min(), sol.nfev)         # R_min / R0 and cost
```

To build the report (needs a LaTeX install, for example TeX Live or MiKTeX):

```
cd report
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

To upload to Overleaf, run `python report/make_overleaf_zip.py` and upload the
`report-overleaf.zip` it writes (sources plus figure PDFs).

## Team

| Name | NetID |
|---|---|
| Naveen Koshy | nkoshy2 |
| Ava Hoppe | arhoppe2 |
| Alison Nazaruk | nazaruk3 |

Roles and contribution estimates: [`docs/contributions.md`](docs/contributions.md).
AI usage: [`docs/ai-usage-log.md`](docs/ai-usage-log.md).

## License

MIT, see [`LICENSE`](LICENSE).
