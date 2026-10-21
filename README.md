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
src/rp/          Model (RHS) and from-scratch integrators (RK4, adaptive embedded RK)
tests/           Unit tests (pytest)
experiments/     One script per report figure
figures/         Figures produced by experiments/ and used in the report
report/          LaTeX report: main.tex, references.bib, sections/
docs/            Project statement, premise notes, AI usage log, contributions
```

No time integrator from an outside library (for example `scipy.integrate.solve_ivp`) is used
to advance the solution. That is a project requirement.

## Install

Placeholder until the code exists. Planned:

```
python -m venv .venv
.venv\Scripts\activate          # Windows; use `source .venv/bin/activate` on macOS/Linux
pip install -r requirements.txt
pip install -e .
```

## Run

Placeholder. Planned:

```
pytest                                    # tests
python experiments/fig_convergence_caseB.py   # one figure, written to figures/
```

To build the report (needs a LaTeX install, for example TeX Live or MiKTeX):

```
cd report
latexmk -pdf main.tex
```

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
