# Premise notes (rubric item 1)

Working notes behind `report/sections/premise.tex`. This file holds the checks, rough
numbers and open questions that don't belong in the report text.

## 1. Motivation in plain language

Rocket turbopumps keep tank pressure low to save tank mass, so the propellant arrives at the
pump only a little above its vapor pressure. The inducer (the first, axial stage) runs with
some vapor on purpose. Near the blade leading edges and tip gap, pressure dips below vapor
pressure and small nuclei grow into bubbles. Those bubbles then move into higher pressure and
collapse hard: fast wall speeds, big local pressures. Near a surface that erodes blades; at
the pump level it costs head and can drive instabilities (rotating cavitation, cavitation
surge).

We study one spherical bubble with the Rayleigh-Plesset (RP) equation. It is the simplest
model that still captures a violent nonlinear collapse. Water at 20 C is a stand-in liquid.

## 2. Research questions

- **Q1 (core): collapse after a pressure jump.** How do R_min/R0, peak wall speed and peak gas
  pressure scale with pressure ratio Pi = p1/p0 and with R0? At what Pi does the gas stop
  mattering until the very end of the collapse?
- **Q2 (core): forced response.** Under p_inf = p0 (1 - eps sin(w t)), where is the resonance
  peak, how far does it move from w0 as eps grows, and at what eps does the response become
  violent?
- **Q3 (stretch): which physics wins.** Switch off gas / viscosity / surface tension one at a
  time across (R0, Pi) and build a regime map.

Changes from the suggested list: Q1 is framed as a pressure step because that both produces
dynamics and matches the "nucleus convected into recovery pressure" picture. Q3 is marked as a
stretch goal so it doesn't eat time from Q1 and Q2.

## 3. The IVP

```
R R'' + (3/2) R'^2 = (1/rho) [ p_v + p_g0 (R0/R)^(3 kappa) - p_inf(t) - 2 sigma/R - 4 mu R'/R ]
R(0) = R0,  R'(0) = 0
p_g0 = p_inf(0) - p_v + 2 sigma/R0      (Case A overrides: p_g0 = 0)
```

State y = [R, R'], so y' = [R', (bracket/rho - 1.5 R'^2) / R].

### Units check (each term times rho should be Pa)

| Term | Units | OK |
|---|---|---|
| rho R R'' | (kg/m^3)(m)(m/s^2) = kg/(m s^2) = Pa | yes |
| rho R'^2 | (kg/m^3)(m^2/s^2) = Pa | yes |
| p_v, p_inf, p_g0 (R0/R)^(3 kappa) | Pa (ratio is dimensionless) | yes |
| 2 sigma / R | (N/m)/m = N/m^2 = Pa | yes |
| 4 mu R' / R | (Pa s)(m/s)/m = Pa | yes |

### Parameters (water, 20 C)

| Symbol | Value | Confidence |
|---|---|---|
| rho | 998.2 kg/m^3 | high |
| mu | 1.002e-3 Pa s | high |
| sigma | 0.0728 N/m | high (IAPWS gives about 72.7 to 72.8 mN/m) |
| p_v | 2.34 kPa (2.339) | high |
| p0 | 101.325 kPa | exact (1 atm) |
| c (not in RP, used for validity check) | about 1482 m/s | high |
| kappa | 1.4 baseline, 1.0 comparison | **modeling choice, see flags** |

All of these are from memory and agree with standard tables to the digits given. **Someone
should still check them against the NIST WebBook and record the access date** in
`references.bib`.

Sweep ranges: R0 from 1 um to 1 mm; Pi from 1.5 to 20; eps from 0.01 to 0.8;
w/w0 from 0.2 to 3.

## 4. Test cases

### Case A: empty cavity (Rayleigh collapse)

p_g0 = p_v = sigma = mu = 0, p_inf = dp constant.

- First integral: R'^2 = (2 dp / 3 rho)(R0^3/R^3 - 1).
- Exact time to reach radius R: t(R) = R0 sqrt(3 rho / 2 dp) * integral from R/R0 to 1 of
  x^(3/2) / sqrt(1 - x^3) dx.
- Collapse time: t_c = (1/3) sqrt(3/2) B(5/6, 1/2) R0 sqrt(rho/dp).

**Constant check.** Evaluated with `math.gamma`: (1/3) sqrt(3/2) B(5/6, 1/2) = **0.914681**.
Brennen (Cavitation and Bubble Dynamics, ch. 2) and Rayleigh (1917) give 0.915. Matches.

**Rough numbers.** R0 = 1 mm, dp = 101.325 kPa, rho = 998.2:
sqrt(rho/dp) = sqrt(0.009852) = 0.09925 s/m, so t_R = 99.3 us and t_c = 0.9147 * 99.3 = 90.8 us.

**Independent cross-check (throwaway, not in repo).** A quick plain-Python RK4 run of
Case A with 20,000 fixed steps to t = 0.8255 t_R, the quadrature time for R = 0.5 R0, gave
R/R0 = 0.500, as expected. (The quadrature used the substitution x = 1 - s^2 to remove the
endpoint singularity; a naive midpoint rule was visibly off in the 3rd digit.) This only
confirms that the equation and the formula agree. The real check comes from our own
integrators in rubric item 3.

**Numerical caution.** R' goes to minus infinity at t_c, so the comparison must stop at some
small R (e.g. R/R0 = 0.05 or 0.01), not at t_c. Near the end R ~ (t_c - t)^(2/5) and the
local time scale R/|R'| = (5/2)(t_c - t).

### Case B: pressure step on a gas bubble

Equilibrium at p0, then p_inf = p1 = Pi p0 for t > 0. Initial acceleration
R''(0+) = (p0 - p1)/(rho R0) < 0.

Inviscid energy (mu = 0, kappa != 1) is conserved:
E = 2 pi rho R^3 R'^2 + (4 pi/3)(p1 - p_v) R^3 + 4 pi sigma R^2
  + (4 pi/3) p_g0 R0^(3 kappa) R^(3(1 - kappa)) / (kappa - 1).
Derivation: d/dt(2 pi rho R^3 R'^2) = 4 pi R^2 R' (p_B - p1), where p_B is the liquid
pressure at the wall, and each potential term's derivative cancels one piece of p_B. I
checked this by hand; please re-derive it once.

**Preliminary throwaway runs** (R0 = 100 um, kappa = 1.4, all physics on, fixed-step RK4
in plain Python, 300k steps over 3 t_R; t_R based on p1 - p_v). They confirm that the case
produces collapse and rebound, and give a rough sense of scale. **Do not quote them in the
report.** Regenerate them with our own code.

| Pi | R_min/R0 | t_min/t_R | max abs R' (m/s) | peak p_g (Pa) |
|---|---|---|---|---|
| 2 | 0.70 | 1.31 | 5.3 | 4.4e5 |
| 5 | 0.42 | 1.11 | 25 | 4.0e6 |
| 10 | 0.26 | 1.02 | 72 | 2.7e7 |

As Pi grows, t_min/t_R heads toward the Rayleigh value 0.915, which is a nice trend for Q1.
At Pi = 10, the wall Mach number is about 72/1482 = 0.05. At Pi = 20 it may approach 0.1,
which is where RP starts to be questionable (see flags).

### Case C: sinusoidal forcing

p_inf = p0 (1 - eps sin(w t)), starting from equilibrium.

Linearize with R = R0 (1 + x):
x'' + 2 beta x' + w0^2 x = (eps p0 / rho R0^2) sin(w t),
w0^2 = (3 kappa p_g0 - 2 sigma/R0) / (rho R0^2), beta = 2 mu / (rho R0^2).

| R0 | kappa | f0 | Q = w0/(2 beta) |
|---|---|---|---|
| 1 mm | 1.4 | 3.25 kHz | about 5,100 |
| 100 um | 1.4 | 32.7 kHz | about 510 |
| 10 um | 1.4 | 342 kHz | about 54 |
| 10 um | 1.0 | 288 kHz | about 45 |

Minnaert check (sigma = p_v = 0, kappa = 1.4, 1 atm): f0 R0 = 3.29 m/s. This matches the
commonly quoted figure of about 3.3 m/s (a 1 mm bubble rings near 3.3 kHz).

**Why R0 = 10 um for Case C.** With only viscous damping, a 1 mm bubble has Q of about 5000,
so the start-up transient takes thousands of cycles to die out before we can read off a
steady amplitude. At 10 um, Q is about 54: an e-fold decay takes about Q/pi = 17 cycles,
and settling to 1e-3 takes roughly 120 cycles. That's affordable.

## 5. Nondimensionalization

R~ = R/R0, t~ = t/t_R, t_R = R0 sqrt(rho/dp).

R~ R~'' + 1.5 R~'^2 = P_g R~^(-3 kappa) - P(t~) - S/R~ - (4/Re) R~'/R~

| Group | Definition | Meaning |
|---|---|---|
| P_g | p_g0/dp | gas cushion strength (Case B: about 1/Pi) |
| P(t~) | (p_inf(t) - p_v)/dp | forcing (1 for Cases A and B) |
| S | 2 sigma/(R0 dp) | surface tension vs driving pressure |
| Re | R0 sqrt(rho dp)/mu | inertia vs viscosity |
| kappa | | gas thermodynamics |
| eps, w t_R | | Case C forcing |

Check: (w0 t_R)^2 = 3 kappa P_g - S.

Crossovers in water at 1 atm: S = 1 at R0 = 1.47 um; Re = 1 at R0 = 0.10 um.

| R0 | S | Re |
|---|---|---|
| 1 mm | 0.0015 | 9,900 |
| 100 um | 0.015 | 990 |
| 10 um | 0.15 | 99 |

Choice of dp: p1 - p_v for Case B, p_inf for Case A, and p0 - p_v for Case C.

## 6. Why numerics

The full RP equation is nonlinear and has no closed form once gas, viscosity, surface tension
or time-varying p_inf are present. The only exact solution is Case A, and even that is
implicit (t as a function of R). The linearization only covers tiny oscillations, which is
exactly what the questions are not about. And the time scale collapses near R_min, which
motivates adaptive stepping.

## 7. Flags: things I'm unsure about or chose for you

1. **kappa = 1.4 constant.** Real bubbles go from near-isothermal (slow, small) to
   near-adiabatic (fast collapse). A fixed kappa is a standard simplification. Running both
   1.0 and 1.4 brackets it. Peak gas pressure and temperature are very sensitive to kappa, so
   present them as model-dependent.
2. **Viscous-only damping in Case C.** RP with only the 4 mu R'/R term leaves out thermal and
   acoustic radiation damping, which dominate for bubbles bigger than a few microns. Our
   resonance peaks will be sharper than in real water. We should either say so plainly or use
   an "effective viscosity". (I believe Brennen and Lauterborn & Kurz both discuss this, but
   I haven't checked the exact treatment.)
3. **Mach number limit.** RP assumes an incompressible liquid. I wrote "about 0.1" as the
   wall Mach number above which results are suspect. That's a rule of thumb, not a sourced
   number. Find a citation or soften it. Keller-Miksis is the usual compressible fix, but it's
   out of scope.
4. **Inertial-cavitation threshold for Q2.** A common criterion is R_max/R0 >= 2. I'm fairly
   sure this appears in the literature (associated with Flynn / Apfel), but I did not put it in
   the report or the bib because I can't confirm the source. Check before using it.
5. **Turbopump details.** I kept the turbopump discussion qualitative (inducer operates with
   some cavitation; erosion, head breakdown, rotating cavitation, cavitation surge). I did not
   include blade-passing frequencies or real inducer pressure ratios because I'd be guessing.
   If you want Pi or w tied to a real pump, find numbers from Brennen's *Hydrodynamics of
   Pumps* or a NASA SP (e.g. the inducer design monograph; I haven't cited it because I'm not
   sure of the number).
6. **Pressure step idealization.** A real nucleus sees a smooth pressure rise over a finite
   time. A step is the cleanest choice for a well-posed IVP with smooth RHS for t > 0. A
   smooth ramp could be a follow-up.
7. **Sweep ranges.** Pi up to 20 and eps up to 0.8 are my picks. Pi = 20 may push past the
   Mach limit. eps < 1 keeps p_inf positive.
8. **References.** `lauterborn2010` and `brennen1994pumps` have TODO notes in the bib file.
   The NIST entry needs an access date.
