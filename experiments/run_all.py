"""Regenerate every figure in figures/ (a few minutes on a laptop)."""

import runpy
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
SCRIPTS = [
    "fig_verification_caseA.py",
    "fig_convergence_caseB.py",
    "fig_verification_caseC.py",
    "fig_stability_stepsize.py",
    "fig_tolerance.py",
    "q1_pressure_sweep.py",
    "q2_resonance.py",
    "q3_physics.py",
]

if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(HERE))
    for name in SCRIPTS:
        print(f"\n=== {name} ===")
        runpy.run_path(str(HERE / name), run_name="__main__")
