---
task: unsupported
engine: pyscf
error_class: support_gap
outcome: workaround
n_atoms: 2
method: B3LYP
basis: 6-31G(d)
cores: 1
mode: slurm
---
Symptom: MAESTRO has no task for isotope-substituted Gibbs thermochemistry; it provides isotope-shifted frequencies only.

Attempts: Ran a direct PySCF B3LYP/6-31G(d) geometry optimization and Hessian on Slurm job 14098 with one CPU on node4. Applied explicit H and D isotope masses in the rigid-rotor/harmonic-oscillator thermochemistry calculation, including the homonuclear symmetry number of 2 for H2 and D2 and 1 for HD.

Result: Completed successfully. For H2 + D2 -> 2 HD at 298.15 K and 1 atm, delta G is -2.919 kJ mol-1. The decomposed contributions are electronic 0.000, translational -0.437, rotational -3.145, and vibrational +0.663 kJ mol-1.

Context: Gas-phase isotope-exchange reaction. Scheduler accounting was unavailable, so no measured wall-time value is recorded.
