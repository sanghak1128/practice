# H2 + D2 -> 2 HD Gibbs free energy

## Result

At 298.15 K and 1 atm in the gas phase,

\[
\Delta G = 2G(\mathrm{HD}) - G(\mathrm{H_2}) - G(\mathrm{D_2})
= -0.00111185\ \mathrm{E_h}
= \mathbf{-2.919\ kJ\ mol^{-1}}.
\]

The negative value indicates that the equilibrium is thermodynamically
favorable toward HD under the stated model.

## Reaction Gibbs-energy decomposition

The values below use the same stoichiometric combination,
\(2\mathrm{HD} - \mathrm{H_2} - \mathrm{D_2}\). The electronic term is
zero because all three isotopologues have the same electronic Hamiltonian
and optimized electronic energy at this level of theory.

| Contribution | Delta G (Eh) | Delta G (kJ mol^-1) |
| --- | ---: | ---: |
| Electronic | 0.0000000000 | 0.000 |
| Translational | -0.0001664508 | -0.437 |
| Rotational | -0.0011979505 | -3.145 |
| Vibrational | +0.0002525508 | +0.663 |
| **Total** | **-0.0011118505** | **-2.919** |

Thus, the favorable free-energy change is dominated by the rotational term;
the vibrational term partly offsets it.

## Method

- Electronic structure: B3LYP/6-31G(d), closed-shell RKS.
- Geometry: H2 optimized to an H-H distance of 0.742788 A.
- Thermochemistry: rigid-rotor/harmonic-oscillator treatment with explicit
  isotope masses (H = 1.00782503223 u; D = 2.01410177812 u).
- Standard state: ideal gas, 298.15 K and 1 atm (101325 Pa).
- Scheduler: Slurm job 14072, requested on idle `node4` with one CPU.

The electronic potential energy surface is the same for the three
isotopologues, so one optimized H2 geometry and Hessian were used; masses and
rotational symmetry numbers were varied in the thermochemistry step.

## Calculated values

| Species | G (Eh) | G (kJ mol^-1) | Harmonic frequency (cm^-1) |
| --- | ---: | ---: | ---: |
| H2 | -1.1768247244 | -3089.7529 | 4453.12 |
| D2 | -1.1814276915 | -3101.8380 | 3150.04 |
| HD | -1.1796821332 | -3097.2550 | 3857.01 |

The geometry optimization converged and each isotopologue has one positive
vibrational mode. Full machine-readable values are in
`h2_d2_hd_gibbs_results.json`.

## Files

- `h2_d2_hd_gibbs.py`: PySCF calculation and isotope-aware thermochemistry.
- `run_h2_d2_hd_gibbs.sbatch`: one-CPU Slurm submission script.
- `h2_d2_hd_gibbs_results.json`: calculation output.

MAESTRO supports ordinary molecular Gibbs thermochemistry but not
isotope-substituted Gibbs thermochemistry; the capability gap is recorded in
`maestro-episodes/2026-09-23-h2-d2-hd-isotope-thermochemistry.md`.
