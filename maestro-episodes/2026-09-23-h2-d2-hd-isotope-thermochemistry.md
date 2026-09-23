---
task: unsupported
engine: none
error_class: support_gap
outcome: support_gap
---
Symptom: The requested reaction Gibbs free energy for H2 + D2 -> 2 HD requires thermochemistry with deuterium isotope masses. ThermoTask accepts molecular geometries but does not expose isotope-mass substitutions; a D element label is rejected when computing charge recommendations.

Attempts: Searched MAESTRO capabilities for free energy, thermochemistry, Gibbs, isotope, and deuterium. ThermoTask produces gibbs_free_energy for standard element masses. IsotopeShiftTask produces isotope_frequencies only.

Result: No MAESTRO task produces isotope-substituted Gibbs free energies, so the requested calculation is outside current MAESTRO support.

Context: Gas-phase H2 + D2 -> 2 HD at 298.15 K and 1 atm; user requested one CPU on an idle Slurm node.
