"""Gas-phase standard-state Gibbs free energy for H2 + D2 -> 2 HD.

Electronic structure is B3LYP/6-31G(d).  The electronic potential energy
surface is identical for these isotopologues, so one H2 geometry optimization
and Hessian are evaluated.  Isotopologue-specific masses are then used in the
rigid-rotor/harmonic-oscillator thermochemistry at 298.15 K and 1 atm.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from pyscf import dft, gto
from pyscf.data import nist
from pyscf.geomopt.geometric_solver import optimize
from pyscf.hessian import thermo

TEMPERATURE_K = 298.15
PRESSURE_PA = 101325.0
H_MASS_U = 1.00782503223
D_MASS_U = 2.01410177812
HARTREE_TO_KJMOL = nist.HARTREE2J * nist.AVOGADRO / 1000.0


def gibbs_from_hessian(mol, hessian, masses_u: np.ndarray, symmetry_number: int):
    """Return RRHO thermal quantities using explicit isotope masses."""
    vib = thermo.harmonic_analysis(mol, hessian, mass=masses_u)
    frequencies_au = np.asarray(vib["freq_au"]).real
    positive = frequencies_au[frequencies_au > 0.0]

    coordinates_m = mol.atom_coords() * nist.BOHR_SI
    masses_kg = masses_u * nist.ATOMIC_MASS
    center_of_mass = np.einsum("i,ij->j", masses_kg, coordinates_m) / masses_kg.sum()
    centered = coordinates_m - center_of_mass
    inertia = np.einsum("i,ij,ik->jk", masses_kg, centered, centered)
    inertia = np.eye(3) * np.trace(inertia) - inertia
    principal_moments = np.sort(np.linalg.eigvalsh(inertia))
    moment_perpendicular = principal_moments[-1]

    boltzmann = nist.BOLTZMANN
    planck = nist.PLANCK
    gas_constant_eh_per_k = boltzmann / nist.HARTREE2J
    e0 = float(mol_energy(mol))
    total_mass_kg = masses_kg.sum()
    q_trans = (
        (2.0 * math.pi * total_mass_kg * boltzmann * TEMPERATURE_K / planck**2) ** 1.5
        * boltzmann * TEMPERATURE_K / PRESSURE_PA
    )
    entropy_trans = gas_constant_eh_per_k * (2.5 + math.log(q_trans))
    enthalpy_trans = 2.5 * gas_constant_eh_per_k * TEMPERATURE_K

    rotational_constant_hz = planck / (8.0 * math.pi**2 * moment_perpendicular)
    q_rot = boltzmann * TEMPERATURE_K / (symmetry_number * planck * rotational_constant_hz)
    entropy_rot = gas_constant_eh_per_k * (1.0 + math.log(q_rot))
    enthalpy_rot = gas_constant_eh_per_k * TEMPERATURE_K

    au_to_hz = math.sqrt(nist.HARTREE2J / (nist.ATOMIC_MASS * nist.BOHR_SI**2)) / (2.0 * math.pi)
    vib_temperatures = positive * au_to_hz * planck / boltzmann
    reduced_temperatures = vib_temperatures / TEMPERATURE_K
    exp_term = np.exp(-reduced_temperatures)
    zpe = 0.5 * gas_constant_eh_per_k * vib_temperatures.sum()
    entropy_vib = gas_constant_eh_per_k * np.sum(
        reduced_temperatures * exp_term / (1.0 - exp_term) - np.log(1.0 - exp_term)
    )
    enthalpy_vib = zpe + gas_constant_eh_per_k * TEMPERATURE_K * np.sum(
        reduced_temperatures * exp_term / (1.0 - exp_term)
    )

    gibbs_components = {
        "electronic": e0,
        "translational": enthalpy_trans - TEMPERATURE_K * entropy_trans,
        "rotational": enthalpy_rot - TEMPERATURE_K * entropy_rot,
        "vibrational": enthalpy_vib - TEMPERATURE_K * entropy_vib,
    }
    gibbs = sum(gibbs_components.values())
    return {
        "electronic_energy_hartree": e0,
        "gibbs_free_energy_hartree": float(gibbs),
        "gibbs_free_energy_kj_mol": float(gibbs * HARTREE_TO_KJMOL),
        "zpe_hartree": float(zpe),
        "gibbs_components_hartree": {
            name: float(value) for name, value in gibbs_components.items()
        },
        "frequency_cm-1": [float(x) for x in np.asarray(vib["freq_wavenumber"]).real],
        "symmetry_number": symmetry_number,
        "masses_u": [float(x) for x in masses_u],
    }


def mol_energy(mol) -> float:
    mf = dft.RKS(mol)
    mf.xc = "B3LYP"
    return mf.kernel()


def main() -> None:
    initial = gto.M(
        atom="H 0 0 0; H 0 0 0.7414",
        basis="6-31G(d)",
        charge=0,
        spin=0,
        unit="Angstrom",
        verbose=4,
    )
    mf = dft.RKS(initial)
    mf.xc = "B3LYP"
    optimized = optimize(mf)
    final_mf = dft.RKS(optimized)
    final_mf.xc = "B3LYP"
    final_mf.kernel()
    hessian = final_mf.Hessian().kernel()

    # thermochemistry calls mol_energy once per isotopologue; pass the same
    # optimized geometry and electronic energy through a lightweight closure.
    global mol_energy
    electronic_energy = float(final_mf.e_tot)
    original_mol_energy = mol_energy
    mol_energy = lambda _mol: electronic_energy
    try:
        species = {
            "H2": gibbs_from_hessian(optimized, hessian, np.array([H_MASS_U, H_MASS_U]), 2),
            "D2": gibbs_from_hessian(optimized, hessian, np.array([D_MASS_U, D_MASS_U]), 2),
            "HD": gibbs_from_hessian(optimized, hessian, np.array([H_MASS_U, D_MASS_U]), 1),
        }
    finally:
        mol_energy = original_mol_energy

    delta_g_hartree = 2.0 * species["HD"]["gibbs_free_energy_hartree"] - species["H2"]["gibbs_free_energy_hartree"] - species["D2"]["gibbs_free_energy_hartree"]
    delta_g_components_hartree = {
        component: (
            2.0 * species["HD"]["gibbs_components_hartree"][component]
            - species["H2"]["gibbs_components_hartree"][component]
            - species["D2"]["gibbs_components_hartree"][component]
        )
        for component in ("electronic", "translational", "rotational", "vibrational")
    }
    results = {
        "reaction": "H2 + D2 -> 2 HD",
        "method": "B3LYP/6-31G(d)",
        "temperature_K": TEMPERATURE_K,
        "pressure_Pa": PRESSURE_PA,
        "optimized_HH_distance_angstrom": float(np.linalg.norm(optimized.atom_coords()[1] - optimized.atom_coords()[0]) * nist.BOHR),
        "species": species,
        "delta_g_hartree": float(delta_g_hartree),
        "delta_g_kj_mol": float(delta_g_hartree * HARTREE_TO_KJMOL),
        "delta_g_components_hartree": {
            name: float(value) for name, value in delta_g_components_hartree.items()
        },
        "delta_g_components_kj_mol": {
            name: float(value * HARTREE_TO_KJMOL)
            for name, value in delta_g_components_hartree.items()
        },
    }
    Path("h2_d2_hd_gibbs_results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
