#!/usr/bin/env python3
"""Extract reactions from a Cantera YAML file.

This script loads the mechanism using Cantera and prints the reaction
string along with Arrhenius parameters A, Ta, and b.

Usage:
    python extract_reactions.py [PATH_TO_YAML]

If no path is given, ``cantera/chem.yaml`` is used.
"""

import os
import sys
import importlib

# Attempt to import the Cantera library. The repository contains a folder
# named ``cantera`` which can shadow the actual Cantera installation. To
# avoid this, remove the local path from ``sys.path`` before importing.
current_dir = os.path.dirname(os.path.abspath(__file__))
local_cantera = os.path.join(current_dir, "cantera")

# Temporarily remove the repository path so that importing ``cantera`` pulls
# in the actual Cantera library instead of the local ``cantera`` directory.
restored = False
if current_dir in sys.path:
    sys.path.remove(current_dir)
    restored = True

try:
    ct = importlib.import_module("cantera")
except ModuleNotFoundError as exc:
    sys.exit(f"Cantera library not found: {exc}")
finally:
    if restored:
        sys.path.insert(0, current_dir)


def main(yaml_path: str) -> None:
    """Print reaction equations and Arrhenius parameters."""
    gas = ct.Solution(yaml_path)
    for idx, rxn in enumerate(gas.reactions()):
        print(f"Reaction {idx + 1}: {rxn.equation}")
        rate = rxn.rate
        # Newer Cantera versions expose Arrhenius parameters via attributes
        # ``pre_exponential_factor``, ``temperature_exponent`` and
        # ``activation_energy``. Fall back to older attribute names if needed.
        try:
            A = rate.pre_exponential_factor
            b = rate.temperature_exponent
            Ea = rate.activation_energy
        except AttributeError:
            A = rate.A
            b = rate.b
            Ea = rate.Ea
        Ta = Ea / ct.gas_constant
        print(f"  A  = {A}")
        print(f"  Ta = {Ta}")
        print(f"  b  = {b}\n")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join("cantera", "chem.yaml")
    main(path)
