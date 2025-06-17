#!/usr/bin/env python3
"""Utility class for accessing reaction data from a Cantera mechanism."""

import os
import sys
import importlib

# Remove the repository directory from ``sys.path`` so that importing
# ``cantera`` loads the actual library rather than the local ``cantera``
# folder shipped with this repository.
current_dir = os.path.dirname(os.path.abspath(__file__))
restored = False
if current_dir in sys.path:
    sys.path.remove(current_dir)
    restored = True

try:
    ct = importlib.import_module("cantera")
except ModuleNotFoundError as exc:
    raise RuntimeError("Cantera library not found") from exc
finally:
    if restored:
        sys.path.insert(0, current_dir)


class canteraToFoam:
    """Read data from a Cantera YAML/CTI file.

    Parameters
    ----------
    yaml_path : str
        Path to the Cantera format file to load.
    """

    def __init__(self, yaml_path: str) -> None:
        self.yaml_path = yaml_path
        self.gas = ct.Solution(yaml_path)

    def get_arrhenius_parameters(self, index: int):
        """Return Arrhenius parameters ``A``, ``Ta`` and ``b`` for a reaction.

        The ``index`` is zero-based.
        """
        rxns = self.gas.reactions()
        if index < 0 or index >= len(rxns):
            raise IndexError("Reaction index out of range")
        rxn = rxns[index]
        rate = rxn.rate
        try:
            A = rate.pre_exponential_factor
            b = rate.temperature_exponent
            Ea = rate.activation_energy
        except AttributeError:
            A = rate.A
            b = rate.b
            Ea = rate.Ea
        Ta = Ea / ct.gas_constant
        return A, Ta, b

    def get_forward_backward_rate(self, index: int):
        """Return forward and backward rate constants for a reaction.

        For irreversible reactions the backward rate constant is returned as
        ``0.0``.
        """
        if index < 0 or index >= self.gas.n_reactions:
            raise IndexError("Reaction index out of range")
        k_f = self.gas.forward_rate_constants[index]
        rxn = self.gas.reaction(index)
        if rxn.reversible:
            k_r = self.gas.reverse_rate_constants[index]
        else:
            k_r = 0.0
        return k_f, k_r

