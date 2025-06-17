#!/usr/bin/env python3
"""Simple demonstration of the :class:`canteraToFoam` class."""

from canteraToFoam import canteraToFoam


def main() -> None:
    ctf = canteraToFoam("./cantera/chem.yaml")

    # Show Arrhenius parameters for the first reaction
    A, Ta, b = ctf.get_arrhenius_parameters(0)
    print(f"Reaction 0 parameters: A={A}, Ta={Ta}, b={b}")

    # Forward and backward rate constants for the same reaction
    k_f, k_r = ctf.get_forward_backward_rate(0)
    print(f"Forward rate constant: {k_f}")
    print(f"Reverse rate constant: {k_r}")


if __name__ == "__main__":
    main()

