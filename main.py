import os
from canteraToFoam import canteraToFoam


ROOT = os.path.dirname(os.path.abspath(__file__))


def main(mech: str = "single", thermo_only: bool = False):
    """Run a small demonstration using the requested mechanism.

    Parameters
    ----------
    mech
        Either "single" or "multi" to select the mechanism.
    thermo_only
        If True, print only the OpenFOAM thermo file.
    """

    base = os.path.join(ROOT, "singleStep" if mech == "single" else "multiStep")
    yaml_name = "chem.yaml" if mech == "single" else "grimech30.yaml"
    yaml_path = os.path.join(base, "cantera", yaml_name)
    try:
        ctf = canteraToFoam(yaml_path)
    except ImportError as exc:
        print(f"Cannot run example: {exc}")
        return

    t_data = ctf.read_transport(os.path.join(base, "chemkin", "transportProperties"))
    thermo_str = ctf.thermo_file_string(transport=t_data)
    if thermo_only:
        print(thermo_str)
        return

    A, Ta, b = ctf.get_arrhenius_parameters(0)
    print(f"First reaction Arrhenius parameters:\n  A = {A}\n  Ta = {Ta}\n  b = {b}")

    orders = ctf.reaction_orders(0)
    print("Forward orders:", orders)

    print("\nGenerated OpenFOAM thermo file:\n")
    print(thermo_str)

    chem_str = ctf.chemistry_file_string()
    print("\nGenerated OpenFOAM chemistry file:\n")
    print(chem_str)


if __name__ == "__main__":
    import sys
    mech = "single"
    thermo_only = False
    if len(sys.argv) > 1:
        mech = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == "thermo":
        thermo_only = True
    main(mech, thermo_only)
