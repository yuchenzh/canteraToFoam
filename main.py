import os
from canteraToFoam import canteraToFoam


ROOT = os.path.dirname(os.path.abspath(__file__))


def main(mech: str = "single", thermo_only: bool = False, make_tests: bool = False):
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

    if make_tests:
        test_dir = os.path.join(ROOT, "test")
        os.makedirs(test_dir, exist_ok=True)

        single_base = os.path.join(ROOT, "singleStep")
        single_yaml = os.path.join(single_base, "cantera", "chem.yaml")
        single_ctf = canteraToFoam(single_yaml)
        t_single = single_ctf.read_transport(os.path.join(single_base, "chemkin", "transportProperties"))
        thermo_single = single_ctf.thermo_file_string(transport=t_single)
        with open(os.path.join(test_dir, "single_thermo"), "w") as f:
            f.write(thermo_single)
        with open(os.path.join(test_dir, "single_reactions"), "w") as f:
            f.write(single_ctf.chemistry_file_string())

        multi_base = os.path.join(ROOT, "multiStep")
        multi_yaml = os.path.join(multi_base, "cantera", "grimech30.yaml")
        multi_ctf = canteraToFoam(multi_yaml)
        t_multi = multi_ctf.read_transport(os.path.join(multi_base, "chemkin", "transportProperties"))
        thermo_multi = multi_ctf.thermo_file_string(transport=t_multi)
        with open(os.path.join(test_dir, "multi_thermo"), "w") as f:
            f.write(thermo_multi)

        lines = []
        for i, rxn in enumerate(multi_ctf.gas.reactions()):
            rtype = multi_ctf.check_reaction_type(i)
            lines.append(f"{i}: {rtype} - {rxn.equation}")
        with open(os.path.join(test_dir, "multi_reaction_types"), "w") as f:
            f.write("\n".join(lines))

    print("\nReaction types:")
    for i, rxn in enumerate(ctf.gas.reactions()):
        rtype = ctf.check_reaction_type(i)
        print(f"  {i}: {rtype} - {rxn.equation}")

    # Show one example of a standard reversible Arrhenius reaction block
    for i in range(len(ctf.gas.reactions())):
        if ctf.check_reaction_type(i) == "reversibleArrhenius":
            example = ctf._reversible_reaction_string(i)
            print("\nExample reversibleArrhenius block:\n")
            print(example)
            break


if __name__ == "__main__":
    import sys
    mech = "single"
    thermo_only = False
    if len(sys.argv) > 1:
        mech = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == "thermo":
        thermo_only = True
    main(mech, thermo_only, make_tests=True)
