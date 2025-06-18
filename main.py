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

    # Show one example of a standard reversible Arrhenius reaction block
    for i in range(len(ctf.gas.reactions())):
        if ctf.check_reaction_type(i) == "reversibleArrhenius":
            example = ctf._reversible_reaction_string(i)
            print("\nExample reversibleArrhenius block:\n")
            print(example)
            break


def run_tests(output_dir: str = "test_outputs") -> None:
    """Generate thermo and reaction files for both mechanisms.

    Parameters
    ----------
    output_dir
        Directory where the generated files are written.

    The following files are created inside *output_dir*:
    ``single_thermo.foam``      thermo data for the single-step mechanism
    ``single_reactions.foam``   reaction blocks for the single-step mechanism
    ``multi_thermo.foam``       thermo data for the multi-step mechanism
    ``multi_reaction_types.txt`` list of reaction types for the multi-step mechanism
    """

    os.makedirs(output_dir, exist_ok=True)

    print("*** Single step mechanism ***")
    single_base = os.path.join(ROOT, "singleStep")
    s_yaml = os.path.join(single_base, "cantera", "chem.yaml")
    ctf_single = canteraToFoam(s_yaml)
    s_transport = ctf_single.read_transport(
        os.path.join(single_base, "chemkin", "transportProperties")
    )
    s_thermo = ctf_single.thermo_file_string(transport=s_transport)
    s_react = ctf_single.chemistry_file_string()
    print("\nThermo file:\n")
    print(s_thermo)
    print("\nReactions file:\n")
    print(s_react)
    with open(os.path.join(output_dir, "single_thermo.foam"), "w") as f:
        f.write(s_thermo)
    with open(os.path.join(output_dir, "single_reactions.foam"), "w") as f:
        f.write(s_react)

    print("\n*** Multi-step mechanism ***")
    multi_base = os.path.join(ROOT, "multiStep")
    m_yaml = os.path.join(multi_base, "cantera", "grimech30.yaml")
    ctf_multi = canteraToFoam(m_yaml)
    m_transport = ctf_multi.read_transport(
        os.path.join(multi_base, "chemkin", "transportProperties")
    )
    m_thermo = ctf_multi.thermo_file_string(transport=m_transport)
    print("\nThermo file:\n")
    print(m_thermo)
    with open(os.path.join(output_dir, "multi_thermo.foam"), "w") as f:
        f.write(m_thermo)
    print("\nReaction types:")
    rtypes = []
    for i in range(len(ctf_multi.gas.reactions())):
        rtype = ctf_multi.check_reaction_type(i)
        rtypes.append(f"{i + 1}: {rtype}")
        print(f"{i + 1}: {rtype}")
    with open(os.path.join(output_dir, "multi_reaction_types.txt"), "w") as f:
        f.write("\n".join(rtypes))


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        run_tests()
    else:
        mech = "single"
        thermo_only = False
        if len(sys.argv) > 1:
            mech = sys.argv[1]
        if len(sys.argv) > 2 and sys.argv[2] == "thermo":
            thermo_only = True
        main(mech, thermo_only)
