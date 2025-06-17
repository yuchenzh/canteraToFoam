from canteraToFoam import canteraToFoam


def main():
    yaml_path = "cantera/chem.yaml"
    try:
        ctf = canteraToFoam(yaml_path)
    except ImportError as exc:
        print(f"Cannot run example: {exc}")
        return

    A, Ta, b = ctf.get_arrhenius_parameters(0)
    print(f"First reaction Arrhenius parameters:\n  A = {A}\n  Ta = {Ta}\n  b = {b}")

    orders = ctf.reaction_orders(0)
    print("Forward orders:", orders)

    t_data = ctf.read_transport("chemkin/transportProperties")
    thermo_str = ctf.thermo_file_string(transport=t_data)
    print("\nGenerated OpenFOAM thermo file:\n")
    print(thermo_str)

    chem_str = ctf.chemistry_file_string()
    print("\nGenerated OpenFOAM chemistry file:\n")
    print(chem_str)


if __name__ == "__main__":
    main()
