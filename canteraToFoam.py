import os
import sys
import importlib

# Ensure we import the actual Cantera package if installed
current_dir = os.path.dirname(os.path.abspath(__file__))
restored = False
if current_dir in sys.path:
    sys.path.remove(current_dir)
    restored = True
try:
    ct = importlib.import_module("cantera")
except ModuleNotFoundError:
    ct = None
finally:
    if restored:
        sys.path.insert(0, current_dir)


class canteraToFoam:
    """Wrapper for extracting reaction info from a Cantera YAML file."""

    def __init__(self, yaml_path: str):
        if ct is None or not hasattr(ct, "Solution"):
            raise ImportError("Cantera is required to use canteraToFoam")
        self.yaml_path = yaml_path
        self.gas = ct.Solution(yaml_path)

    def get_arrhenius_parameters(self, index: int):
        """Return Arrhenius A, Ta, and b for the reaction at *index*."""
        rxn = self.gas.reactions()[index]
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

    def reaction_orders(self, index: int):
        """Return the forward reaction orders (FORD entries) as a dict."""
        rxn = self.gas.reactions()[index]
        return getattr(rxn, "orders", {})

    def _set_state(self, T: float, concentrations: dict):
        """Internal helper to set the gas state from species concentrations."""
        conc_list = [concentrations.get(sp, 0.0) for sp in self.gas.species_names]
        total = sum(conc_list)
        mole_fractions = [0.0] * len(conc_list)
        P = ct.one_atm
        if total > 0:
            mole_fractions = [c / total for c in conc_list]
            P = total * ct.gas_constant * T
        self.gas.TPX = T, P, mole_fractions

    def forward_rate(self, index: int, T: float, composition: dict):
        """Return the forward rate of progress for the given state."""
        self._set_state(T, composition)
        return self.gas.forward_rates_of_progress[index]

    def backward_rate(self, index: int, T: float, composition: dict):
        """Return the backward rate of progress for the given state."""
        self._set_state(T, composition)
        rxn = self.gas.reactions()[index]
        if rxn.reversible:
            return self.gas.reverse_rates_of_progress[index]
        return 0.0

    # --- Transport reader ---------------------------------------------------

    def read_transport(self, path: str) -> dict:
        """Read an OpenFOAM-style ``transportProperties`` file."""
        data = {}
        current = None
        if not os.path.exists(path):
            return data
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("//"):
                    continue
                if line.endswith("{") and not line.startswith("transport"):
                    current = line.strip("{").strip().strip('"')
                    data.setdefault(current, {})
                elif line.startswith("As"):
                    try:
                        data[current]["As"] = float(line.split()[1].strip(";"))
                    except (IndexError, ValueError, TypeError):
                        pass
                elif line.startswith("Ts"):
                    try:
                        data[current]["Ts"] = float(line.split()[1].strip(";"))
                    except (IndexError, ValueError, TypeError):
                        pass
        return data

    # --- OpenFOAM thermo file helpers -------------------------------------------------

    def _species_thermo_string(
        self,
        name: str,
        As: float = 1.67212e-06,
        Ts: float = 170.672,
    ) -> str:
        """Return the OpenFOAM thermo entry for ``name``."""

        sp = self.gas.species(name)
        thermo = sp.thermo
        t_range = thermo.input_data["temperature-ranges"]
        low_coeffs, high_coeffs = thermo.input_data["data"]

        lines = [f"{name}", "{"]
        lines += [
            "    specie",
            "    {",
            f"        molWeight       {sp.molecular_weight};",
            "    }",
            "    thermodynamics",
            "    {",
            f"        Tlow            {t_range[0]};",
            f"        Thigh           {t_range[2]};",
            f"        Tcommon         {t_range[1]};",
            "        highCpCoeffs    ( "
            + " ".join(str(c) for c in high_coeffs)
            + " );",
            "        lowCpCoeffs     ( "
            + " ".join(str(c) for c in low_coeffs)
            + " );",
            "    }",
            "    transport",
            "    {",
            f"        As              {As};",
            f"        Ts              {Ts};",
            "    }",
            "    elements",
            "    {",
        ]
        for el, amt in sp.composition.items():
            lines.append(f"        {el}               {int(amt)};")
        lines += ["    }", "}"]
        return "\n".join(lines)

    def thermo_file_string(self, species_list=None, transport=None) -> str:
        """Return a string representing an OpenFOAM ``thermos`` file."""

        if species_list is None:
            species_list = self.gas.species_names

        lines = [
            f"species         {len(species_list)} ( {' '.join(species_list)} );",
        ]

        for name in species_list:
            As = 1.67212e-06
            Ts = 170.672
            if transport and name in transport:
                As = transport[name].get("As", As)
                Ts = transport[name].get("Ts", Ts)
            elif transport and "*" in transport:
                As = transport["*"].get("As", As)
                Ts = transport["*"].get("Ts", Ts)
            lines.append("")
            lines.append(self._species_thermo_string(name, As, Ts))

        return "\n".join(lines)

    # --- Reaction writers ---------------------------------------------------

    def check_reaction_type(self, index: int) -> str:
        """Return the OpenFOAM reaction type for the reaction at *index*."""
        rxn = self.gas.reactions()[index]
        if getattr(rxn, "orders", None):
            return "irreversibleArrhenius"
        return "arrhenius"

    def _arrhenius_coeffs(self, index: int):
        return self.get_arrhenius_parameters(index)

    def _reaction_equation_string(self, rxn) -> str:
        left = []
        orders = getattr(rxn, "orders", {})
        for sp, nu in rxn.reactants.items():
            coeff = f"{nu}" if nu != 1 else ""
            order = orders.get(sp)
            if order is not None:
                left.append(f"{coeff}{sp}^{order}")
            else:
                left.append(f"{coeff}{sp}")
        right = []
        for sp, nu in rxn.products.items():
            coeff = f"{nu}" if nu != 1 else ""
            right.append(f"{coeff}{sp}")
        return " + ".join(left) + " = " + " + ".join(right)

    def _irreversible_reaction_string(self, index: int) -> str:
        rxn = self.gas.reactions()[index]
        A, Ta, b = self._arrhenius_coeffs(index)
        eq = self._reaction_equation_string(rxn)
        lines = [
            f"un-named-reaction-{index}",
            "{",
            "    type            irreversibleArrhenius;",
            f"    reaction        \"{eq}\";",
            f"    A               {A};",
            f"    beta            {b};",
            f"    Ta              {Ta};",
            "}",
        ]
        return "\n".join(lines)

    def chemistry_file_string(self, indices=None) -> str:
        lines = ["reactions", "{"]
        if indices is None:
            indices = range(len(self.gas.reactions()))
        for i in indices:
            typ = self.check_reaction_type(i)
            if typ == "irreversibleArrhenius":
                block = self._irreversible_reaction_string(i)
            else:
                block = ""  # unhandled types
            for ln in block.splitlines():
                lines.append("    " + ln)
        lines.append("}")
        lines.append("Tlow            200;")
        lines.append("Thigh           6000;")
        return "\n".join(lines)
