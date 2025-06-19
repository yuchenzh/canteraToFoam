import os
import sys
import importlib
import numpy as np
from copy import deepcopy

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
        self.third_body_dict = self._construct_dummy_third_body_dict()
        
    def _construct_dummy_third_body_dict(self):
        returnDict = {}
        for sp in self.gas.species_names:
            returnDict[sp] = 1
        return returnDict

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
        """Read an OpenFOAM-style ``transportProperties`` file.

        The returned dictionary contains an entry for every species in the
        mechanism.  Species explicitly listed in ``transportProperties`` retain
        their specified ``As`` and ``Ts`` values while all others take the
        values from the ``".*"`` entry.
        """

        raw = {}
        current = None

        if not os.path.exists(path):
            return {}

        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("//"):
                    continue

                # species name lines
                if line.startswith('"'):
                    # handle both '"H2"' and '"H2"{' forms
                    current = line.strip('"')
                    if current.endswith('{'):
                        current = current[:-1].strip()
                    raw.setdefault(current, {})
                    continue

                if line == '{' or line == 'transport' or line.startswith('transport'):
                    # ignore standalone braces and transport keywords
                    continue

                if line.startswith("As") and current is not None:
                    try:
                        raw[current]["As"] = float(line.split()[1].strip(";"))
                    except (IndexError, ValueError, TypeError):
                        pass
                    continue

                if line.startswith("Ts") and current is not None:
                    try:
                        raw[current]["Ts"] = float(line.split()[1].strip(";"))
                    except (IndexError, ValueError, TypeError):
                        pass
                    continue

                # reset current after closing a species block
                if line == "}":
                    current = None

        default = raw.get(".*") or raw.get("*") or {"As": 1.67212e-06, "Ts": 170.672}

        data = {}
        for sp in self.gas.species_names:
            if sp in raw:
                data[sp] = {
                    "As": raw[sp].get("As", default.get("As")),
                    "Ts": raw[sp].get("Ts", default.get("Ts")),
                }
            else:
                data[sp] = default.copy()

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

    def is_reversible(self, index: int) -> bool:
        """Return ``True`` if the reaction is reversible."""

        rxn = self.gas.reactions()[index]
        return getattr(rxn, "reversible", False)

    def is_falloff(self, index: int) -> bool:
        """Return ``True`` if the reaction is a falloff reaction."""

        rtype = getattr(self.gas.reactions()[index], "reaction_type", "").lower()
        return rtype in {"falloff-lindemann", "falloff-troe"}

    def is_troe(self, index: int) -> bool:
        """Return ``True`` if the reaction uses Troe falloff."""

        rtype = getattr(self.gas.reactions()[index], "reaction_type", "").lower()
        return rtype == "falloff-troe"

    def is_lindemann(self, index: int) -> bool:
        """Return ``True`` if the reaction uses the Lindemann form."""

        rtype = getattr(self.gas.reactions()[index], "reaction_type", "").lower()
        return rtype == "falloff-lindemann"

    def check_reaction_type(self, index: int) -> str:
        """Return the OpenFOAM reaction type for the reaction at *index*.

        The mapping is based on :attr:`Reaction.reaction_type` provided by
        Cantera. Only a few common types are recognized.
        """

        rxn = self.gas.reactions()[index]
        rtype = getattr(rxn, "reaction_type", "").lower()

        # Reactions with explicit orders are treated as irreversible Arrhenius
        # since OpenFOAM cannot represent arbitrary forward orders.
        if getattr(rxn, "orders", None):
            return "irreversibleArrhenius"

        if self.is_plog(index):
            return "pressureDependentArrhenius"

        if rtype == "three-body-arrhenius":
            return (
                "reversibleThirdBodyArrhenius"
                if self.is_reversible(index)
                else "irreversibleThirdBodyArrhenius"
            )

        if not self.is_reversible(index):
            return "irreversibleArrhenius"

        if not self.is_falloff(index):
            return "reversibleArrhenius"

        if self.is_troe(index):
            return (
                "reversibleArrheniusTroeFallOff"
                if self.is_reversible(index)
                else "irreversibleArrheniusTroeFallOff"
            )

        if self.is_lindemann(index):
            return (
                "reversibleArrheniusLindemannFallOff"
                if self.is_reversible(index)
                else "irreversibleArrheniusLindemannFallOff"
            )

        raise ValueError(f"Unknown reaction type: {rtype}")

    def _arrhenius_coeffs(self, index: int):
        return self.get_arrhenius_parameters(index)

    def _arrhenius_from_rate(self, rate):
        """Return Arrhenius parameters from a Cantera ``Rate`` object."""
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

    def _plog_coeffs(self, index: int) -> list:
        """Return pressure-dependent Arrhenius coefficients for a reaction."""

        rxn = self.gas.reactions()[index]
        if not isinstance(rxn.rate, ct.reaction.PlogRate):
            raise ValueError("Reaction does not use PLOG form")

        coeffs = []
        for P, arr in rxn.rate.rates:
            A, Ta, b = self._arrhenius_from_rate(arr)
            coeffs.append((P, A, b, Ta))
        return coeffs

    def is_plog(self, index: int) -> bool:
        """Return ``True`` if the reaction uses pressure-dependent Arrhenius rates."""

        rxn = self.gas.reactions()[index]
        return isinstance(getattr(rxn, "rate", None), ct.reaction.PlogRate)

    def _plog_block(self, coeffs) -> list:
        """Return the ``ArrheniusData`` block for pressure-dependent rates."""

        lines = ["    ArrheniusData", "    (", "        // convert MPa to Pa"]
        for P, A, b, Ta in coeffs:
            P_pa = P * 1e6
            lines.append(
                f"        ({P_pa:.6g}  {A:.6g} {b:.6g} {Ta:.6g} )   // PLOG /p A beta Ta/"
            )
        lines.append("    );")
        return lines

    def _plog_reaction_string(self, index: int) -> str:
        """Return an OpenFOAM block for a pressure-dependent Arrhenius reaction."""

        coeffs = self._plog_coeffs(index)
        A0, Ta0, b0 = coeffs[0][1], coeffs[0][3], coeffs[0][2]
        lines = self._reaction_block_start(
            index, "pressureDependentArrhenius", A0, b0, Ta0
        )
        lines += self._plog_block(coeffs)
        lines.append(self._block_end())
        return "\n".join(lines)

    def _reaction_equation_string(self, rxn) -> str:
        left = []
        orders = getattr(rxn, "orders", {})
        for sp, nu in rxn.reactants.items():
            coeff = f"{round(nu)}" if (np.abs(nu-1)>1e-6) else ""
            order = orders.get(sp)
            if order is not None:
                left.append(f"{coeff}{sp}^{order}")
            else:
                left.append(f"{coeff}{sp}")
        right = []
        for sp, nu in rxn.products.items():
            coeff = f"{round(nu)}" if (np.abs(nu-1)>1e-6) else ""
            right.append(f"{coeff}{sp}")
        return " + ".join(left) + " = " + " + ".join(right)

    def _reaction_block_start(self, index: int, typ: str, A: float, b: float, Ta: float):
        rxn = self.gas.reactions()[index]
        eq = self._reaction_equation_string(rxn)
        lines = [
            f"un-named-reaction-{index}",
            "{",
            f"    type            {typ};",
            f"    reaction        \"{eq}\";",
            f"    A               {A:.6e};",
            f"    beta            {b:.6f};",
            f"    Ta              {Ta:.6f};",
        ]
        return lines

    def _block_end(self) -> str:
        return "}"

    def _irreversible_reaction_string(self, index: int) -> str:
        A, Ta, b = self._arrhenius_coeffs(index)
        lines = self._reaction_block_start(index, "irreversibleArrhenius", A, b, Ta)
        lines.append(self._block_end())
        return "\n".join(lines)

    def _reversible_reaction_string(self, index: int) -> str:
        """Return an OpenFOAM reaction block for a reversible Arrhenius reaction."""

        A, Ta, b = self._arrhenius_coeffs(index)
        lines = self._reaction_block_start(index, "reversibleArrhenius", A, b, Ta)
        lines.append(self._block_end())
        return "\n".join(lines)

    # --- Extended reaction writers -----------------------------------------


    def _third_body_coeff_lines(self, effs: dict) -> list:
        # start with the dummy third body efficiencies
        effs_true = deepcopy(self.third_body_dict)
        
        for sp, values in effs.items():
            effs_true[sp] = values
            
        effs = effs_true
        indent = f"    "
        lines = [indent + str(len(effs))]
        lines.append(indent  + "(")
        for sp, val in effs.items():
            lines.append(indent*2  +f"({sp} {val})")
        lines += ["    )", "    ;"]
     
        return lines

    def _third_body_efficiencies_block(self, effs: dict) -> list:
        #lines = ["    thirdBodyEfficiencies", "    {"]
        lines = []
        inner = self._third_body_coeff_lines(effs)
        lines.extend(["        " + l.strip() if i > 0 else l for i, l in enumerate(inner)])
        lines.append("    }")
        return lines

    def _arrhenius_sub_block(self, name: str, rate) -> list:
        A, Ta, b = self._arrhenius_from_rate(rate)
        lines = [f"    {name}", "    {", f"        A               {A};", f"        beta            {b};", f"        Ta              {Ta};", "    }"]
        return lines

    def _troe_params_block(self, coeffs) -> list:
        """Return the ``F`` block for Troe falloff coefficients."""

        # Some mechanisms omit the ``Tss`` value
        alpha, tsss, ts, *rest = coeffs
        tss = rest[0] if rest else None

        lines = [
            "    F",
            "    {",
            f"        alpha           {alpha};",
            f"        Tsss            {tsss};",
            f"        Ts              {ts};",
        ]
        if tss is not None:
            lines.append(f"        Tss             {tss};")
        lines.append("    }")
        return lines

    def _third_body_reaction_string(self, index: int) -> str:
        rxn = self.gas.reactions()[index]
        A, Ta, b = self._arrhenius_coeffs(index)
        typ = (
            "reversibleThirdBodyArrhenius" if self.is_reversible(index) else "irreversibleThirdBodyArrhenius"
        )
        lines = self._reaction_block_start(index, typ, A, b, Ta)
        lines += self._third_body_coeff_lines(rxn.third_body.efficiencies)
        lines.append(self._block_end())
        return "\n".join(lines)

    def _lindemann_reaction_string(self, index: int) -> str:
        rxn = self.gas.reactions()[index]
        typ = (
            "reversibleArrheniusLindemannFallOff" if self.is_reversible(index) else "irreversibleArrheniusLindemannFallOff"
        )
        eq = self._reaction_equation_string(rxn)
        lines = [f"un-named-reaction-{index}", "{", f"    type            {typ};", f"    reaction        \"{eq}\";"]
        rate = rxn.rate
        lines += self._arrhenius_sub_block("k0", rate.low_rate)
        lines += self._arrhenius_sub_block("kInf", rate.high_rate)
        lines += ["    F", "    {", "    }"]
        lines += self._third_body_efficiencies_block(rxn.third_body.efficiencies)
        lines.append("}")
        return "\n".join(lines)

    def _troe_reaction_string(self, index: int) -> str:
        rxn = self.gas.reactions()[index]
        typ = (
            "reversibleArrheniusTroeFallOff" if self.is_reversible(index) else "irreversibleArrheniusTroeFallOff"
        )
        eq = self._reaction_equation_string(rxn)
        lines = [f"un-named-reaction-{index}", "{", f"    type            {typ};", f"    reaction        \"{eq}\";"]
        rate = rxn.rate
        lines += self._arrhenius_sub_block("k0", rate.low_rate)
        lines += self._arrhenius_sub_block("kInf", rate.high_rate)
        lines += self._troe_params_block(rate.falloff_coeffs)
        lines += self._third_body_efficiencies_block(rxn.third_body.efficiencies)
        lines.append("}")
        return "\n".join(lines)

    def chemistry_file_string(self, indices=None) -> str:
        lines = ["reactions", "{"]
        if indices is None:
            indices = range(len(self.gas.reactions()))
        for i in indices:
            typ = self.check_reaction_type(i)
            if typ == "irreversibleArrhenius":
                block = self._irreversible_reaction_string(i)
            elif typ == "reversibleArrhenius":
                block = self._reversible_reaction_string(i)
            elif typ in ("reversibleThirdBodyArrhenius", "irreversibleThirdBodyArrhenius"):
                block = self._third_body_reaction_string(i)
            elif typ in (
                "reversibleArrheniusLindemannFallOff",
                "irreversibleArrheniusLindemannFallOff",
            ):
                block = self._lindemann_reaction_string(i)
            elif typ in (
                "reversibleArrheniusTroeFallOff",
                "irreversibleArrheniusTroeFallOff",
            ):
                block = self._troe_reaction_string(i)
            elif typ == "pressureDependentArrhenius":
                block = self._plog_reaction_string(i)
            else:
                block = ""
            for ln in block.splitlines():
                lines.append("    " + ln)
        lines.append("}")
        lines.append("Tlow            200;")
        lines.append("Thigh           6000;")
        return "\n".join(lines)

    