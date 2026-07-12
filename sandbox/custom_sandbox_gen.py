"""Generate custom PFLOTRAN reaction sandboxes from a reaction definition.

Writes a timestamped folder under sandbox/ containing:
  - hanford.dat (base DB + any new basis species)
  - reaction_sandbox_<name>.F90
  - reaction_sandbox_block.in
  - reaction_def.json
  - README.md

Kinetic form (v1): power-law over reactants with optional water-activity
inhibition, matching the AWINHIBIT family.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

SOLVENT_NAMES = frozenset({"H2O", "H2O(l)"})
SANDBOX_DIR = Path(__file__).resolve().parent
DEFAULT_HANFORD = SANDBOX_DIR / "hanford.dat"

_NAME_RE = re.compile(r"^[a-z][a-z0-9_]*$")
_QUOTED_NAME_RE = re.compile(r"^'([^']+)'")


@dataclass
class NewSpecies:
    """Basis-species line fields for hanford.dat."""

    name: str
    a0: float = 3.0
    charge: float = 0.0
    molecular_weight: float = 0.0


@dataclass
class ReactionDef:
    """User-facing definition of a custom reaction sandbox."""

    name: str
    reactants: dict[str, float]
    products: dict[str, float]
    description: str = ""
    new_species: list[NewSpecies] = field(default_factory=list)
    rate_constant: float = 1.0e-10
    aw_threshold: float = 0.5
    inhibition_type: str = "THRESHOLD"  # THRESHOLD | SMOOTHSTEP
    activation_energy: float | None = None
    reference_temperature: float | None = None

    def validate(self) -> None:
        if not _NAME_RE.match(self.name):
            raise ValueError(
                "name must be lowercase alphanumeric/underscore, "
                f"starting with a letter (got {self.name!r})"
            )
        if not self.reactants:
            raise ValueError("reactants must be non-empty")
        if not self.products:
            raise ValueError("products must be non-empty")
        for side, mapping in (
            ("reactants", self.reactants),
            ("products", self.products),
        ):
            for sp, coeff in mapping.items():
                if coeff <= 0:
                    raise ValueError(f"{side}[{sp!r}] coefficient must be > 0")
        inhib = self.inhibition_type.upper()
        if inhib not in {"THRESHOLD", "SMOOTHSTEP"}:
            raise ValueError("inhibition_type must be THRESHOLD or SMOOTHSTEP")
        self.inhibition_type = inhib
        if self.activation_energy is not None and self.reference_temperature is None:
            raise ValueError(
                "reference_temperature is required when activation_energy is set"
            )


@dataclass
class HanfordParsed:
    """Parsed sections of a hanford.dat thermodynamic database."""

    text: str
    basis: set[str]
    aqueous: set[str]
    gases: set[str]
    minerals: set[str]
    surface: set[str]

    @property
    def all_species(self) -> set[str]:
        return (
            self.basis
            | self.aqueous
            | self.gases
            | self.minerals
            | self.surface
            | SOLVENT_NAMES
        )


def _pascal(name: str) -> str:
    return "".join(part.capitalize() for part in name.split("_"))


def _keyword(name: str) -> str:
    return name.upper().replace("_", "")


def species_slug(species: str) -> str:
    """Map a chemistry name to a safe Fortran identifier fragment."""
    s = species.lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = s.strip("_")
    if not s:
        s = "sp"
    if s[0].isdigit():
        s = "sp_" + s
    return s


def _fmt_d(value: float) -> str:
    """Format a float as a Fortran double literal."""
    if float(value).is_integer():
        return f"{int(value)}.d0"
    text = f"{value:.10g}".replace("e", "d").replace("E", "d")
    if "d" not in text:
        text = text + "d0"
    return text


def reaction_equation(reactants: dict[str, float], products: dict[str, float]) -> str:
    def _side(mapping: dict[str, float]) -> str:
        parts = []
        for name, coeff in mapping.items():
            if coeff == 1:
                parts.append(name)
            else:
                c = int(coeff) if float(coeff).is_integer() else coeff
                parts.append(f"{c} {name}")
        return " + ".join(parts)

    return f"{_side(reactants)} -> {_side(products)}"


def parse_hanford_dat(path: str | Path) -> HanfordParsed:
    """Parse hanford.dat into named-species sets per section."""
    path = Path(path)
    text = path.read_text()
    lines = text.splitlines()

    sections: list[list[str]] = []
    current: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("'null'"):
            sections.append(current)
            current = []
        else:
            current.append(line)
    if current:
        sections.append(current)

    def names_from(section_lines: list[str], skip_temp: bool = False) -> set[str]:
        out: set[str] = set()
        for i, line in enumerate(section_lines):
            if skip_temp and i == 0 and "temperature points" in line:
                continue
            m = _QUOTED_NAME_RE.match(line.strip())
            if m:
                out.add(m.group(1))
        return out

    basis = names_from(sections[0], skip_temp=True) if sections else set()
    aqueous = names_from(sections[1]) if len(sections) > 1 else set()
    gases = names_from(sections[2]) if len(sections) > 2 else set()
    minerals = names_from(sections[3]) if len(sections) > 3 else set()
    surface = names_from(sections[4]) if len(sections) > 4 else set()

    return HanfordParsed(
        text=text,
        basis=basis,
        aqueous=aqueous,
        gases=gases,
        minerals=minerals,
        surface=surface,
    )


def species_status(
    names: list[str] | set[str], parsed: HanfordParsed
) -> list[dict[str, str]]:
    """Classify each species relative to the database."""
    rows = []
    for name in sorted(names, key=str.lower):
        if name in SOLVENT_NAMES:
            status = "solvent"
        elif name in parsed.basis:
            status = "in_basis"
        elif name in parsed.aqueous:
            status = "in_aqueous_complexes"
        elif name in parsed.gases:
            status = "in_gases"
        elif name in parsed.minerals:
            status = "in_minerals"
        elif name in parsed.surface:
            status = "in_surface"
        else:
            status = "missing"
        rows.append({"species": name, "status": status})
    return rows


def tracked_species(defn: ReactionDef) -> list[str]:
    """Species that need Fortran indices (exclude solvent)."""
    names = list(defn.reactants) + list(defn.products)
    seen: list[str] = []
    for n in names:
        if n in SOLVENT_NAMES:
            continue
        if n not in seen:
            seen.append(n)
    return seen


def add_basis_species(dat_text: str, species: list[NewSpecies] | NewSpecies) -> str:
    """Insert new basis species lines before the first 'null' sentinel."""
    if isinstance(species, NewSpecies):
        species = [species]
    if not species:
        return dat_text

    lines = dat_text.splitlines(keepends=True)
    insert_at = None
    for i, line in enumerate(lines):
        if line.strip().startswith("'null'"):
            insert_at = i
            break
    if insert_at is None:
        raise ValueError("Could not find first 'null' sentinel in hanford.dat")

    new_lines = []
    for sp in species:
        new_lines.append(
            f"'{sp.name}' {_fmt_plain(sp.a0)} {_fmt_plain(sp.charge)} "
            f"{_fmt_plain(sp.molecular_weight)}\n"
        )
    return "".join(lines[:insert_at] + new_lines + lines[insert_at:])


def _fmt_plain(value: float) -> str:
    if float(value).is_integer():
        return f"{float(value):.1f}"
    return f"{value:g}"


def _symbols(defn: ReactionDef) -> dict[str, str]:
    pascal = _pascal(defn.name)
    keyword = _keyword(defn.name)
    return {
        "name": defn.name,
        "pascal": pascal,
        "keyword": keyword,
        "module": f"Reaction_Sandbox_{pascal}_class",
        "type": f"reaction_sandbox_{defn.name}_type",
        "create": f"{pascal}Create",
        "read": f"{pascal}Read",
        "setup": f"{pascal}Setup",
        "evaluate": f"{pascal}Evaluate",
        "destroy": f"{pascal}Destroy",
        "param_threshold": f"{keyword}_THRESHOLD_INHIBITION",
        "param_smoothstep": f"{keyword}_SMOOTHSTEP_INHIBITION",
        "file": f"reaction_sandbox_{defn.name}.F90",
    }


def render_f90(defn: ReactionDef) -> str:
    """Emit a full AWINHIBIT-shaped reaction sandbox module."""
    defn.validate()
    sym = _symbols(defn)
    tracked = tracked_species(defn)
    if not tracked:
        raise ValueError("reaction has no non-solvent aqueous species to track")

    # Ensure unique Fortran field names
    slugs: dict[str, str] = {}
    used: set[str] = set()
    for sp in tracked:
        base = species_slug(sp)
        slug = base
        n = 2
        while slug in used:
            slug = f"{base}_{n}"
            n += 1
        used.add(slug)
        slugs[sp] = slug

    eq = reaction_equation(defn.reactants, defn.products)
    desc = defn.description or eq

    index_decls = "\n".join(f"    PetscInt :: i_{slugs[sp]}" for sp in tracked)
    index_inits = "\n".join(
        f"  {sym['create']}%i_{slugs[sp]} = UNINITIALIZED_INTEGER" for sp in tracked
    )

    setup_lookups = []
    for sp in tracked:
        setup_lookups.append(
            f"  word = '{sp}'\n"
            f"  this%i_{slugs[sp]} = &\n"
            f"    ReactionAuxGetPriSpecIDFromName(word,reaction,option)"
        )
    setup_block = "\n\n".join(setup_lookups)

    conc_decls = ", ".join(f"C_{slugs[sp]}" for sp in tracked)
    conc_assigns = []
    for sp in tracked:
        s = slugs[sp]
        conc_assigns.append(
            f"  C_{s} = rt_auxvar%pri_molal(this%i_{s}) * &\n"
            f"         rt_auxvar%pri_act_coef(this%i_{s}) * molality_to_molarity"
        )
    conc_block = "\n".join(conc_assigns)

    rate_terms = []
    for sp, coeff in defn.reactants.items():
        if sp in SOLVENT_NAMES:
            continue
        s = slugs[sp]
        if coeff == 1:
            rate_terms.append(f"C_{s}")
        else:
            rate_terms.append(f"(C_{s}**{_fmt_d(coeff)})")
    if not rate_terms:
        raise ValueError("need at least one non-solvent reactant for the rate law")
    rate_product = " * ".join(rate_terms)

    residual_lines = []
    for sp, coeff in defn.reactants.items():
        if sp in SOLVENT_NAMES:
            residual_lines.append(f"  ! {coeff} {sp} consumed (solvent, not tracked)")
            continue
        s = slugs[sp]
        residual_lines.append(f"  ! {coeff} {sp} consumed")
        if coeff == 1:
            residual_lines.append(
                f"  Residual(this%i_{s}) = Residual(this%i_{s}) + reaction_rate"
            )
        else:
            residual_lines.append(
                f"  Residual(this%i_{s}) = Residual(this%i_{s}) + "
                f"{_fmt_d(coeff)} * reaction_rate"
            )
    for sp, coeff in defn.products.items():
        if sp in SOLVENT_NAMES:
            residual_lines.append(f"  ! {coeff} {sp} produced (solvent, not tracked)")
            continue
        s = slugs[sp]
        residual_lines.append(f"  ! {coeff} {sp} produced")
        if coeff == 1:
            residual_lines.append(
                f"  Residual(this%i_{s}) = Residual(this%i_{s}) - reaction_rate"
            )
        else:
            residual_lines.append(
                f"  Residual(this%i_{s}) = Residual(this%i_{s}) - "
                f"{_fmt_d(coeff)} * reaction_rate"
            )
    residual_block = "\n".join(residual_lines)

    # Keep error_string short (< 32 chars for some PFLOTRAN paths)
    err = f"CHEMISTRY,RXN_SANDBOX,{sym['keyword']}"
    if len(err) > 48:
        err = f"RXN_SANDBOX,{sym['keyword']}"

    return f"""\
module {sym['module']}

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  PetscInt, parameter :: {sym['param_threshold']} = 1
  PetscInt, parameter :: {sym['param_smoothstep']} = 2

  type, public, &
    extends(reaction_sandbox_base_type) :: {sym['type']}

    ! Water activity inhibition parameters
    PetscReal :: aw_threshold
    PetscInt :: inhibition_type

    ! Reaction parameters: {eq}
    PetscReal :: rate_constant
    PetscReal :: activation_energy
    PetscReal :: reference_temperature

    ! Species indices
{index_decls}

  contains
    procedure, public :: ReadInput  => {sym['read']}
    procedure, public :: Setup      => {sym['setup']}
    procedure, public :: Evaluate   => {sym['evaluate']}
    procedure, public :: Destroy    => {sym['destroy']}
  end type {sym['type']}

  public :: {sym['create']}

contains

! ************************************************************************** !

function {sym['create']}()
  !
  ! Allocates {sym['pascal']} reaction object
  ! {desc}
  !

  implicit none

  class({sym['type']}), pointer :: {sym['create']}

  allocate({sym['create']})

  {sym['create']}%aw_threshold = 0.5d0
  {sym['create']}%inhibition_type = {sym['param_threshold']}

  {sym['create']}%rate_constant = UNINITIALIZED_DOUBLE
  {sym['create']}%activation_energy = UNINITIALIZED_DOUBLE
  {sym['create']}%reference_temperature = UNINITIALIZED_DOUBLE

{index_inits}

  nullify({sym['create']}%next)

end function {sym['create']}

! ************************************************************************** !

subroutine {sym['read']}(this,input,option)
  !
  ! Reads input deck for {sym['pascal']} parameters
  !

  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class({sym['type']}) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string
  error_string = '{err}'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       trim(error_string))
    call StringToUpper(word)

    select case(trim(word))
      case('WATER_ACTIVITY_THRESHOLD')
        call InputReadDouble(input,option,this%aw_threshold)
        call InputErrorMsg(input,option,'water_activity_threshold',error_string)
        if (this%aw_threshold < 0.d0 .or. this%aw_threshold > 1.d0) then
          option%io_buffer = 'WATER_ACTIVITY_THRESHOLD must be between 0 and 1'
          call PrintErrMsg(option)
        endif

      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant',error_string)
        call InputReadAndConvertUnits(input,this%rate_constant,'mol/m^3-sec',&
                        trim(error_string)//',rate_constant',option)

      case('INHIBITION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('THRESHOLD')
            this%inhibition_type = {sym['param_threshold']}
          case('SMOOTHSTEP')
            this%inhibition_type = {sym['param_smoothstep']}
          case default
            error_string = trim(error_string) // ',INHIBITION_TYPE'
            call InputKeywordUnrecognized(input,word,error_string ,option)
        end select

      case('ACTIVATION_ENERGY')
        call InputReadDouble(input,option,this%activation_energy)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%activation_energy,'j/mol',&
                          trim(error_string)//',activation energy',option)

      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,this%reference_temperature)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%reference_temperature,'C',&
                          trim(error_string)//',reference temperature',option)

      case default
        call InputKeywordUnrecognized(input,word,error_string ,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine {sym['read']}

! ************************************************************************** !

subroutine {sym['setup']}(this,reaction,option)
  !
  ! Sets up the {sym['pascal']} reaction
  !

  use Option_module
  use Utility_module
  use Reaction_Aux_module

  implicit none

  class({sym['type']}) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT must be provided'
    call PrintErrMsg(option)
  endif

  ! Species indices: {eq}
{setup_block}

  if (Initialized(this%activation_energy) .and. &
      UnInitialized(this%reference_temperature)) then
    option%io_buffer = 'REFERENCE_TEMPERATURE required with ACTIVATION_ENERGY'
    call PrintErrMsg(option)
  endif

end subroutine {sym['setup']}

! ************************************************************************** !

subroutine {sym['evaluate']}(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
  !
  ! Evaluates: {eq}
  ! Inhibited below water activity threshold
  !

  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  class({sym['type']}) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  PetscReal :: molality_to_molarity
  PetscReal :: water_activity, aw_inhibition, tempreal
  PetscReal :: rate_constant
  PetscReal :: reaction_rate

  PetscReal :: {conc_decls}

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! m^3 -> L

  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3

  water_activity = exp(rt_auxvar%ln_act_h2o)

  rate_constant = this%rate_constant
  if (Initialized(this%activation_energy)) then
    rate_constant = rate_constant * Arrhenius(this%activation_energy, &
                                            global_auxvar%temp, &
                                            this%reference_temperature)
  endif

{conc_block}

  select case(this%inhibition_type)
    case({sym['param_smoothstep']})
      call ReactionInhibitionSmoothstep(water_activity, this%aw_threshold, &
                                        0.05d0, aw_inhibition, tempreal)
      aw_inhibition = 1.d0 - aw_inhibition
    case({sym['param_threshold']})
      if (water_activity < this%aw_threshold) then
        aw_inhibition = 0.d0
      else
        aw_inhibition = 1.d0
      endif
  end select

  reaction_rate = rate_constant * {rate_product} * aw_inhibition
  reaction_rate = reaction_rate * L_water  ! Convert to mol/sec

{residual_block}

  if (compute_derivative) then
    option%io_buffer = 'REACTION_SANDBOX {sym['keyword']} needs NUMERICAL_JACOBIAN'
    call PrintErrMsg(option)
  endif

end subroutine {sym['evaluate']}

! ************************************************************************** !

subroutine {sym['destroy']}(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !

  implicit none
  class({sym['type']}) :: this

end subroutine {sym['destroy']}

end module {sym['module']}
"""


def render_in_snippet(defn: ReactionDef) -> str:
    """Emit a REACTION_SANDBOX block for pasting into a .in deck."""
    defn.validate()
    sym = _symbols(defn)
    lines = [
        "REACTION_SANDBOX",
        f"  {sym['keyword']}",
        f"    WATER_ACTIVITY_THRESHOLD {defn.aw_threshold}",
        f"    RATE_CONSTANT {defn.rate_constant}  ! mol/m^3-sec",
        f"    INHIBITION_TYPE {defn.inhibition_type}",
    ]
    if defn.activation_energy is not None:
        lines.append(f"    ACTIVATION_ENERGY {defn.activation_energy}  ! J/mol")
        lines.append(
            f"    REFERENCE_TEMPERATURE {defn.reference_temperature}  ! C"
        )
    lines.extend(
        [
            "  /",
            "/",
            "",
            "! Also required in NUMERICAL_METHODS TRANSPORT / NEWTON_SOLVER:",
            "!   NUMERICAL_JACOBIAN",
            "!",
            f"! Reaction: {reaction_equation(defn.reactants, defn.products)}",
            "! Point DATABASE at this folder's hanford.dat",
        ]
    )
    return "\n".join(lines) + "\n"


def render_readme(defn: ReactionDef, out_dir: Path) -> str:
    sym = _symbols(defn)
    eq = reaction_equation(defn.reactants, defn.products)
    rel = out_dir.name
    return f"""# Custom sandbox: {sym['keyword']}

Generated reaction: `{eq}`

## Files

| File | Role |
|------|------|
| `{sym['file']}` | Fortran reaction module |
| `hanford.dat` | Thermodynamic DB (base + any new species) |
| `reaction_sandbox_block.in` | Snippet to paste into your `.in` deck |
| `reaction_def.json` | Reproducible definition used to generate this folder |

## Install into PFLOTRAN

From the repo root:

```bash
python3 scripts/patch_pflotran_sandboxes.py \\
  --pflotran-src /path/to/pflotran/src/pflotran \\
  --sandbox-dir sandbox/ \\
  --extra-sandbox-dir sandbox/{rel}

cd /path/to/pflotran/src/pflotran && make clean && make pflotran
```

## Use in an input deck

1. Set `DATABASE` to this folder's `hanford.dat`.
2. Paste the contents of `reaction_sandbox_block.in` into the chemistry section.
3. List every tracked species under `PRIMARY_SPECIES`.
4. Enable `NUMERICAL_JACOBIAN` under `NUMERICAL_METHODS TRANSPORT` / `NEWTON_SOLVER`.
"""


def _defn_to_jsonable(defn: ReactionDef) -> dict[str, Any]:
    data = asdict(defn)
    return data


def reaction_def_from_dict(data: dict[str, Any]) -> ReactionDef:
    new_species = [NewSpecies(**ns) for ns in data.get("new_species", [])]
    return ReactionDef(
        name=data["name"],
        reactants=dict(data["reactants"]),
        products=dict(data["products"]),
        description=data.get("description", ""),
        new_species=new_species,
        rate_constant=float(data.get("rate_constant", 1.0e-10)),
        aw_threshold=float(data.get("aw_threshold", 0.5)),
        inhibition_type=str(data.get("inhibition_type", "THRESHOLD")),
        activation_energy=data.get("activation_energy"),
        reference_temperature=data.get("reference_temperature"),
    )


def write_custom_sandbox(
    defn: ReactionDef,
    out_root: str | Path | None = None,
    hanford_path: str | Path | None = None,
    timestamp: str | None = None,
) -> Path:
    """Write a custom_<timestamp>/ folder with all generated artifacts."""
    defn.validate()
    out_root = Path(out_root) if out_root else SANDBOX_DIR
    hanford_path = Path(hanford_path) if hanford_path else DEFAULT_HANFORD

    parsed = parse_hanford_dat(hanford_path)
    all_rxn = set(defn.reactants) | set(defn.products)
    audit = species_status(all_rxn, parsed)
    missing = {row["species"] for row in audit if row["status"] == "missing"}
    provided = {sp.name for sp in defn.new_species}
    still_missing = missing - provided
    if still_missing:
        raise ValueError(
            "Species missing from hanford.dat and not in new_species: "
            + ", ".join(sorted(still_missing))
        )

    # Warn if species only appear as complexes (need to be primary in deck)
    for row in audit:
        if row["status"] == "in_aqueous_complexes":
            # Allowed if user will list as primary; still need DB entry — complexes
            # are fine for lookup existence, but Setup requires primary species.
            pass

    ts = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = out_root / f"custom_{ts}"
    out_dir.mkdir(parents=True, exist_ok=False)

    dat_text = add_basis_species(parsed.text, defn.new_species)
    sym = _symbols(defn)

    (out_dir / "hanford.dat").write_text(dat_text)
    (out_dir / sym["file"]).write_text(render_f90(defn))
    (out_dir / "reaction_sandbox_block.in").write_text(render_in_snippet(defn))
    (out_dir / "reaction_def.json").write_text(
        json.dumps(_defn_to_jsonable(defn), indent=2) + "\n"
    )
    (out_dir / "README.md").write_text(render_readme(defn, out_dir))

    return out_dir


def search_species(parsed: HanfordParsed, query: str, limit: int = 30) -> list[str]:
    """Case-insensitive substring search across all species names."""
    q = query.lower()
    hits = sorted(
        (n for n in parsed.all_species if q in n.lower()),
        key=str.lower,
    )
    return hits[:limit]
