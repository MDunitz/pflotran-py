#!/usr/bin/env python3
"""Patch a PFLOTRAN source tree to include pflotran-py custom reaction sandboxes.

Copies the three AWINHIBIT*.F90 modules from sandbox/ into PFLOTRAN's src
directory and updates reaction_sandbox.F90, pflotran_object_files.txt, and
pflotran_dependencies.txt. Idempotent: safe to run multiple times.

Optionally install one or more generated custom sandboxes via
--extra-sandbox-dir (e.g. sandbox/custom_YYYYMMDD_HHMMSS/).

Usage:
    python3 scripts/patch_pflotran_sandboxes.py \\
        --pflotran-src /scratch/pflotran/src/pflotran \\
        --sandbox-dir /work/sandbox

    python3 scripts/patch_pflotran_sandboxes.py \\
        --pflotran-src /scratch/pflotran/src/pflotran \\
        --sandbox-dir /work/sandbox \\
        --extra-sandbox-dir /work/sandbox/custom_20260712_143022
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

SANDBOX_FILES = [
    "reaction_sandbox_awinhibit.F90",
    "reaction_sandbox_awinhibitacetate.F90",
    "reaction_sandbox_awinhibitmethyl.F90",
]

USE_LINES = [
    "  use Reaction_Sandbox_AWInhibit_class",
    "  use Reaction_Sandbox_AWInhibitAcetate_class",
    "  use Reaction_Sandbox_AWInhibitMethyl_class",
]

CASE_LINES = [
    "      case('AWINHIBIT')",
    "        new_sandbox => AWInhibitCreate()",
    "      case('AWINHIBITACETATE')",
    "        new_sandbox => AWInhibitAcetateCreate()",
    "      case('AWINHIBITMETHYL')",
    "        new_sandbox => AWInhibitMethylCreate()",
]

OBJECT_LINES = [
    "\t${common_src}reaction_sandbox_awinhibit.o \\",
    "\t${common_src}reaction_sandbox_awinhibitacetate.o \\",
    "\t${common_src}reaction_sandbox_awinhibitmethyl.o \\",
]

DEP_BLOCK = """\
reaction_sandbox_awinhibit.o : \\
  reaction_sandbox_base.o \\
  reactive_transport_aux.o \\
  global_aux.o \\
  reaction_aux.o
reaction_sandbox_awinhibitacetate.o : \\
  reaction_sandbox_base.o \\
  reactive_transport_aux.o \\
  global_aux.o \\
  reaction_aux.o
reaction_sandbox_awinhibitmethyl.o : \\
  reaction_sandbox_base.o \\
  reactive_transport_aux.o \\
  global_aux.o \\
  reaction_aux.o
"""

_F90_STEM_RE = re.compile(r"^reaction_sandbox_(.+)\.F90$")
_MODULE_RE = re.compile(
    r"^\s*module\s+(Reaction_Sandbox_(\w+)_class)\s*$", re.IGNORECASE
)
_CREATE_RE = re.compile(r"^\s*public\s*::\s*(\w+Create)\s*$", re.IGNORECASE)


def copy_sandbox_sources(sandbox_dir: Path, pflotran_src: Path) -> None:
    for name in SANDBOX_FILES:
        shutil.copy2(sandbox_dir / name, pflotran_src / name)


def patch_reaction_sandbox_f90(path: Path) -> None:
    text = path.read_text()
    if "Reaction_Sandbox_AWInhibit_class" in text:
        return

    anchor = "  use Reaction_Sandbox_Radon_class"
    if anchor not in text:
        raise RuntimeError(f"Could not find use-anchor in {path}")
    text = text.replace(
        anchor,
        anchor + "\n" + "\n".join(USE_LINES),
        1,
    )

    anchor = "      case default"
    if anchor not in text:
        raise RuntimeError(f"Could not find case-anchor in {path}")
    text = text.replace(
        anchor,
        "\n".join(CASE_LINES) + "\n" + anchor,
        1,
    )
    path.write_text(text)


def patch_object_files(path: Path) -> None:
    text = path.read_text()
    if "reaction_sandbox_awinhibit.o" in text:
        return

    anchor = "\t${common_src}reaction_sandbox_ufd_wp.o \\"
    if anchor not in text:
        raise RuntimeError(f"Could not find object-file anchor in {path}")
    text = text.replace(
        anchor,
        anchor + "\n" + "\n".join(OBJECT_LINES),
        1,
    )
    path.write_text(text)


def patch_dependencies(path: Path) -> None:
    text = path.read_text()
    if "reaction_sandbox_awinhibit.o" in text:
        return

    anchor = "  reaction_sandbox_ufd_wp.o \\"
    if anchor not in text:
        raise RuntimeError(f"Could not find dependency anchor in {path}")
    text = text.replace(
        anchor,
        anchor
        + "\n  reaction_sandbox_awinhibit.o\\\n"
        + "  reaction_sandbox_awinhibitacetate.o\\\n"
        + "  reaction_sandbox_awinhibitmethyl.o\\",
        1,
    )
    text = text.rstrip() + "\n" + DEP_BLOCK + "\n"
    path.write_text(text)


def _discover_extra_modules(extra_dir: Path) -> list[dict[str, str]]:
    """Find generated reaction_sandbox_*.F90 modules in a custom_* folder."""
    modules: list[dict[str, str]] = []
    for path in sorted(extra_dir.glob("reaction_sandbox_*.F90")):
        stem_m = _F90_STEM_RE.match(path.name)
        if not stem_m:
            continue
        stem = stem_m.group(1)
        if stem in {"awinhibit", "awinhibitacetate", "awinhibitmethyl", "template", "aq"}:
            continue

        text = path.read_text()
        module_name = None
        create_name = None
        for line in text.splitlines():
            if module_name is None:
                m = _MODULE_RE.match(line)
                if m:
                    module_name = m.group(1)
            if create_name is None:
                m = _CREATE_RE.match(line)
                if m:
                    create_name = m.group(1)
            if module_name and create_name:
                break
        if not module_name or not create_name:
            raise RuntimeError(
                f"Could not parse module/create symbols from {path}"
            )

        # Keyword: uppercase stem with underscores removed (matches generator)
        keyword = stem.upper().replace("_", "")
        modules.append(
            {
                "file": path.name,
                "stem": stem,
                "module": module_name,
                "create": create_name,
                "keyword": keyword,
                "path": str(path),
            }
        )
    if not modules:
        raise RuntimeError(f"No reaction_sandbox_*.F90 found in {extra_dir}")
    return modules


def patch_extra_sandbox(pflotran_src: Path, extra_dir: Path) -> list[str]:
    """Copy and register generated modules from a custom_* directory."""
    modules = _discover_extra_modules(extra_dir)
    installed: list[str] = []

    rs_path = pflotran_src / "reaction_sandbox.F90"
    obj_path = pflotran_src / "pflotran_object_files.txt"
    dep_path = pflotran_src / "pflotran_dependencies.txt"

    rs_text = rs_path.read_text()
    obj_text = obj_path.read_text()
    dep_text = dep_path.read_text()

    for mod in modules:
        src = Path(mod["path"])
        shutil.copy2(src, pflotran_src / mod["file"])

        use_line = f"  use {mod['module']}"
        if use_line not in rs_text:
            anchor = "  use Reaction_Sandbox_Radon_class"
            if anchor not in rs_text:
                # Fall back: after AWINHIBIT methyl use if present
                anchor = "  use Reaction_Sandbox_AWInhibitMethyl_class"
            if anchor not in rs_text:
                raise RuntimeError(f"Could not find use-anchor in {rs_path}")
            rs_text = rs_text.replace(anchor, anchor + "\n" + use_line, 1)

        case_block = (
            f"      case('{mod['keyword']}')\n"
            f"        new_sandbox => {mod['create']}()"
        )
        if f"case('{mod['keyword']}')" not in rs_text:
            anchor = "      case default"
            if anchor not in rs_text:
                raise RuntimeError(f"Could not find case-anchor in {rs_path}")
            rs_text = rs_text.replace(anchor, case_block + "\n" + anchor, 1)

        obj_line = f"\t${{common_src}}reaction_sandbox_{mod['stem']}.o \\"
        if f"reaction_sandbox_{mod['stem']}.o" not in obj_text:
            # Prefer after awinhibitmethyl if present, else ufd_wp
            for anchor in (
                "\t${common_src}reaction_sandbox_awinhibitmethyl.o \\",
                "\t${common_src}reaction_sandbox_ufd_wp.o \\",
            ):
                if anchor in obj_text:
                    obj_text = obj_text.replace(anchor, anchor + "\n" + obj_line, 1)
                    break
            else:
                raise RuntimeError(f"Could not find object-file anchor in {obj_path}")

        obj_name = f"reaction_sandbox_{mod['stem']}.o"
        if obj_name not in dep_text:
            # Add to reaction_sandbox.o dependency list
            for anchor in (
                "  reaction_sandbox_awinhibitmethyl.o\\",
                "  reaction_sandbox_awinhibitmethyl.o \\",
                "  reaction_sandbox_ufd_wp.o \\",
            ):
                if anchor in dep_text:
                    dep_text = dep_text.replace(
                        anchor,
                        anchor + f"\n  {obj_name}\\",
                        1,
                    )
                    break
            else:
                raise RuntimeError(
                    f"Could not find reaction_sandbox.o dep anchor in {dep_path}"
                )

            # Modules that extend BioHill need that object as a dependency.
            src_text = Path(mod["path"]).read_text()
            extra_deps = ""
            if "Reaction_Sandbox_BioHill_class" in src_text:
                extra_deps = "  reaction_sandbox_biohill.o \\\n"
            dep_block = (
                f"{obj_name} : \\\n"
                f"{extra_deps}"
                f"  reaction_sandbox_base.o \\\n"
                f"  reactive_transport_aux.o \\\n"
                f"  global_aux.o \\\n"
                f"  reaction_aux.o\n"
            )
            if dep_block.strip() not in dep_text:
                dep_text = dep_text.rstrip() + "\n" + dep_block + "\n"

        installed.append(mod["file"])

    rs_path.write_text(rs_text)
    obj_path.write_text(obj_text)
    dep_path.write_text(dep_text)
    return installed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pflotran-src",
        required=True,
        help="PFLOTRAN src/pflotran directory (contains reaction_sandbox.F90)",
    )
    parser.add_argument(
        "--sandbox-dir",
        required=True,
        help="Directory containing reaction_sandbox_awinhibit*.F90",
    )
    parser.add_argument(
        "--extra-sandbox-dir",
        action="append",
        default=[],
        help=(
            "Generated custom_* folder with reaction_sandbox_*.F90 "
            "(may be repeated)"
        ),
    )
    args = parser.parse_args()

    pflotran_src = Path(args.pflotran_src)
    sandbox_dir = Path(args.sandbox_dir)

    copy_sandbox_sources(sandbox_dir, pflotran_src)
    patch_reaction_sandbox_f90(pflotran_src / "reaction_sandbox.F90")
    patch_object_files(pflotran_src / "pflotran_object_files.txt")
    patch_dependencies(pflotran_src / "pflotran_dependencies.txt")
    print(f"Patched PFLOTRAN at {pflotran_src}")

    for extra in args.extra_sandbox_dir:
        installed = patch_extra_sandbox(pflotran_src, Path(extra))
        print(f"Installed extra sandbox(es) from {extra}: {', '.join(installed)}")


if __name__ == "__main__":
    main()
