#!/usr/bin/env python3
"""Patch a PFLOTRAN source tree to include pflotran-py custom reaction sandboxes.

Copies the three AWINHIBIT*.F90 modules from sandbox/ into PFLOTRAN's src
directory and updates reaction_sandbox.F90, pflotran_object_files.txt, and
pflotran_dependencies.txt. Idempotent: safe to run multiple times.

Usage:
    python3 scripts/patch_pflotran_sandboxes.py \\
        --pflotran-src /scratch/pflotran/src/pflotran \\
        --sandbox-dir /work/sandbox
"""

import argparse
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
    args = parser.parse_args()

    pflotran_src = Path(args.pflotran_src)
    sandbox_dir = Path(args.sandbox_dir)

    copy_sandbox_sources(sandbox_dir, pflotran_src)
    patch_reaction_sandbox_f90(pflotran_src / "reaction_sandbox.F90")
    patch_object_files(pflotran_src / "pflotran_object_files.txt")
    patch_dependencies(pflotran_src / "pflotran_dependencies.txt")
    print(f"Patched PFLOTRAN at {pflotran_src}")


if __name__ == "__main__":
    main()
