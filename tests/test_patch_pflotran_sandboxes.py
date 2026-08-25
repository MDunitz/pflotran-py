"""Unit tests for scripts/patch_pflotran_sandboxes.py helpers."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT = _REPO_ROOT / "scripts" / "patch_pflotran_sandboxes.py"


def _load_patch_module():
    spec = importlib.util.spec_from_file_location("patch_pflotran_sandboxes", _SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


patch = _load_patch_module()


def _write_f90(
    directory: Path,
    stem: str,
    *,
    module: str | None = None,
    create: str | None = None,
    extra: str = "",
) -> Path:
    module = module or f"Reaction_Sandbox_{stem.title().replace('_', '')}_class"
    create = create or f"{stem.title().replace('_', '')}Create"
    path = directory / f"reaction_sandbox_{stem}.F90"
    path.write_text(
        f"module {module}\n"
        f"  public :: {create}\n"
        f"{extra}"
        f"end module {module}\n"
    )
    return path


def test_discover_extra_modules_finds_custom_and_formats_fields(tmp_path: Path):
    _write_f90(
        tmp_path,
        "tempbiohill",
        module="Reaction_Sandbox_TempBioHill_class",
        create="TempBioHillCreate",
        extra="  use Reaction_Sandbox_BioHill_class\n",
    )
    _write_f90(
        tmp_path,
        "my_custom",
        module="Reaction_Sandbox_MyCustom_class",
        create="MyCustomCreate",
    )
    # Skipped stems (bundled / template) must not appear in the result.
    _write_f90(tmp_path, "awinhibit")
    _write_f90(tmp_path, "awinhibitacetate")
    _write_f90(tmp_path, "awinhibitmethyl")
    _write_f90(tmp_path, "template")
    _write_f90(tmp_path, "aq")
    (tmp_path / "readme.txt").write_text("ignore me")

    modules = patch._discover_extra_modules(tmp_path)

    assert [m["stem"] for m in modules] == ["my_custom", "tempbiohill"]
    by_stem = {m["stem"]: m for m in modules}

    assert by_stem["tempbiohill"] == {
        "file": "reaction_sandbox_tempbiohill.F90",
        "stem": "tempbiohill",
        "module": "Reaction_Sandbox_TempBioHill_class",
        "create": "TempBioHillCreate",
        "keyword": "TEMPBIOHILL",
        "path": str(tmp_path / "reaction_sandbox_tempbiohill.F90"),
    }
    assert by_stem["my_custom"]["keyword"] == "MYCUSTOM"
    assert by_stem["my_custom"]["module"] == "Reaction_Sandbox_MyCustom_class"
    assert by_stem["my_custom"]["create"] == "MyCustomCreate"


def test_discover_extra_modules_requires_at_least_one(tmp_path: Path):
    _write_f90(tmp_path, "awinhibit")
    with pytest.raises(RuntimeError, match="No reaction_sandbox_"):
        patch._discover_extra_modules(tmp_path)


def test_discover_extra_modules_requires_parseable_symbols(tmp_path: Path):
    bad = tmp_path / "reaction_sandbox_broken.F90"
    bad.write_text("! no module line\n")
    with pytest.raises(RuntimeError, match="Could not parse"):
        patch._discover_extra_modules(tmp_path)


def test_build_module_dep_block_adds_biohill_when_needed():
    plain = patch._build_module_dep_block("reaction_sandbox_foo.o", "module Foo\n")
    assert "reaction_sandbox_biohill.o" not in plain
    assert plain.startswith("reaction_sandbox_foo.o : \\")

    bio = patch._build_module_dep_block(
        "reaction_sandbox_tempbiohill.o",
        "use Reaction_Sandbox_BioHill_class\n",
    )
    assert "reaction_sandbox_biohill.o \\" in bio
    assert patch._needs_biohill_dep("use Reaction_Sandbox_BioHill_class")
    assert not patch._needs_biohill_dep("use Reaction_Sandbox_Base_class")
