"""Smoke test: the packaged generator is importable as a library."""


def test_generator_importable():
    from pflotran_py.generator import PFLOTRANGenerator
    assert PFLOTRANGenerator is not None
