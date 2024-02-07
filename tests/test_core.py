from __future__ import annotations


def test_import():
    import l200geom  # noqa: F401


def test_construct():
    import l200geom.core

    l200geom.core.construct(use_detailed_fiber_model=False)
    l200geom.core.construct(use_detailed_fiber_model=True)
