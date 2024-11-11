from __future__ import annotations


def test_construct_with_sis(tmp_path):
    from l200geom import core

    cfg = {
        "sis": {
            "1": {
                "sis_z": 8250,
                "sources": [None, None, None, "Th228"],
            },
            "2": None,
            "3": None,
            "4": None,
        },
    }

    registry = core.construct(use_detailed_fiber_model=False, config=cfg)
    assert "source_inner" in registry.physicalVolumeDict
