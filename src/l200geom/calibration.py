from __future__ import annotations

import logging

import numpy as np
from pyg4ometry import geant4

from . import core, hpge_strings

log = logging.getLogger(__name__)


def place_calibration_system(b: core.InstrumentationData) -> None:
    """Construct LEGEND-200 calibration system."""
    # place calibration tubes.
    if len(b.special_metadata.calibration) == 0:
        return

    calib_tubes = {}
    calib_tube_length = []
    calib_tube_xy = np.empty((2, len(b.special_metadata.calibration)))
    for i, tube in b.special_metadata.calibration.items():
        idx = int(i) - 1
        if tube.length_in_mm not in calib_tubes:
            calib_tubes[tube.length_in_mm] = hpge_strings._get_nylon_mini_shroud(
                tube.tube_radius_in_mm, tube.length_in_mm, True, b.materials, b.registry
            )
            calib_tube_length.append(tube.length_in_mm)

        phi = np.deg2rad(tube.angle_in_deg)
        calib_tube_xy[:, idx] = np.array([tube.radius_in_mm * np.cos(phi), -tube.radius_in_mm * np.sin(phi)])
        nms_pv = geant4.PhysicalVolume(
            [0, 0, 0],
            [*calib_tube_xy[:, idx], b.top_plate_z_pos - tube.length_in_mm / 2],
            calib_tubes[tube.length_in_mm],
            f"calibration_tube_{i}",
            b.mother_lv,
            b.registry,
        )
        hpge_strings._add_nms_surfaces(nms_pv, b.mother_pv, b.materials, b.registry)

    # check if we have one shared length of calibration tubes.
    calib_tube_length = (
        calib_tube_length[0] if all(length == calib_tube_length[0] for length in calib_tube_length) else None
    )
