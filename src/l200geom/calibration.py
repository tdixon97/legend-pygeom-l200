from __future__ import annotations

import logging

import numpy as np
from pyg4ometry import geant4

from . import core, hpge_strings

log = logging.getLogger(__name__)


def place_calibration_system(b: core.InstrumentationData) -> None:
    """Construct LEGEND-200 calibration system."""
    # place calibration tubes.
    # note: the length is just a rough guess from MaGe, the radius is from slides of A. Lubashevskiy.
    calib_tube_length = 1400
    calib_tube = hpge_strings._get_nylon_mini_shroud(19, calib_tube_length, True, b.materials, b.registry)
    calib_tube_z = b.top_plate_z_pos - calib_tube_length / 2

    # all positions from CAD model.
    calib_tube_r = 155  # mm
    calib_tube_phi = np.deg2rad(np.array([158.57, 261.43, 338.57, 81.43]))
    calib_tube_xy = np.array([calib_tube_r * np.cos(calib_tube_phi), -calib_tube_r * np.sin(calib_tube_phi)])
    for i in range(calib_tube_phi.shape[0]):
        nms_pv = geant4.PhysicalVolume(
            [0, 0, 0],
            [*calib_tube_xy[:, i], calib_tube_z],
            calib_tube,
            f"calibration_tube_{i+1}",
            b.mother_lv,
            b.registry,
        )
        hpge_strings._add_nms_surfaces(nms_pv, b.mother_pv, b.materials, b.registry)
