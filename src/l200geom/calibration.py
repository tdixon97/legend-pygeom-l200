from __future__ import annotations

import logging
import math

import numpy as np
from pyg4ometry import geant4

from . import core, hpge_strings, materials

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


def place_ta_absorber(
    z0: float,
    mothervolume: geant4.LogicalVolume,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> None:
    """Construct LEGEND-200 HPGe strings.

    Parameters
    ----------
    z0
        The position of the top most slot of the strings in the z-axis.
    mothervolume
        pyg4ometry Geant4 LogicalVolume instance in which the strings
        are to be placed.
    registry
        pyg4ometry Geant4 registry instance.
    """
    if registry is None:
        msg = "registry cannot be None"
        raise ValueError(msg)

    source_height = 17.6  # mm
    source_outside = 10.6  # mm
    source_inside = source_height - source_outside
    source_radius = 6.4 / 2
    absorber_height = 37.5
    ta_absorber_outer = geant4.solid.Tubs(
        "ta_absorber_outer",
        0,
        16.2, # estimate!
        absorber_height,
        0,
        2 * math.pi,
        registry
    )
    ta_absorber_inner = geant4.solid.Tubs(
        "ta_absorber_inner",
        0,
        source_radius,
        source_height - source_outside,
        0,
        2 * math.pi,
        registry
    )
    ta_absorber = geant4.solid.Subtraction("ta_absorber", ta_absorber_outer, ta_absorber_inner, [[0, 0, 0], [0, 0, (absorber_height - source_inside)/2]], registry)
    ta_absorber_lv =geant4.LogicalVolume(ta_absorber, materials.metal_tantalum, "ta_absorber", registry)

    # all positions from MaGe, might be incorrect!
    geant4.PhysicalVolume(
        [0, 0, 0],
        # [-121.304, 96.48977, calib_tube_z],
        [-118.4352, 93.870075, z0],  # to fix overlaps, not the physical position.
        ta_absorber_lv,
        "ta_absorber_3",
        mothervolume,
        registry,
    )

    source_outer = geant4.solid.Tubs(
        "source_outer",
        0,
        source_radius,
        source_height,
        0,
        2 * math.pi,
        registry
    )
    source_outer = geant4.LogicalVolume(source_outer, materials.metal_steel, "source_outer", registry)
    source_z = z0 + absorber_height/2 + source_height/2 - source_inside
    geant4.PhysicalVolume(
        [0, 0, 0],
        # [-121.304, 96.48977, calib_tube_z],
        [-118.4352, 93.870075, source_z],  # to fix overlaps, not the physical position.
        source_outer,
        "source_outer",
        mothervolume,
        registry,
    )
    source_inner = geant4.solid.Tubs(
        "source_inner",
        0,
        2,
        4,
        0,
        2 * math.pi,
        registry
    )
    source_inner = geant4.LogicalVolume(source_inner, materials.metal_caps_gold, "source_inner", registry)
    source_inner.pygeom_color_rgba = (1, 1, 0, 1)
    source_inner_z = source_height/2 - 4/2 - 1.2
    geant4.PhysicalVolume(
        [0, 0, 0],
        # [-121.304, 96.48977, calib_tube_z],
        [0, 0, source_inner_z],  # to fix overlaps, not the physical position.
        source_inner,
        "source_inner",
        source_outer,
        registry,
    )

    cu_aborber_thickness = 2  # mm
    cu_absorber_outer = geant4.solid.Tubs(
        "cu_absorber_outer",
        0,
        source_radius + cu_aborber_thickness,
        source_outside + cu_aborber_thickness,
        0,
        2 * math.pi,
        registry
    )
    cu_absorber_inner = geant4.solid.Tubs(
        "cu_absorber_inner",
        0,
        source_radius,
        source_outside,
        0,
        2 * math.pi,
        registry
    )
    cu_absorber = geant4.solid.Subtraction("cu_absorber", cu_absorber_outer, cu_absorber_inner, [[0, 0, 0], [0, 0, -cu_aborber_thickness/2]], registry)
    cu_absorber = geant4.LogicalVolume(cu_absorber, materials.metal_copper, "cu_absorber", registry)
    cu_absorber_z = z0 + absorber_height/2 + cu_aborber_thickness/2 + source_outside/2
    if False:
        geant4.PhysicalVolume(
            [0, 0, 0],
            # [-121.304, 96.48977, calib_tube_z],
            [-118.4352, 93.870075, cu_absorber_z],  # to fix overlaps, not the physical position.
            cu_absorber,
            "cu_absorber",
            mothervolume,
            registry,
        )

    peek_outside = geant4.solid.Box("peek_outside", 33.1, 9, 25, registry)
    peek_inside = geant4.solid.Box("peek_inside", 14.5, 9, 15, registry)
    peek_holder = geant4.solid.Subtraction("peek_holder", peek_outside, peek_inside, [[0, 0, 0], [0, 0, -10/2]], registry)
    peek_holder = geant4.LogicalVolume(peek_holder, materials.peek, "peek_holder", registry)
    peek_holder_z = z0 + absorber_height/2 + 25/2
    geant4.PhysicalVolume(
        [0, 0, 0],
        # [-121.304, 96.48977, calib_tube_z],
        [-118.4352, 93.870075, peek_holder_z],  # to fix overlaps, not the physical position.
        peek_holder,
        "peek_holder",
        mothervolume,
        registry,
    )
