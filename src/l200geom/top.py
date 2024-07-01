from __future__ import annotations

import logging
from importlib import resources

import pyg4ometry
from pyg4ometry import geant4

from . import materials

log = logging.getLogger(__name__)


def place_top_plate(
    z0: float,
    mothervolume: geant4.LogicalVolume,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> None:
    """Construct LEGEND-200 HPGe strings.

    Parameters
    ----------
    z0
        The z coordinate of the top face of the array top plate.
    mothervolume
        pyg4ometry Geant4 LogicalVolume instance in which the strings
        are to be placed.
    registry
        pyg4ometry Geant4 registry instance.
    """
    if registry is None:
        msg = "registry cannot be None"
        raise ValueError(msg)

    plate_file = resources.files("l200geom") / "models" / "TopPlate.stl"
    plate = pyg4ometry.stl.Reader(
        plate_file, solidname="top_plate", centre=False, registry=registry
    ).getSolid()
    plate = geant4.LogicalVolume(plate, materials.metal_copper, "top_plate", registry)
    plate.pygeom_color_rgba = (0.72, 0.45, 0.2, 0.2)

    geant4.PhysicalVolume([0, 0, 0], [0, 0, z0], plate, plate.name, mothervolume, registry)
