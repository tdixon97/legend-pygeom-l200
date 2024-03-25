from __future__ import annotations

import json
import math
from importlib import resources
from pathlib import Path

import pyg4ometry
from legendhpges import make_hpge
from legendmeta import AttrsDict
from pyg4ometry import geant4

from . import materials


def place_hpge_strings(
    channelmap: str | dict | AttrsDict,
    string_config: str | dict | AttrsDict,
    z0: float,
    mothervolume: geant4.LogicalVolume,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> None:
    """Construct LEGEND-200 HPGe strings.

    Parameters
    ----------
    channelmap
        LEGEND-200 HPGe channel map containing germanium detectors
        configuration in the string and their geometry.
    string_config
        LEGEND-200 germanium detector string configuration file.
        Used to reconstruct the spatial position of each string.
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

    if channelmap is None:
        msg = "configuration metadata file cannot be None"
        raise ValueError(msg)

    if string_config is None:
        msg = "string configuration cannot be None"
        raise ValueError(msg)

    if not isinstance(channelmap, (dict, AttrsDict)):
        with Path(channelmap).open() as jfile:
            ch_map = AttrsDict(json.load(jfile))
    else:
        ch_map = AttrsDict(channelmap)

    if not isinstance(string_config, (dict, AttrsDict)):
        with Path(string_config).open() as jfile:
            hpge_string_config = AttrsDict(json.load(jfile))
    else:
        hpge_string_config = AttrsDict(string_config)

    ch_map = ch_map.map("system", unique=False).geds.values()

    for hpge_meta in ch_map:
        # Temporary fix for gedet with null enrichment value
        if hpge_meta.production.enrichment is None:
            continue

        hpge_string_id = str(hpge_meta.location.string)
        hpge_string = hpge_string_config.hpge_string[hpge_string_id]
        hpge_unit_id_in_string = hpge_meta.location.position

        x_pos = hpge_string.radius_in_mm * math.cos(math.pi * hpge_string.angle_in_deg / 180)

        y_pos = hpge_string.radius_in_mm * math.sin(math.pi * hpge_string.angle_in_deg / 180)

        z_pos = (
            z0
            - sum(hpge_string.hpge_unit_heights_in_mm[: hpge_unit_id_in_string - 1])
            - hpge_string.hpge_unit_heights_in_mm[hpge_unit_id_in_string - 1] / 2
        )

        hpge = make_hpge(hpge_meta, registry)

        geant4.PhysicalVolume(
            [
                0,
                0,
                0,
            ],
            [x_pos, y_pos, z_pos],
            hpge,
            hpge_meta.name,
            mothervolume,
            registry,
        )


_pen_plate_cache = {}
_minishroud_cache = {}


def _get_pen_plate(
    size: str,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> geant4.LogicalVolume:
    if size not in ["small", "medium", "large", "xlarge"]:
        msg = f"Invalid PEN-plate size {size}"
        raise ValueError(msg)

    if size not in _pen_plate_cache:
        pen_file = resources.files("l200geom") / "models" / f"BasePlate_{size}.stl"
        pen_solid = pyg4ometry.stl.Reader(
            pen_file, solidname=f"pen_{size}", centre=False, registry=registry
        ).getSolid()
        _pen_plate_cache[size] = geant4.LogicalVolume(pen_solid, materials.pen, f"pen_{size}", registry)

    return _pen_plate_cache[size]


def _get_nylon_mini_shroud(
    radius: int,
    length: int,
    materials: materials.OpticalMaterialRegistry,
    registry: geant4.Registry,
) -> geant4.LogicalVolume:
    """Create a nylon/TPB funnel of the given outer dimensions, which will be closed at the top/bottom.

    .. note:: this can also be used for calibration tubes.
    """
    shroud_name = f"minishroud_{radius}x{length}"
    if shroud_name not in _minishroud_cache:
        shroudThickness = 0.1  # mm
        outer = geant4.solid.Tubs(f"{shroud_name}_outer", 0, radius, length, 0, 2 * math.pi, registry)
        inner = geant4.solid.Tubs(
            f"{shroud_name}_inner",
            0,
            radius - shroudThickness,
            length - 2 * shroudThickness,
            0,
            2 * math.pi,
            registry,
        )
        # subtract the slightly smaller solid from the larger one, to get a hollow and closed volume.
        shroud = geant4.solid.Subtraction("minishroud", outer, inner, [[0, 0, 0], [0, 0, 0]], registry)
        _minishroud_cache[shroud_name] = geant4.LogicalVolume(
            shroud,
            materials.tpb_on_nylon,
            shroud_name,
            registry,
        )

    # TODO: implement optical surfaces
    return _minishroud_cache[shroud_name]
