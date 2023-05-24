from __future__ import annotations

from pyg4ometry import geant4


def construct() -> geant4.Registry:
    """Construct the LEGEND-200 geometry and return the pyg4ometry Registry containing the world volume."""

    reg = geant4.Registry()

    # Create the world volume
    world_material = geant4.MaterialPredefined('G4_Galactic')
    world = geant4.solid.Box('world', 20, 20, 20, reg, 'm')
    world_lv = geant4.LogicalVolume(world, world_material, 'world', reg)
    reg.setWorld(world_lv)

    return reg

