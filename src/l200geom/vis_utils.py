from __future__ import annotations

import logging

import pyg4ometry.geant4 as g4
from pyg4ometry import visualisation

log = logging.getLogger(__name__)


def _color_recursive(lv: g4.LogicalVolume, viewer: visualisation.ViewerBase) -> None:
    if hasattr(lv, "pygeom_color_rgba"):
        for vis in viewer.instanceVisOptions[lv.name]:
            if lv.pygeom_color_rgba is False:
                vis.alpha = 0
                vis.visible = False
            else:
                vis.colour = lv.pygeom_color_rgba[0:3]
                vis.alpha = lv.pygeom_color_rgba[3]
                vis.visible = vis.alpha > 0

    for pv in lv.daughterVolumes:
        if pv.type == "placement":
            _color_recursive(pv.logicalVolume, viewer)


def visualize(registry: g4.Registry) -> None:
    v = visualisation.VtkViewerColouredNew()
    v.addLogicalVolume(registry.worldVolume)

    _color_recursive(registry.worldVolume, v)

    v.buildPipelinesAppend()
    v.addAxes(length=5000)
    v.view()
