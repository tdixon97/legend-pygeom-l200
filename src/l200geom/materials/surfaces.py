"""Subpackage to provide all implemented optical surfaces and their properties."""

from __future__ import annotations

import legendoptics.copper
import legendoptics.tetratex
import pyg4ometry.geant4 as g4


class OpticalSurfaceRegistry:
    """Register and define optical surfaces.

    Note on Models
    --------------

    * UNIFIED model:
        `value` is the `sigma_alpha` parameter, the stddev of the newly chosen facet normal direction.
        For details on this model and its parameters, see `UNIFIED model diagram`_.
    * GLISUR model:
        `value` as smoothness, in range [0,1]

    .. _UNIFIED model diagram: https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/_images/UNIFIED_model_diagram.png
    """

    def __init__(self, reg: g4.Registry):
        self.g4_registry = reg
        self._model = "unified"

    @property
    def to_copper(self) -> g4.Material:
        """Reflective surface for copper structure."""
        if hasattr(self, "_to_copper"):
            return self._wlsr_tpb_to_tetratex

        self._wlsr_tpb_to_tetratex = g4.solid.OpticalSurface(
            "surface_to_copper",
            finish="ground",
            model=self._model,
            surf_type="dielectric_metal",
            value=0.9,
            registry=self.g4_registry,
        )

        legendoptics.copper.pyg4_copper_attach_reflectivity(
            self._to_copper,
            self.g4_registry,
        )

        return self._to_copper

    @property
    def wlsr_tpb_to_tetratex(self) -> g4.Material:
        """Reflective surface Tetratex diffuse reflector."""
        if hasattr(self, "_wlsr_tpb_to_tetratex"):
            return self._wlsr_tpb_to_tetratex

        self._wlsr_tpb_to_tetratex = g4.solid.OpticalSurface(
            "surface_wlsr_tpb_to_tetratex",
            finish="groundfrontpainted",
            model=self._model,
            surf_type="dielectric_dielectric",
            value=0.9,
            registry=self.g4_registry,
        )

        legendoptics.tetratex.pyg4_tetratex_attach_reflectivity(
            self._wlsr_tpb_to_tetratex,
            self.g4_registry,
        )

        return self._wlsr_tpb_to_tetratex
