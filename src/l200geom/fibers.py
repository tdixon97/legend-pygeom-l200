from __future__ import annotations

import math

from legendmeta.jsondb import AttrsDict, JsonDB
from pyg4ometry import geant4 as g4

from . import materials


def place_fiber_modules(
    fiber_metadata: JsonDB,
    ch_map: AttrsDict,
    mother_lv: g4.LogicalVolume,
    materials: materials.OpticalMaterialRegistry,
    registry: g4.Registry,
    use_detailed_fiber_model: bool = False,
) -> None:
    # Unroll the provided metadata into a structure better suited for the next steps.
    # The geometry here is based on physical modules and not on channels.
    modules = {}
    ch_map = ch_map.map("system", unique=False).spms
    for ch in ch_map.values():
        mod = modules.get(ch.location.fiber)

        if mod is None:
            # initialize a new module object if we don't have one yet.
            tpb_thickness = fiber_metadata[ch.location.fiber].geometry.tpb.thickness_in_nm
            module_type = fiber_metadata[ch.location.fiber].type
            mod = {
                "type": module_type,
                "top": None,
                "bottom": None,
                "tpb_thickness": tpb_thickness,
            }
            modules[ch.location.fiber] = mod

        assert mod[ch.location.position] is None
        mod[ch.location.position] = ch.name

    factory = ModuleFactorySingleFibers if use_detailed_fiber_model else ModuleFactorySegment

    ob_factory = factory(
        radius_mm=590 / 2,
        fiber_length_mm=1900,
        fiber_count_per_module=90,
        number_of_modules=20,
        registry=registry,
        materials=materials,
    )
    ib_factory = factory(
        radius_mm=269 / 2,
        fiber_length_mm=1400,
        fiber_count_per_module=90,
        number_of_modules=9,
        registry=registry,
        materials=materials,
    )

    for mod_name, mod in modules.items():
        m0, m1 = int(mod_name[2:5]), int(mod_name[5:8]) - 1
        assert m0 == m1
        module_num = (m0 - 1) / 2
        if mod["type"] == "outer":
            ob_factory.create_module(mod_name, mod, mother_lv, module_num, z_displacement=+1000)
        if mod["type"] == "inner":
            ib_factory.create_module(mod_name, mod, mother_lv, module_num, z_displacement=+1250)


class ModuleFactoryBase:
    FIBER_DIM = 1

    SIPM_HEIGHT = 1  # mm, dummy
    SIPM_OUTER_EXTRA = 0.2  # mm
    SIPM_GAP = 0.05  # mm
    SIPM_OVERLAP = 0.1  # mm

    def __init__(
        self,
        radius_mm: float,
        fiber_length_mm: float,
        fiber_count_per_module: int,
        number_of_modules: int,
        materials: materials.OpticalMaterialRegistry,
        registry: g4.Registry,
    ):
        self.radius = radius_mm
        self.fiber_length = fiber_length_mm
        self.fiber_count_per_module = fiber_count_per_module
        self.number_of_modules = number_of_modules
        self.materials = materials
        self.registry = registry

    def _cached_sipm_volumes(self) -> None:
        v_suffix = f"_r{self.radius}_nmod{self.number_of_modules}"
        v_name = f"sipm{v_suffix}"
        if v_name in self.registry.solidDict:
            return

        sipm_dim = self.FIBER_DIM + 0.01  # +0.01 to fit round->square
        fiber_segment = 2 * math.pi / self.number_of_modules

        sipm = g4.solid.Tubs(
            v_name,
            self.radius - sipm_dim / 2,
            self.radius + sipm_dim / 2,
            self.SIPM_HEIGHT,
            0,
            fiber_segment,
            self.registry,
            "mm",
        )
        self.sipm_lv = g4.LogicalVolume(sipm, self.materials.metal_silicon, v_name, self.registry)

        sipm_outer1 = g4.solid.Tubs(
            f"sipm_outer1{v_suffix}",
            self.radius - sipm_dim / 2 - self.SIPM_OUTER_EXTRA,
            self.radius + sipm_dim / 2 + self.SIPM_OUTER_EXTRA,
            self.SIPM_HEIGHT + self.SIPM_OUTER_EXTRA + self.SIPM_OVERLAP,
            0,
            fiber_segment,
            self.registry,
            "mm",
        )
        sipm_outer2 = g4.solid.Tubs(
            f"sipm_outer2{v_suffix}",
            self.radius - sipm_dim / 2,
            self.radius + sipm_dim / 2,
            self.SIPM_HEIGHT + 2 * self.SIPM_GAP + self.SIPM_OVERLAP,
            0,
            fiber_segment,
            self.registry,
            "mm",
        )
        sipm_outer_top = g4.solid.Subtraction(
            f"sipm_outer_top{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [0, 0, -self.SIPM_OUTER_EXTRA / 2]],
            self.registry,
        )
        sipm_outer_bottom = g4.solid.Subtraction(
            f"sipm_outer_bottom{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [0, 0, +self.SIPM_OUTER_EXTRA / 2]],
            self.registry,
        )
        self.sipm_outer_top_lv = g4.LogicalVolume(
            sipm_outer_top,
            self.materials.metal_copper,
            f"sipm_outer_top{v_suffix}",
            self.registry,
        )
        self.sipm_outer_bottom_lv = g4.LogicalVolume(
            sipm_outer_bottom,
            self.materials.metal_copper,
            f"sipm_outer_bottom{v_suffix}",
            self.registry,
        )

        sipm_outer_end = g4.solid.Box(
            f"sipm_outer_end{v_suffix}",
            sipm_dim + self.SIPM_OUTER_EXTRA * 2,
            self.SIPM_OUTER_EXTRA,
            self.SIPM_HEIGHT + self.SIPM_OUTER_EXTRA + self.SIPM_OVERLAP,
            self.registry,
        )
        g4.LogicalVolume(
            sipm_outer_end,
            self.materials.metal_copper,
            f"sipm_outer_end{v_suffix}",
            self.registry,
        )

    def _create_sipm(
        self,
        module_num: int,
        fibers: list,
        is_top: bool,
        mother_lv: g4.LogicalVolume,
        sipm_name: str,
        z_displacement: float,
    ) -> None:
        z = +self.fiber_length / 2 + self.SIPM_HEIGHT / 2 + self.SIPM_GAP  # add small gap
        z_outer = (
            z + self.SIPM_OUTER_EXTRA / 2 - self.SIPM_OVERLAP / 2 - self.SIPM_GAP
            if is_top
            else -z - self.SIPM_OUTER_EXTRA / 2 + self.SIPM_OVERLAP / 2 + self.SIPM_GAP
        )
        z = z if is_top else -z
        z += z_displacement
        z_outer += z_displacement
        start_angle = 2 * math.pi / self.number_of_modules * module_num
        g4.PhysicalVolume(
            [0, 0, -start_angle],
            [0, 0, z],
            self.sipm_lv,
            sipm_name,
            mother_lv,
            self.registry,
        )
        g4.PhysicalVolume(
            [0, 0, -start_angle],
            [0, 0, z_outer],
            self.sipm_outer_top_lv if is_top else self.sipm_outer_bottom_lv,
            f"{sipm_name}_wrap",
            mother_lv,
            self.registry,
        )

        # TODO: Add border surfaces to all fibers in this segment.
        # g4.BorderSurface(f'bsurface_lar_sipm{module_num}_{int(is_top)}', lar_phys, sipm_phys,
        #    (mats.surface_sipm_top if is_top else mats.surface_sipm_bottom), self.registry)


class ModuleFactorySingleFibers(ModuleFactoryBase):
    def _cached_fiber_volumes(self) -> None:
        """Create solids, logical and physical volumes for the fibers, as specified by the parameters of this instance."""
        v_suffix = f"_l{self.fiber_length}"
        if f"fiber_cl2{v_suffix}" in self.registry.solidDict:
            return

        fiber_thickness_cl1 = 0.04 * self.FIBER_DIM  # (BCF-91A document)
        fiber_thickness_cl2 = 0.02 * self.FIBER_DIM  # (BCF-91A document)
        # create solids
        self.fiber_cl2 = g4.solid.Box(
            f"fiber_cl2{v_suffix}",
            self.FIBER_DIM,
            self.FIBER_DIM,
            self.fiber_length,
            self.registry,
            "mm",
        )
        dim_cl1 = self.FIBER_DIM - fiber_thickness_cl1
        self.fiber_cl1 = g4.solid.Box(
            f"fiber_cl1{v_suffix}",
            dim_cl1,
            dim_cl1,
            self.fiber_length,
            self.registry,
            "mm",
        )
        dim_core = self.FIBER_DIM - fiber_thickness_cl1 - fiber_thickness_cl2
        self.fiber_core = g4.solid.Box(
            f"fiber_core{v_suffix}",
            dim_core,
            dim_core,
            self.fiber_length,
            self.registry,
            "mm",
        )

        self.fiber_cl2_lv = g4.LogicalVolume(
            self.fiber_cl2,
            self.materials.pmma_out,
            f"fiber_cl2{v_suffix}",
            self.registry,
        )
        self.fiber_cl1_lv = g4.LogicalVolume(
            self.fiber_cl1, self.materials.pmma, f"fiber_cl1{v_suffix}", self.registry
        )
        self.fiber_core_lv = g4.LogicalVolume(
            self.fiber_core,
            self.materials.ps_fibers,
            f"fiber_core{v_suffix}",
            self.registry,
        )

        self.fiber_cl1_pv = g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            self.fiber_cl1_lv,
            f"fiber_cl1{v_suffix}",
            self.fiber_cl2_lv,
            self.registry,
        )
        self.fiber_core_pv = g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            self.fiber_core_lv,
            f"fiber_core{v_suffix}",
            self.fiber_cl1_lv,
            self.registry,
        )

    def _cached_tpb_coating_volume(self, tpb_thickness_nm: float) -> g4.LogicalVolume:
        """Create and cache a TPB coating layer of the specified thickness.

        The TPB-Layer is dependent on the module (i.e. the applied thickness varies slightly),
        so we cannot cache it globally on this instance.
        """
        v_suffix = f"_l{self.fiber_length}_tpb{tpb_thickness_nm}"
        v_name = f"fiber_coating{v_suffix}"
        if v_name in self.registry.solidDict:
            return self.registry.logicalVolumeDict[v_name]

        coating_dim = self.FIBER_DIM + tpb_thickness_nm / 1e6
        coating = g4.solid.Box(v_name, coating_dim, coating_dim, self.fiber_length, self.registry, "mm")
        coating_lv = g4.LogicalVolume(coating, self.materials.tpb_on_fibers, v_name, self.registry)
        g4.PhysicalVolume(
            [0, 0, 0], [0, 0, 0], self.fiber_cl2_lv, f"fiber_cl2{v_suffix}", coating_lv, self.registry
        )

        coating_lv.pygeom_color_rgba = [0, 1, 0, 1]

        return coating_lv

    def create_module(
        self, mod_name: str, mod, mother_lv: g4.LogicalVolume, module_num, z_displacement: float
    ) -> None:
        self._cached_fiber_volumes()
        self._cached_sipm_volumes()
        coating_lv = self._cached_tpb_coating_volume(mod["tpb_thickness"])

        start_angle = 2 * math.pi / self.number_of_modules * module_num
        fibers = []
        for n in range(0, self.fiber_count_per_module):
            th = start_angle + 2 * math.pi / self.number_of_modules / self.fiber_count_per_module * (n + 0.5)
            x = self.radius * math.cos(th)
            y = self.radius * math.sin(th)
            fibers.append(
                g4.PhysicalVolume(
                    [0, 0, -th],
                    [x, y, z_displacement],
                    coating_lv,
                    f"fiber_{mod_name}_{n}",
                    mother_lv,
                    self.registry,
                )
            )

        # create SiPMs and attach to fibers
        self._create_sipm(module_num, fibers, True, mother_lv, mod["top"], z_displacement)
        self._create_sipm(module_num, fibers, False, mother_lv, mod["bottom"], z_displacement)


class ModuleFactorySegment(ModuleFactoryBase):
    def _cached_fiber_volumes(self) -> None:
        """Create solids, logical and physical volumes for the fibers, as specified by the parameters of this instance."""
        v_suffix = f"_l{self.fiber_length}"
        if f"fiber_cl2{v_suffix}" in self.registry.solidDict:
            return

        fiber_thickness_cl1 = 0.04 * self.FIBER_DIM  # (BCF-91A document)
        fiber_thickness_cl2 = 0.02 * self.FIBER_DIM  # (BCF-91A document)
        # create solids
        angle = 2 * math.pi / self.number_of_modules
        self.fiber_cl2 = g4.solid.Tubs(
            f"fiber_cl2{v_suffix}",
            self.radius - self.FIBER_DIM,
            self.radius,
            self.fiber_length,
            0,
            angle,
            self.registry,
            "mm",
        )
        dim_cl1 = self.FIBER_DIM - fiber_thickness_cl1
        self.fiber_cl1 = g4.solid.Tubs(
            f"fiber_cl1{v_suffix}",
            self.radius - dim_cl1,
            self.radius + dim_cl1 - self.FIBER_DIM,
            self.fiber_length,
            0,
            angle,
            self.registry,
            "mm",
        )
        dim_core = self.FIBER_DIM - fiber_thickness_cl1 - fiber_thickness_cl2
        self.fiber_core = g4.solid.Tubs(
            f"fiber_core{v_suffix}",
            self.radius - dim_core,
            self.radius + dim_core - self.FIBER_DIM,
            self.fiber_length,
            0,
            angle,
            self.registry,
            "mm",
        )

        self.fiber_cl2_lv = g4.LogicalVolume(
            self.fiber_cl2,
            self.materials.pmma_out,
            f"fiber_cl2{v_suffix}",
            self.registry,
        )
        self.fiber_cl1_lv = g4.LogicalVolume(
            self.fiber_cl1, self.materials.pmma, f"fiber_cl1{v_suffix}", self.registry
        )
        self.fiber_core_lv = g4.LogicalVolume(
            self.fiber_core,
            self.materials.ps_fibers,
            f"fiber_core{v_suffix}",
            self.registry,
        )

        self.fiber_cl1_pv = g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            self.fiber_cl1_lv,
            f"fiber_cl1{v_suffix}",
            self.fiber_cl2_lv,
            self.registry,
        )
        self.fiber_core_pv = g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            self.fiber_core_lv,
            f"fiber_core{v_suffix}",
            self.fiber_cl1_lv,
            self.registry,
        )

    def _cached_tpb_coating_volume(self, tpb_thickness_nm: float) -> g4.LogicalVolume:
        """Create and cache a TPB coating layer of the specified thickness.

        The TPB-Layer is dependent on the module (i.e. the applied thickness varies slightly),
        so we cannot cache it globally on this instance.
        """
        v_suffix = f"_l{self.fiber_length}_tpb{tpb_thickness_nm}"
        v_name = f"fiber_coating{v_suffix}"
        if v_name in self.registry.solidDict:
            return self.registry.logicalVolumeDict[v_name]

        coating_dim = self.FIBER_DIM + tpb_thickness_nm / 1e6
        coating = g4.solid.Tubs(
            v_name,
            self.radius - coating_dim,
            self.radius + coating_dim - self.FIBER_DIM,
            self.fiber_length,
            0,
            2 * math.pi / self.number_of_modules,
            self.registry,
            "mm",
        )
        coating_lv = g4.LogicalVolume(coating, self.materials.tpb_on_fibers, v_name, self.registry)
        g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0],
            self.fiber_cl2_lv,
            f"fiber_cl2{v_suffix}",
            coating_lv,
            self.registry,
        )

        coating_lv.pygeom_color_rgba = [0, 1, 0, 1]

        return coating_lv

    def create_module(
        self,
        mod_name: str,
        mod,
        mother_lv: g4.LogicalVolume,
        module_num,
        z_displacement: float,
    ) -> None:
        self._cached_fiber_volumes()
        self._cached_sipm_volumes()
        coating_lv = self._cached_tpb_coating_volume(mod["tpb_thickness"])

        start_angle = 2 * math.pi / self.number_of_modules * (module_num)
        fibers = []

        th = start_angle
        fibers.append(
            g4.PhysicalVolume(
                [0, 0, -th],
                [0, 0, z_displacement],
                coating_lv,
                f"fiber_{mod_name}_s",
                mother_lv,
                self.registry,
            )
        )

        # create SiPMs and attach to fibers
        self._create_sipm(module_num, fibers, True, mother_lv, mod["top"], z_displacement)
        self._create_sipm(module_num, fibers, False, mother_lv, mod["bottom"], z_displacement)
