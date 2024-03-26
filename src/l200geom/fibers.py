from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from legendmeta import AttrsDict, TextDB
from pyg4ometry import geant4 as g4
from scipy.spatial.transform import Rotation

from . import materials


def place_fiber_modules(
    fiber_metadata: TextDB,
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
            mod = FiberModuleData(
                barrel=fiber_metadata[ch.location.fiber].type,
                tpb_thickness=fiber_metadata[ch.location.fiber].geometry.tpb.thickness_in_nm,
            )
            modules[ch.location.fiber] = mod

        assert getattr(mod, f"channel_{ch.location.position}") is None
        setattr(mod, f"channel_{ch.location.position}", ch.name)

    factory = ModuleFactorySingleFibers if use_detailed_fiber_model else ModuleFactorySegment

    # note: actually the radius is only 150mm and another short straight segment of 50mm is following after
    # the bend. to simplify things here, those two are combined to one bent shape, to have at least the same
    # covered solid angle.
    ob_radius = 150 + 50
    ob_factory = factory(
        radius_mm=590 / 2,
        fiber_length_mm=1900 - math.pi * ob_radius / 2,
        bend_radius_mm=ob_radius,
        fiber_count_per_module=81,
        number_of_modules=20,
        registry=registry,
        materials=materials,
    )
    ib_factory = factory(
        radius_mm=269 / 2,
        fiber_length_mm=1400,
        bend_radius_mm=None,
        fiber_count_per_module=81,
        number_of_modules=9,
        registry=registry,
        materials=materials,
    )

    for mod_name, mod in modules.items():
        m0, m1 = int(mod_name[2:5]), int(mod_name[5:8]) - 1
        assert m0 == m1
        module_num = int((m0 - 1) / 2)
        if mod.barrel == "outer":
            ob_factory.create_module(mod_name, mod, mother_lv, module_num, z_displacement=+1000)
        if mod.barrel == "inner":
            ib_factory.create_module(mod_name, mod, mother_lv, module_num, z_displacement=+1250)


@dataclass
class FiberModuleData:
    barrel: str
    tpb_thickness: float
    channel_top: str | None = None
    channel_bottom: str | None = None


class ModuleFactoryBase(ABC):
    FIBER_DIM = 1  # mm
    FIBER_THICKNESS_CL1 = 0.04 * FIBER_DIM  # (BCF-91A document)
    FIBER_THICKNESS_CL2 = 0.02 * FIBER_DIM  # (BCF-91A document)

    # The "SiPM" is only a dummy implementation here, the fibers are not bent to the real geometry of 3x3
    # fibers coupled to one physical SiPM!
    SIPM_HEIGHT = 1  # mm, dummy
    # There is a LAr gap between the fiber end and SiPM:
    SIPM_GAP = 0.05  # mm
    # Because of this gap, we cannot use surfaces between fibers and the "SiPM". To stop stray light from
    # entering the "SiPM" from the other sides, we add an outer envelope that blocks light. On each outer
    # side of the SiPM volume, an additional solid of size "outer extra" is added:
    SIPM_OUTER_EXTRA = 0.2  # mm
    # To also stop stray light from directions more close to the fibers, the envelope extends a bit more
    # along the fibers:
    SIPM_OVERLAP = 0.1  # mm
    SIPM_GAP_SIDE = 0.01  # mm, for fitting problems with round "SiPMs" and square fibers.

    def __init__(
        self,
        radius_mm: float,
        fiber_length_mm: float,
        fiber_count_per_module: int,
        bend_radius_mm: float | None,
        number_of_modules: int,
        materials: materials.OpticalMaterialRegistry,
        registry: g4.Registry,
    ):
        """
        Create a fiber module factory.

        Parameters
        ----------
        radius_mm
            radius of the fiber barrel
        fiber_length_mm
            length of the straight section of this fiber module
        fiber_count_per_module
            number of single fibers per module
        bend_radius_mm
            radius of the bottom bend, or None if the fibers are not bent at the bottom end.
        number_of_modules
            number of modules that cover the full circle
        """
        self.radius = radius_mm
        self.fiber_length = fiber_length_mm
        self.fiber_count_per_module = fiber_count_per_module
        self.bend_radius_mm = bend_radius_mm
        self.number_of_modules = number_of_modules
        self.materials = materials
        self.registry = registry

    def _cached_sipm_volumes(self) -> None:
        """Creates (dummy) SiPM volumes for use at the top/bottom of straight fiber sections."""
        v_suffix = f"_r{self.radius}_nmod{self.number_of_modules}"
        v_name = f"sipm{v_suffix}"
        if v_name in self.registry.solidDict:
            return

        sipm_dim = self.FIBER_DIM + self.SIPM_GAP_SIDE  # GAP_SIDE to fit round->square
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

        # TODO: implement partial modules with end envelopes for SiPM.
        # sipm_outer_end = g4.solid.Box(
        #    f"sipm_outer_end{v_suffix}",
        #    sipm_dim + self.SIPM_OUTER_EXTRA * 2,
        #    self.SIPM_OUTER_EXTRA,
        #    self.SIPM_HEIGHT + self.SIPM_OUTER_EXTRA + self.SIPM_OVERLAP,
        #    self.registry,
        # )
        # g4.LogicalVolume(
        #    sipm_outer_end,
        #    self.materials.metal_copper,
        #    f"sipm_outer_end{v_suffix}",
        #    self.registry,
        # )

    @abstractmethod
    def create_module(
        self,
        mod_name: str,
        mod: FiberModuleData,
        mother_lv: g4.LogicalVolume,
        module_num: int,
        z_displacement: float,
    ) -> None:
        raise NotImplementedError()

    def _create_sipm(
        self,
        module_num: int,
        fibers: list[g4.PhysicalVolume],
        is_top: bool,
        mother_lv: g4.LogicalVolume,
        sipm_name: str,
        z_displacement: float,
    ) -> None:
        """Creates a (dummy) SiPM physical volume for use at the top/bottom of straight fiber sections."""
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
    def _cached_sipm_volumes_bend(self) -> None:
        """Creates (dummy) SiPM volumes for use at the bottom of bent fiber sections."""
        v_suffix = "_bend"  # this dummy SiPM is for one single fiber, so we do not need any attributes here.
        v_name = f"sipm{v_suffix}"
        if v_name in self.registry.solidDict:
            return

        sipm_dim = self.FIBER_DIM + self.SIPM_GAP_SIDE  # GAP_SIDE to fit round->square

        sipm = g4.solid.Box(
            v_name,
            self.SIPM_HEIGHT,
            sipm_dim,
            sipm_dim,
            self.registry,
            "mm",
        )
        self.sipm_lv_bend = g4.LogicalVolume(sipm, self.materials.metal_silicon, v_name, self.registry)

        sipm_outer1 = g4.solid.Box(
            f"sipm_outer1{v_suffix}",
            self.SIPM_HEIGHT + self.SIPM_OUTER_EXTRA + self.SIPM_OVERLAP,
            sipm_dim + 2 * self.SIPM_OUTER_EXTRA,
            sipm_dim + 2 * self.SIPM_OUTER_EXTRA,
            self.registry,
            "mm",
        )
        sipm_outer2 = g4.solid.Box(
            f"sipm_outer2{v_suffix}",
            self.SIPM_HEIGHT + 2 * self.SIPM_GAP + self.SIPM_OVERLAP,
            sipm_dim,
            sipm_dim,
            self.registry,
            "mm",
        )
        sipm_outer_bottom = g4.solid.Subtraction(
            f"sipm_outer_bottom{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [-self.SIPM_OUTER_EXTRA / 2, 0, 0]],
            self.registry,
        )
        self.sipm_outer_bottom_lv_bend = g4.LogicalVolume(
            sipm_outer_bottom,
            self.materials.metal_copper,
            f"sipm_outer_bottom{v_suffix}",
            self.registry,
        )

    def _cached_fiber_volumes(self) -> None:
        """Create solids, logical and physical volumes for the fibers, as specified by the parameters of this instance."""
        v_suffix = f"_l{self.fiber_length}_b{self.bend_radius_mm}"
        if f"fiber_cl2{v_suffix}" in self.registry.solidDict:
            return

        # create solids
        self.fiber_cl2 = g4.solid.Box(
            f"fiber_cl2{v_suffix}",
            self.FIBER_DIM,
            self.FIBER_DIM,
            self.fiber_length,
            self.registry,
            "mm",
        )
        dim_cl1 = self.FIBER_DIM - self.FIBER_THICKNESS_CL1
        self.fiber_cl1 = g4.solid.Box(
            f"fiber_cl1{v_suffix}",
            dim_cl1,
            dim_cl1,
            self.fiber_length,
            self.registry,
            "mm",
        )
        dim_core = self.FIBER_DIM - self.FIBER_THICKNESS_CL1 - self.FIBER_THICKNESS_CL2
        self.fiber_core = g4.solid.Box(
            f"fiber_core{v_suffix}",
            dim_core,
            dim_core,
            self.fiber_length,
            self.registry,
            "mm",
        )
        if self.bend_radius_mm is not None:
            self.fiber_cl2_bend = g4.solid.Tubs(
                f"fiber_cl2_bend{v_suffix}",
                self.bend_radius_mm - self.FIBER_DIM / 2,
                self.bend_radius_mm + self.FIBER_DIM / 2,
                self.FIBER_DIM,
                0,
                math.pi / 2,
                self.registry,
                "mm",
            )
            self.fiber_cl1_bend = g4.solid.Tubs(
                f"fiber_cl1_bend{v_suffix}",
                self.bend_radius_mm - dim_cl1 / 2,
                self.bend_radius_mm + dim_cl1 / 2,
                dim_cl1,
                0,
                math.pi / 2,
                self.registry,
                "mm",
            )
            self.fiber_core_bend = g4.solid.Tubs(
                f"fiber_core_bend{v_suffix}",
                self.bend_radius_mm - dim_core / 2,
                self.bend_radius_mm + dim_core / 2,
                dim_core,
                0,
                math.pi / 2,
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
        if self.bend_radius_mm is not None:
            self.fiber_cl2_bend_lv = g4.LogicalVolume(
                self.fiber_cl2_bend,
                self.materials.pmma_out,
                f"fiber_cl2_bend{v_suffix}",
                self.registry,
            )
            self.fiber_cl1_bend_lv = g4.LogicalVolume(
                self.fiber_cl1_bend, self.materials.pmma, f"fiber_cl1_bend{v_suffix}", self.registry
            )
            self.fiber_core_bend_lv = g4.LogicalVolume(
                self.fiber_core_bend,
                self.materials.ps_fibers,
                f"fiber_core_bend{v_suffix}",
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
        if self.bend_radius_mm is not None:
            self.fiber_cl1_bend_pv = g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                self.fiber_cl1_bend_lv,
                f"fiber_cl1_bend{v_suffix}",
                self.fiber_cl2_bend_lv,
                self.registry,
            )
            self.fiber_core_bend_pv = g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, 0],
                self.fiber_core_bend_lv,
                f"fiber_core_bend{v_suffix}",
                self.fiber_cl1_bend_lv,
                self.registry,
            )

    def _cached_tpb_coating_volume(self, tpb_thickness_nm: float, bend: bool = False) -> g4.LogicalVolume:
        """Create and cache a TPB coating layer of the specified thickness.

        The TPB-Layer is dependent on the module (i.e. the applied thickness varies slightly),
        so we cannot cache it globally on this instance.
        """
        v_suffix = f"{'_bend' if bend else ''}_l{self.fiber_length}_tpb{tpb_thickness_nm}"
        v_name = f"fiber_coating{v_suffix}"
        if v_name in self.registry.solidDict:
            return self.registry.logicalVolumeDict[v_name]

        coating_dim = self.FIBER_DIM + tpb_thickness_nm / 1e6
        if not bend:
            coating = g4.solid.Box(v_name, coating_dim, coating_dim, self.fiber_length, self.registry, "mm")
            inner_lv = self.fiber_cl2_lv
        else:
            coating = g4.solid.Tubs(
                v_name,
                self.bend_radius_mm - coating_dim / 2,
                self.bend_radius_mm + coating_dim / 2,
                coating_dim,
                0,
                math.pi / 2,
                self.registry,
                "mm",
            )
            inner_lv = self.fiber_cl2_bend_lv
        coating_lv = g4.LogicalVolume(coating, self.materials.tpb_on_fibers, v_name, self.registry)
        g4.PhysicalVolume([0, 0, 0], [0, 0, 0], inner_lv, f"fiber_cl2{v_suffix}", coating_lv, self.registry)

        coating_lv.pygeom_color_rgba = [0, 1, 0, 1]

        return coating_lv

    def create_module(
        self,
        mod_name: str,
        mod: FiberModuleData,
        mother_lv: g4.LogicalVolume,
        module_num: int,
        z_displacement: float,
    ) -> None:
        if module_num < 0 or module_num >= self.number_of_modules:
            msg = f"invalid module number {module_num} for a maximum of {self.number_of_modules}-1 modules."
            raise ValueError(msg)

        self._cached_fiber_volumes()
        self._cached_sipm_volumes()
        coating_lv = self._cached_tpb_coating_volume(mod.tpb_thickness, bend=False)
        if self.bend_radius_mm is not None:
            self._cached_sipm_volumes_bend()
            coating_lv_bend = self._cached_tpb_coating_volume(mod.tpb_thickness, bend=True)

        start_angle = 2 * math.pi / self.number_of_modules * module_num
        fibers = []
        for n in range(self.fiber_count_per_module):
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
            if self.bend_radius_mm is not None:
                rotvec = Rotation.from_euler("XZ", [math.pi / 2, -th]).as_euler("xyz")
                x2 = (self.radius - self.bend_radius_mm) * math.cos(th)
                y2 = (self.radius - self.bend_radius_mm) * math.sin(th)
                fibers.append(
                    g4.PhysicalVolume(
                        list(rotvec),
                        [x2, y2, z_displacement - self.fiber_length / 2],
                        coating_lv_bend,
                        f"fiber_bend_{mod_name}_{n}",
                        mother_lv,
                        self.registry,
                    )
                )
                # TODO: for bent modules, implement something not to have overlaps

                # add per-fiber SiPMs (I do not know any other way...)
                sipm_placement_r = self.radius - self.bend_radius_mm - self.SIPM_GAP - self.SIPM_HEIGHT / 2
                x2 = sipm_placement_r * math.cos(th)
                y2 = sipm_placement_r * math.sin(th)
                z = z_displacement - self.fiber_length / 2 - self.bend_radius_mm
                g4.PhysicalVolume(
                    [0, 0, -th],
                    [x2, y2, z],
                    self.sipm_lv_bend,
                    f"{mod.channel_bottom}_{n}",
                    mother_lv,
                    self.registry,
                )

                sipm_placement_outer_r = (
                    sipm_placement_r + self.SIPM_OUTER_EXTRA / 2 - self.SIPM_OVERLAP / 2 - self.SIPM_GAP
                )
                x2 = sipm_placement_outer_r * math.cos(th)
                y2 = sipm_placement_outer_r * math.sin(th)
                g4.PhysicalVolume(
                    [0, 0, -th],
                    [x2, y2, z],
                    self.sipm_outer_bottom_lv_bend,
                    f"{mod.channel_bottom}_{n}_wrap",
                    mother_lv,
                    self.registry,
                )

        # create SiPMs and attach to fibers
        self._create_sipm(module_num, fibers, True, mother_lv, mod.channel_top, z_displacement)
        if self.bend_radius_mm is None:
            self._create_sipm(module_num, fibers, False, mother_lv, mod.channel_bottom, z_displacement)


class ModuleFactorySegment(ModuleFactoryBase):
    def _get_bend_polycone(
        self, inner_r: float, outer_r: float
    ) -> tuple[np.typing.ArrayLike, np.typing.ArrayLike]:
        """In the segmented model, there is no fundamental shape for the fiber bent available, so we
        use a polycone as a replacement."""
        delta_r_mm = (outer_r - inner_r) / 2
        bend_r_outer = self.bend_radius_mm + delta_r_mm
        bend_r_inner = self.bend_radius_mm - delta_r_mm

        angles = np.linspace(0, np.pi / 2, 100)
        z1 = bend_r_outer * np.sin(angles)
        r1 = bend_r_outer * np.cos(angles)
        z2 = bend_r_inner * np.sin(angles)
        r2 = bend_r_inner * np.cos(angles)

        z = self.bend_radius_mm - np.concatenate((z1, np.flip(z2)))
        # offset by the radius at the inner end of the bend.
        r = (outer_r - bend_r_outer) + np.concatenate((r1, np.flip(r2)))

        return z, r

    def _cached_sipm_volumes_bend(self) -> None:
        """Creates (dummy) SiPM volumes for use at the bottom of bent fiber sections."""
        v_suffix = f"_bend{self.bend_radius_mm}_r{self.radius}_nmod{self.number_of_modules}"
        v_name = f"sipm{v_suffix}"
        if v_name in self.registry.solidDict:
            return

        sipm_dim = self.FIBER_DIM + self.SIPM_GAP_SIDE  # GAP_SIDE to fit round->square
        fiber_segment = 2 * math.pi / self.number_of_modules
        # radius of the inner circle at the bottom, already including the small gap between fibers and SiPM.
        inner_radius = self.radius - self.bend_radius_mm - self.SIPM_GAP

        sipm = g4.solid.Tubs(
            v_name,
            inner_radius - self.SIPM_HEIGHT,
            inner_radius,
            sipm_dim,
            0,
            fiber_segment,
            self.registry,
            "mm",
        )
        self.sipm_lv_bend = g4.LogicalVolume(sipm, self.materials.metal_silicon, v_name, self.registry)

        sipm_outer1 = g4.solid.Tubs(
            f"sipm_outer1{v_suffix}",
            inner_radius - self.SIPM_HEIGHT - self.SIPM_OUTER_EXTRA,
            inner_radius + self.SIPM_OVERLAP,
            sipm_dim + 2 * self.SIPM_OUTER_EXTRA,
            0,
            fiber_segment,
            self.registry,
            "mm",
        )
        sipm_outer2 = g4.solid.Tubs(
            f"sipm_outer2{v_suffix}",
            inner_radius - self.SIPM_HEIGHT,
            inner_radius + 2 * self.SIPM_GAP + self.SIPM_OVERLAP,
            sipm_dim,
            0,
            fiber_segment,
            self.registry,
            "mm",
        )
        sipm_outer_bottom = g4.solid.Subtraction(
            f"sipm_outer_bottom{v_suffix}",
            sipm_outer1,
            sipm_outer2,
            [[0, 0, 0], [0, 0, 0]],
            self.registry,
        )
        self.sipm_outer_bottom_lv_bend = g4.LogicalVolume(
            sipm_outer_bottom,
            self.materials.metal_copper,
            f"sipm_outer_bottom{v_suffix}",
            self.registry,
        )

    def _cached_fiber_volumes(self) -> None:
        """Create solids, logical and physical volumes for the fibers, as specified by the parameters of this instance."""
        v_suffix = f"_l{self.fiber_length}"
        if f"fiber_cl2{v_suffix}" in self.registry.solidDict:
            return

        # create solids
        angle = 2 * np.pi / self.number_of_modules
        dim_cl2 = self.FIBER_DIM
        self.fiber_cl2 = g4.solid.Tubs(
            f"fiber_cl2{v_suffix}",
            self.radius - dim_cl2 / 2,
            self.radius + dim_cl2 / 2,
            self.fiber_length,
            0,
            angle,
            self.registry,
            "mm",
        )
        dim_cl1 = self.FIBER_DIM - self.FIBER_THICKNESS_CL1
        self.fiber_cl1 = g4.solid.Tubs(
            f"fiber_cl1{v_suffix}",
            self.radius - dim_cl1 / 2,
            self.radius + dim_cl1 / 2,
            self.fiber_length,
            0,
            angle,
            self.registry,
            "mm",
        )
        dim_core = self.FIBER_DIM - self.FIBER_THICKNESS_CL1 - self.FIBER_THICKNESS_CL2
        self.fiber_core = g4.solid.Tubs(
            f"fiber_core{v_suffix}",
            self.radius - dim_core / 2,
            self.radius + dim_core / 2,
            self.fiber_length,
            0,
            angle,
            self.registry,
            "mm",
        )
        if self.bend_radius_mm is not None:
            z, r = self._get_bend_polycone(self.radius - dim_cl2 / 2, self.radius + dim_cl2 / 2)
            self.fiber_cl2_bend = g4.solid.GenericPolycone(
                f"fiber_cl2_bend{v_suffix}", 0, angle, r, z, self.registry, "mm"
            )
            z, r = self._get_bend_polycone(self.radius - dim_cl1 / 2, self.radius + dim_cl1 / 2)
            self.fiber_cl1_bend = g4.solid.GenericPolycone(
                f"fiber_cl1_bend{v_suffix}", 0, angle, r, z, self.registry, "mm"
            )
            z, r = self._get_bend_polycone(self.radius - dim_core / 2, self.radius + dim_core / 2)
            self.fiber_core_bend = g4.solid.GenericPolycone(
                f"fiber_core_bend{v_suffix}", 0, angle, r, z, self.registry, "mm"
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
        if self.bend_radius_mm is not None:
            self.fiber_cl2_bend_lv = g4.LogicalVolume(
                self.fiber_cl2_bend,
                self.materials.pmma_out,
                f"fiber_cl2_bend{v_suffix}",
                self.registry,
            )
            self.fiber_cl1_bend_lv = g4.LogicalVolume(
                self.fiber_cl1_bend, self.materials.pmma, f"fiber_cl1_bend{v_suffix}", self.registry
            )
            self.fiber_core_bend_lv = g4.LogicalVolume(
                self.fiber_core_bend,
                self.materials.ps_fibers,
                f"fiber_core_bend{v_suffix}",
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
        if self.bend_radius_mm is not None:
            self.fiber_cl1_bend_pv = g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, self.FIBER_DIM - dim_cl1],
                self.fiber_cl1_bend_lv,
                f"fiber_cl1_bend{v_suffix}",
                self.fiber_cl2_bend_lv,
                self.registry,
            )
            self.fiber_core_bend_pv = g4.PhysicalVolume(
                [0, 0, 0],
                [0, 0, self.FIBER_DIM - dim_core],
                self.fiber_core_bend_lv,
                f"fiber_core_bend{v_suffix}",
                self.fiber_cl1_bend_lv,
                self.registry,
            )

    def _cached_tpb_coating_volume(self, tpb_thickness_nm: float, bend: bool = False) -> g4.LogicalVolume:
        """Create and cache a TPB coating layer of the specified thickness.

        The TPB-Layer is dependent on the module (i.e. the applied thickness varies slightly),
        so we cannot cache it globally on this instance.
        """
        v_suffix = f"{'_bend' if bend else ''}_l{self.fiber_length}_tpb{tpb_thickness_nm}"
        v_name = f"fiber_coating{v_suffix}"
        if v_name in self.registry.solidDict:
            return self.registry.logicalVolumeDict[v_name]

        coating_dim = self.FIBER_DIM + tpb_thickness_nm / 1e6
        if not bend:
            coating = g4.solid.Tubs(
                v_name,
                self.radius - coating_dim / 2,
                self.radius + coating_dim / 2,
                self.fiber_length,
                0,
                2 * math.pi / self.number_of_modules,
                self.registry,
                "mm",
            )
            inner_lv = self.fiber_cl2_lv
            z_displacement = 0
        else:
            angle = 2 * np.pi / self.number_of_modules
            z, r = self._get_bend_polycone(self.radius - coating_dim / 2, self.radius + coating_dim / 2)
            coating = g4.solid.GenericPolycone(v_name, 0, angle, r, z, self.registry, "mm")
            inner_lv = self.fiber_cl2_bend_lv
            z_displacement = self.FIBER_DIM - coating_dim
        coating_lv = g4.LogicalVolume(coating, self.materials.tpb_on_fibers, v_name, self.registry)
        g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, z_displacement],
            inner_lv,
            f"fiber_cl2{v_suffix}",
            coating_lv,
            self.registry,
        )

        coating_lv.pygeom_color_rgba = [0, 1, 0, 1]

        return coating_lv

    def create_module(
        self,
        mod_name: str,
        mod: FiberModuleData,
        mother_lv: g4.LogicalVolume,
        module_num: int,
        z_displacement: float,
    ) -> None:
        if module_num < 0 or module_num >= self.number_of_modules:
            msg = f"invalid module number {module_num} for a maximum of {self.number_of_modules}-1 modules."
            raise ValueError(msg)

        self._cached_fiber_volumes()
        self._cached_sipm_volumes()
        coating_lv = self._cached_tpb_coating_volume(mod.tpb_thickness, bend=False)
        if self.bend_radius_mm is not None:
            self._cached_sipm_volumes_bend()
            coating_lv_bend = self._cached_tpb_coating_volume(mod.tpb_thickness, bend=True)

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
        if self.bend_radius_mm is not None:
            fibers.append(
                g4.PhysicalVolume(
                    [0, 0, -th],
                    [0, 0, z_displacement - self.fiber_length / 2 - self.bend_radius_mm],
                    coating_lv_bend,
                    f"fiber_bend_{mod_name}_s",
                    mother_lv,
                    self.registry,
                )
            )

        # create SiPMs and attach to fibers
        self._create_sipm(module_num, fibers, True, mother_lv, mod.channel_top, z_displacement)
        if self.bend_radius_mm is None:
            self._create_sipm(module_num, fibers, False, mother_lv, mod.channel_bottom, z_displacement)
        else:
            start_angle = 2 * math.pi / self.number_of_modules * module_num
            z = z_displacement - self.fiber_length / 2 - self.bend_radius_mm
            g4.PhysicalVolume(
                [0, 0, -start_angle],
                [0, 0, z],
                self.sipm_lv_bend,
                mod.channel_bottom,
                mother_lv,
                self.registry,
            )
            g4.PhysicalVolume(
                [0, 0, -start_angle],
                [0, 0, z],
                self.sipm_outer_bottom_lv_bend,
                f"{mod.channel_bottom}_wrap",
                mother_lv,
                self.registry,
            )
