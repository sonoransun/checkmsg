"""3D voxel grid carrying material identity + density for muon forward projection."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from checkmsg.refdata.muon_data import Material, get_material


@dataclass
class VoxelGrid:
    """A Cartesian voxel grid representing a 3D composite subject.

    Each voxel carries:
      - a `Material` (class + bulk properties Z_eff, A_eff, X_0, I)
      - a density in g/cc (overrides the material's bulk density; 0 means void)

    All physical positions are expressed in millimetres. The grid origin is the
    centre of voxel (0, 0, 0).
    """

    materials: np.ndarray  # dtype=object, shape (Nx, Ny, Nz)
    densities: np.ndarray  # dtype=float, shape (Nx, Ny, Nz), g/cc
    spacing_mm: tuple[float, float, float] = (1.0, 1.0, 1.0)
    origin_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.materials.shape != self.densities.shape:
            raise ValueError(
                f"materials shape {self.materials.shape} does not match densities {self.densities.shape}"
            )
        if self.materials.ndim != 3:
            raise ValueError("VoxelGrid requires 3-D arrays")

    @classmethod
    def filled(cls, shape: tuple[int, int, int], material_name: str = "vacuum",
               *, spacing_mm: tuple[float, float, float] = (1.0, 1.0, 1.0),
               origin_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)) -> VoxelGrid:
        """Construct a uniform grid filled with a single material at its bulk density."""
        material = get_material(material_name)
        materials = np.empty(shape, dtype=object)
        materials.fill(material)
        densities = np.full(shape, material.density_g_cc, dtype=float)
        return cls(materials=materials, densities=densities,
                   spacing_mm=spacing_mm, origin_mm=origin_mm)

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.materials.shape  # type: ignore[return-value]

    def bounds_mm(self) -> tuple[np.ndarray, np.ndarray]:
        """Lower and upper corners of the grid in millimetres."""
        spacing = np.asarray(self.spacing_mm, dtype=float)
        origin = np.asarray(self.origin_mm, dtype=float)
        nx, ny, nz = self.shape
        lower = origin - 0.5 * spacing
        upper = origin + (np.array([nx, ny, nz]) - 0.5) * spacing
        return lower, upper

    def set_box(self, lo: tuple[int, int, int], hi: tuple[int, int, int],
                material_name: str, density_g_cc: float | None = None) -> None:
        """Fill an axis-aligned box [lo, hi) (inclusive of lo, exclusive of hi) with a material."""
        material = get_material(material_name)
        rho = material.density_g_cc if density_g_cc is None else density_g_cc
        ix, iy, iz = slice(lo[0], hi[0]), slice(lo[1], hi[1]), slice(lo[2], hi[2])
        self.materials[ix, iy, iz] = material
        self.densities[ix, iy, iz] = rho

    def set_sphere(self, center_voxel: tuple[int, int, int], radius_voxels: float,
                   material_name: str, density_g_cc: float | None = None) -> None:
        """Fill a sphere centred on the given voxel index with a material."""
        material = get_material(material_name)
        rho = material.density_g_cc if density_g_cc is None else density_g_cc
        cx, cy, cz = center_voxel
        nx, ny, nz = self.shape
        xs = np.arange(nx)[:, None, None]
        ys = np.arange(ny)[None, :, None]
        zs = np.arange(nz)[None, None, :]
        dist2 = (xs - cx) ** 2 + (ys - cy) ** 2 + (zs - cz) ** 2
        mask = dist2 <= radius_voxels ** 2
        self.materials[mask] = material
        self.densities[mask] = rho

    def density_array(self) -> np.ndarray:
        """Plain density array (g/cc), useful for transmission / scatter-density baselines."""
        return self.densities.copy()

    def x0_mass_array(self) -> np.ndarray:
        """For each voxel, return X_0 (g/cm²); used to weight Highland scattering."""
        out = np.empty(self.shape, dtype=float)
        for ix in range(self.shape[0]):
            for iy in range(self.shape[1]):
                for iz in range(self.shape[2]):
                    m: Material = self.materials[ix, iy, iz]
                    out[ix, iy, iz] = m.X_0_g_cm2
        return out

    def z_eff_array(self) -> np.ndarray:
        """Effective atomic number per voxel — useful for scattering-density visualisation."""
        out = np.empty(self.shape, dtype=float)
        for ix in range(self.shape[0]):
            for iy in range(self.shape[1]):
                for iz in range(self.shape[2]):
                    m: Material = self.materials[ix, iy, iz]
                    out[ix, iy, iz] = m.Z_eff
        return out
