from __future__ import annotations

import argparse
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class GeometryParams:
    gap: float
    T: float
    H_wall: float
    W_wall: float
    p: float
    t: float
    L: float


@dataclass(frozen=True)
class MeshingParams:
    background_cell_size: float = 0.005
    target_surface_cell_size: float | None = None
    max_local_cells: int = 2_000_000
    max_global_cells: int = 10_000_000
    min_refinement_cells: int = 10
    n_cells_between_levels: int = 2
    resolve_feature_angle_deg: float = 30.0
    # Feature snapping is enabled by default because it is required to preserve
    # sharp fin-wall edges reliably for this geometry.
    # NOTE: When using --write-only, you must run the feature tool before snappyHexMesh.
    feature_snap: bool = True
    feature_tool: str = "surfaceFeatures"  # or "surfaceFeatureExtract"
    feature_extract_included_angle_deg: float = 150.0


def _expected_stl_files(tri_surface_dir: Path) -> list[Path]:
    """STLs that this project generates and snappy should consume."""

    stls = [tri_surface_dir / "solid_left.stl", tri_surface_dir / "solid_right.stl"]
    missing = [p for p in stls if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing expected STL files: " + ", ".join(str(p) for p in missing)
        )
    return stls


def _max_count_with_gap(total_length: float, thickness: float, gap: float) -> int:
    if total_length <= 0:
        return 0
    if thickness <= 0:
        return 0
    if gap < 0:
        return 0
    return int(math.floor((total_length + gap) / (thickness + gap) + 1e-12))


def _blockmesh_cell_count(length: float, cell_size: float, min_cells: int = 10) -> int:
    if length <= 0:
        return min_cells
    if cell_size <= 0:
        return min_cells
    return max(min_cells, int(math.ceil(length / cell_size)))


def generate_blockmeshdict_text(geom: GeometryParams, mesh: MeshingParams) -> str:
    if geom.L > geom.gap:
        raise ValueError(
            "Invalid geometry for fin protrusion length: require L <= gap to avoid penetrating the opposite wall."
        )

    if geom.H_wall <= 0 or geom.W_wall <= 0 or geom.T <= 0 or geom.gap <= 0:
        raise ValueError("Invalid geometry: gap, T, H_wall, and W_wall must be > 0")

    num_y = _max_count_with_gap(geom.H_wall, geom.t, geom.p)
    num_z = _max_count_with_gap(geom.W_wall, geom.t, geom.p)
    used_y = num_y * geom.t + (num_y - 1) * geom.p if num_y > 0 else 0.0
    used_z = num_z * geom.t + (num_z - 1) * geom.p if num_z > 0 else 0.0
    y_start = (geom.H_wall - used_y) / 2 if geom.H_wall > 0 else 0.0
    z_start = (geom.W_wall - used_z) / 2 if geom.W_wall > 0 else 0.0
    pitch = geom.t + geom.p

    def _unique_sorted(vals: list[float], eps: float = 1e-12) -> list[float]:
        out: list[float] = []
        for v in sorted(vals):
            if not out or abs(v - out[-1]) > eps:
                out.append(v)
        return out

    x_planes = _unique_sorted(
        [
            -geom.gap / 2 - geom.T,
            -geom.gap / 2,
            -geom.gap / 2 + geom.L,
            geom.gap / 2 - geom.L,
            geom.gap / 2,
            geom.gap / 2 + geom.T,
        ]
    )

    y_planes_raw = [0.0, geom.H_wall]
    z_planes_raw = [0.0, geom.W_wall]
    if num_y > 0 and num_z > 0:
        for iy in range(num_y):
            y0 = y_start + iy * pitch
            y1 = y0 + geom.t
            y_planes_raw.extend([y0, y1])
        for iz in range(num_z):
            z0 = z_start + iz * pitch
            z1 = z0 + geom.t
            z_planes_raw.extend([z0, z1])

    y_planes = _unique_sorted(y_planes_raw)
    z_planes = _unique_sorted(z_planes_raw)

    nxp = len(x_planes)
    nyp = len(y_planes)
    nzp = len(z_planes)
    if nxp < 2 or nyp < 2 or nzp < 2:
        raise ValueError("Failed to build valid blockMesh plane set")

    def v_idx(ix: int, iy: int, iz: int) -> int:
        return ix + nxp * (iy + nyp * iz)

    def _tile_side(y: float, z: float) -> str | None:
        if num_y <= 0 or num_z <= 0 or pitch <= 0:
            return None
        if not (y_start <= y <= y_start + used_y and z_start <= z <= z_start + used_z):
            return None

        iy = int((y - y_start) / pitch)
        iz = int((z - z_start) / pitch)
        iy = min(max(iy, 0), num_y - 1)
        iz = min(max(iz, 0), num_z - 1)

        y0 = y_start + iy * pitch
        z0 = z_start + iz * pitch
        in_fin = (y0 <= y < y0 + geom.t) and (z0 <= z < z0 + geom.t)
        if not in_fin:
            return None
        return "left" if (iy + iz) % 2 == 0 else "right"

    def _zone_for_cell(x_mid: float, y_mid: float, z_mid: float) -> str:
        left_wall_inner = -geom.gap / 2
        right_wall_inner = geom.gap / 2
        left_tip = -geom.gap / 2 + geom.L
        right_tip = geom.gap / 2 - geom.L

        if x_mid <= left_wall_inner:
            return "solid_left"
        if x_mid >= right_wall_inner:
            return "solid_right"

        side = _tile_side(y_mid, z_mid)
        if side == "left" and x_mid <= left_tip:
            return "solid_left"
        if side == "right" and x_mid >= right_tip:
            return "solid_right"
        return "fluid"

    vertex_lines: list[str] = []
    vid = 0
    for iz in range(nzp):
        z = z_planes[iz]
        for iy in range(nyp):
            y = y_planes[iy]
            for ix in range(nxp):
                x = x_planes[ix]
                vertex_lines.append(f"    ({x:.9f} {y:.9f} {z:.9f})  // {vid}")
                vid += 1

    block_lines: list[str] = []
    for ix in range(nxp - 1):
        x0 = x_planes[ix]
        x1 = x_planes[ix + 1]
        nx = _blockmesh_cell_count(x1 - x0, mesh.background_cell_size, min_cells=1)
        x_mid = 0.5 * (x0 + x1)
        for iy in range(nyp - 1):
            y0 = y_planes[iy]
            y1 = y_planes[iy + 1]
            ny = _blockmesh_cell_count(y1 - y0, mesh.background_cell_size, min_cells=1)
            y_mid = 0.5 * (y0 + y1)
            for iz in range(nzp - 1):
                z0 = z_planes[iz]
                z1 = z_planes[iz + 1]
                nz = _blockmesh_cell_count(
                    z1 - z0, mesh.background_cell_size, min_cells=1
                )
                z_mid = 0.5 * (z0 + z1)

                zone = _zone_for_cell(x_mid, y_mid, z_mid)
                v000 = v_idx(ix, iy, iz)
                v100 = v_idx(ix + 1, iy, iz)
                v110 = v_idx(ix + 1, iy + 1, iz)
                v010 = v_idx(ix, iy + 1, iz)
                v001 = v_idx(ix, iy, iz + 1)
                v101 = v_idx(ix + 1, iy, iz + 1)
                v111 = v_idx(ix + 1, iy + 1, iz + 1)
                v011 = v_idx(ix, iy + 1, iz + 1)
                block_lines.append(
                    "    hex "
                    + f"({v000} {v100} {v110} {v010} {v001} {v101} {v111} {v011}) "
                    + f"{zone} ({nx} {ny} {nz}) simpleGrading (1 1 1)"
                )

    faces: list[tuple[int, int, int, int]] = []

    # x-min / x-max
    for iy in range(nyp - 1):
        for iz in range(nzp - 1):
            faces.append(
                (
                    v_idx(0, iy, iz),
                    v_idx(0, iy, iz + 1),
                    v_idx(0, iy + 1, iz + 1),
                    v_idx(0, iy + 1, iz),
                )
            )
            faces.append(
                (
                    v_idx(nxp - 1, iy + 1, iz),
                    v_idx(nxp - 1, iy + 1, iz + 1),
                    v_idx(nxp - 1, iy, iz + 1),
                    v_idx(nxp - 1, iy, iz),
                )
            )

    # y-min / y-max
    for ix in range(nxp - 1):
        for iz in range(nzp - 1):
            faces.append(
                (
                    v_idx(ix + 1, 0, iz),
                    v_idx(ix + 1, 0, iz + 1),
                    v_idx(ix, 0, iz + 1),
                    v_idx(ix, 0, iz),
                )
            )
            faces.append(
                (
                    v_idx(ix, nyp - 1, iz),
                    v_idx(ix, nyp - 1, iz + 1),
                    v_idx(ix + 1, nyp - 1, iz + 1),
                    v_idx(ix + 1, nyp - 1, iz),
                )
            )

    # z-min / z-max
    for ix in range(nxp - 1):
        for iy in range(nyp - 1):
            faces.append(
                (
                    v_idx(ix, iy, 0),
                    v_idx(ix, iy + 1, 0),
                    v_idx(ix + 1, iy + 1, 0),
                    v_idx(ix + 1, iy, 0),
                )
            )
            faces.append(
                (
                    v_idx(ix, iy, nzp - 1),
                    v_idx(ix + 1, iy, nzp - 1),
                    v_idx(ix + 1, iy + 1, nzp - 1),
                    v_idx(ix, iy + 1, nzp - 1),
                )
            )

    face_lines = [f"            ({a} {b} {c} {d})" for (a, b, c, d) in faces]
    vertices_block = "\n".join(vertex_lines)
    blocks_block = "\n".join(block_lines)
    boundary_faces_block = "\n".join(face_lines)

    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// Generated by openfoam_automation.py

scale 1;

vertices
(
{vertices_block}
);

blocks
(
{blocks_block}
);

edges
();

boundary
(
    background
    {{
        type patch;
        faces
        (
{boundary_faces_block}
        );
    }}
);

mergePatchPairs
();

// ************************************************************************* //
"""


def _surface_refinement_level(
    background_cell_size: float, target_cell_size: float
) -> int:
    if background_cell_size <= 0 or target_cell_size <= 0:
        return 0
    ratio = background_cell_size / target_cell_size
    if ratio <= 1.0:
        return 0
    return int(math.ceil(math.log(ratio, 2)))


def _safe_location_in_mesh(
    geom: GeometryParams, mesh: MeshingParams
) -> tuple[float, float, float]:
    """Pick a point strictly inside the kept *fluid* region.

    Requirements:
    - Must be inside the blockMesh bounds.
    - Must be outside all solid STL volumes (walls and fins).
    - Must not lie on/too close to surfaces/mesh boundaries (numerical stability).
    """

    def _clamp(v: float, lo: float, hi: float) -> float:
        return min(max(v, lo), hi)

    def _grid_params() -> tuple[int, int, float, float, float, float, float]:
        num_y = _max_count_with_gap(geom.H_wall, geom.t, geom.p)
        num_z = _max_count_with_gap(geom.W_wall, geom.t, geom.p)
        used_y = num_y * geom.t + (num_y - 1) * geom.p if num_y > 0 else 0.0
        used_z = num_z * geom.t + (num_z - 1) * geom.p if num_z > 0 else 0.0
        y_start = (geom.H_wall - used_y) / 2 if geom.H_wall > 0 else 0.0
        z_start = (geom.W_wall - used_z) / 2 if geom.W_wall > 0 else 0.0
        pitch = geom.t + geom.p
        return num_y, num_z, used_y, used_z, y_start, z_start, pitch

    def _cell_indices(
        y: float, z: float, y_start: float, z_start: float, pitch: float
    ) -> tuple[int, int]:
        # y/z are expected to be within the used band when this is called.
        iy = int((y - y_start) / pitch)
        iz = int((z - z_start) / pitch)
        return iy, iz

    def _in_used_band(
        y: float, z: float, y_start: float, z_start: float, used_y: float, used_z: float
    ) -> bool:
        return (y_start <= y <= y_start + used_y) and (z_start <= z <= z_start + used_z)

    def _inside_fin_footprint(
        y: float,
        z: float,
        y_start: float,
        z_start: float,
        used_y: float,
        used_z: float,
        pitch: float,
    ) -> tuple[bool, int, int]:
        if pitch <= 0:
            return False, -1, -1
        if not _in_used_band(y, z, y_start, z_start, used_y, used_z):
            return False, -1, -1
        iy, iz = _cell_indices(y, z, y_start, z_start, pitch)
        # Guard against y/z at the upper boundary due to float rounding.
        iy = min(max(iy, 0), max(0, num_y - 1))
        iz = min(max(iz, 0), max(0, num_z - 1))
        y0 = y_start + iy * pitch
        z0 = z_start + iz * pitch
        in_y = (y0 <= y) and (y < y0 + geom.t)
        in_z = (z0 <= z) and (z < z0 + geom.t)
        return (in_y and in_z), iy, iz

    def _point_inside_wall_solid(x: float, y: float, z: float) -> bool:
        if not (0.0 <= y <= geom.H_wall and 0.0 <= z <= geom.W_wall):
            return False
        if (-geom.gap / 2 - geom.T) <= x <= (-geom.gap / 2):
            return True
        if (geom.gap / 2) <= x <= (geom.gap / 2 + geom.T):
            return True
        return False

    def _point_inside_fin_solid(
        x: float,
        y: float,
        z: float,
        y_start: float,
        z_start: float,
        used_y: float,
        used_z: float,
        pitch: float,
    ) -> bool:
        in_fp, iy, iz = _inside_fin_footprint(
            y, z, y_start, z_start, used_y, used_z, pitch
        )
        if not in_fp:
            return False
        if (iy + iz) % 2 == 0:
            # Left fin.
            return (-geom.gap / 2) <= x <= (-geom.gap / 2 + geom.L)
        # Right fin.
        return (geom.gap / 2 - geom.L) <= x <= (geom.gap / 2)

    num_y, num_z, used_y, used_z, y_start, z_start, pitch = _grid_params()

    # Conservative epsilon: keep away from mesh boundaries and surfaces.
    # OpenFOAM guidance is to stay away from faces/edges; 10-20% of cell size is a good rule.
    cell_eps = max(1e-6, 0.2 * mesh.background_cell_size)
    x_lo = -geom.gap / 2 + cell_eps
    x_hi = geom.gap / 2 - cell_eps
    y_lo = 0.0 + cell_eps
    y_hi = geom.H_wall - cell_eps
    z_lo = 0.0 + cell_eps
    z_hi = geom.W_wall - cell_eps

    if not (x_lo < x_hi and y_lo < y_hi and z_lo < z_hi):
        raise ValueError(
            "No room for a safe locationInMesh point: domain too small vs background_cell_size; reduce background_cell_size or increase dimensions."
        )

    # Candidate sets for y and z.
    y_candidates: list[float] = []
    z_candidates: list[float] = []

    # 1) Prefer a y/z gap between fins.
    if geom.p > 0 and pitch > 0:
        if num_y >= 2:
            y_candidates.append(y_start + geom.t + 0.5 * geom.p)
        if num_z >= 2:
            z_candidates.append(z_start + geom.t + 0.5 * geom.p)

    # 2) Use leftover margins outside the fin band (if any).
    y_margin_bot = y_start
    y_margin_top = geom.H_wall - (y_start + used_y)
    # Only use margin candidates if they can stay at least cell_eps away from the fin band.
    if y_margin_bot > 2 * cell_eps:
        y_candidates.append(0.5 * y_margin_bot)
    if y_margin_top > 2 * cell_eps:
        y_candidates.append(geom.H_wall - 0.5 * y_margin_top)

    z_margin_bot = z_start
    z_margin_top = geom.W_wall - (z_start + used_z)
    if z_margin_bot > 2 * cell_eps:
        z_candidates.append(0.5 * z_margin_bot)
    if z_margin_top > 2 * cell_eps:
        z_candidates.append(geom.W_wall - 0.5 * z_margin_top)

    # 3) Always include centered values as a last resort.
    y_candidates.append(geom.H_wall / 2)
    z_candidates.append(geom.W_wall / 2)

    # Deduplicate while preserving order.
    def _uniq(vals: list[float]) -> list[float]:
        out: list[float] = []
        for v in vals:
            if all(abs(v - w) > 1e-12 for w in out):
                out.append(v)
        return out

    y_candidates = _uniq([_clamp(y, y_lo, y_hi) for y in y_candidates])
    z_candidates = _uniq([_clamp(z, z_lo, z_hi) for z in z_candidates])

    clearance = geom.gap - geom.L
    for y_cand in y_candidates:
        for z_cand in z_candidates:
            y_try = y_cand
            z_try = z_cand
            # If (y,z) is not in any fin footprint, x=0 is safe (even when fins cross at x=0).
            in_fp, iy, iz = _inside_fin_footprint(
                y_try, z_try, y_start, z_start, used_y, used_z, pitch
            )
            if not in_fp:
                x = 0.0
            else:
                # Move (y,z) away from fin footprint boundaries.
                if geom.t <= 2 * cell_eps:
                    continue

                y0 = y_start + iy * pitch
                z0 = z_start + iz * pitch
                y_try = _clamp(y0 + 0.5 * geom.t, y_lo, y_hi)
                z_try = _clamp(z0 + 0.5 * geom.t, z_lo, z_hi)

                if clearance <= 2 * cell_eps:
                    # In a fully tiled y-z plane and L ~ gap, there may be no fluid region.
                    continue
                if (iy + iz) % 2 == 0:
                    # Left fin exists in this cell: pick x near the right wall inner face.
                    x = x_hi
                else:
                    # Right fin exists in this cell: pick x near the left wall inner face.
                    x = x_lo

            # Validate point is outside solids.
            if _point_inside_wall_solid(x, y_try, z_try):
                continue
            if _point_inside_fin_solid(
                x, y_try, z_try, y_start, z_start, used_y, used_z, pitch
            ):
                continue

            return (x, y_try, z_try)

    raise ValueError(
        "Failed to find a safe locationInMesh point in the kept fluid region; try increasing p (create y/z gaps), or ensure gap > L (provide clearance), or adjust H_wall/W_wall so fins do not tile the entire cross-section."
    )


def _safe_inside_point_in_solid(
    geom: GeometryParams, mesh: MeshingParams, *, side: str
) -> tuple[float, float, float]:
    """Pick a point strictly inside a wall solid volume (not on any faces)."""

    if side not in {"left", "right"}:
        raise ValueError("side must be 'left' or 'right'")

    cell_eps = max(1e-6, 0.2 * mesh.background_cell_size)
    if geom.T <= 2 * cell_eps:
        raise ValueError(
            "No safe solid inside point: wall thickness T too small vs background_cell_size; reduce background_cell_size or increase T."
        )

    if side == "left":
        x = -geom.gap / 2 - 0.5 * geom.T
        x = min(max(x, -geom.gap / 2 - geom.T + cell_eps), -geom.gap / 2 - cell_eps)
    else:
        x = geom.gap / 2 + 0.5 * geom.T
        x = min(max(x, geom.gap / 2 + cell_eps), geom.gap / 2 + geom.T - cell_eps)

    y = min(max(0.5 * geom.H_wall, cell_eps), geom.H_wall - cell_eps)
    z = min(max(0.5 * geom.W_wall, cell_eps), geom.W_wall - cell_eps)

    return (x, y, z)


def generate_snappyhexmeshdict_text(
    geom: GeometryParams, mesh: MeshingParams, tri_surface_dir: Path
) -> str:
    stl_files = _expected_stl_files(tri_surface_dir)

    target_size = mesh.target_surface_cell_size
    if target_size is None:
        target_size = max(geom.t / 2.0, 1e-6)

    level = _surface_refinement_level(mesh.background_cell_size, target_size)
    level_pair = f"({level} {level})"

    fluid_location = _safe_location_in_mesh(geom, mesh)
    solid_left_location = _safe_inside_point_in_solid(geom, mesh, side="left")
    solid_right_location = _safe_inside_point_in_solid(geom, mesh, side="right")

    inside_points = [fluid_location, solid_left_location, solid_right_location]
    inside_points_str = (
        "(\n"
        + "\n".join(
            [f"        ({p[0]:.10f} {p[1]:.10f} {p[2]:.10f})" for p in inside_points]
        )
        + "\n    )"
    )

    geometry_entries: list[str] = []
    refinement_entries: list[str] = []
    feature_entries: list[str] = []
    for f in stl_files:
        name = f.stem
        geometry_entries.append(
            f'    {f.name}\n    {{\n        type triSurfaceMesh;\n        file "{f.name}";\n        name {name};\n    }}\n'
        )
        if name.startswith("solid_"):
            refinement_entries.append(
                "    "
                + name
                + "\n    {\n        level "
                + level_pair
                + ";\n        faceZone "
                + name
                + ";\n        cellZone "
                + name
                + ";\n        mode inside;\n    }\n"
            )
        else:
            refinement_entries.append(
                "    "
                + name
                + "\n    {\n        level "
                + level_pair
                + ";\n        patchInfo { type wall; }\n    }\n"
            )
        if mesh.feature_snap:
            ext = (
                "extendedFeatureEdgeMesh"
                if mesh.feature_tool == "surfaceFeatures"
                else "eMesh"
            )
            feature_entries.append(
                f'        {{ file "{name}.{ext}"; level {level}; }}\n'
            )

    geometry_block = "".join(geometry_entries)
    refinement_block = "".join(refinement_entries)
    features_block = "".join(feature_entries)

    snap_feature_block = ""
    if mesh.feature_snap:
        # Note: In OpenFOAM, leaving out nFeatureSnapIter disables feature snapping.
        snap_feature_block = (
            "    nFeatureSnapIter 15;\n"
            "    implicitFeatureSnap false;\n"
            "    explicitFeatureSnap true;\n"
            "    multiRegionFeatureSnap true;\n"
        )

    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}
// Generated by openfoam_automation.py

castellatedMesh true;
snap            true;
addLayers       false;

geometry
{{
{geometry_block}}}

castellatedMeshControls
{{
    maxLocalCells {mesh.max_local_cells};
    maxGlobalCells {mesh.max_global_cells};
    minRefinementCells {mesh.min_refinement_cells};
    nCellsBetweenLevels {mesh.n_cells_between_levels};

    features
    (
{features_block}    );

    refinementSurfaces
    {{
{refinement_block}    }}

    resolveFeatureAngle {mesh.resolve_feature_angle_deg:.1f};

    refinementRegions
    {{
    }}

    // Keep multiple disconnected regions (fluid + solids)
    insidePoints {inside_points_str};
    allowFreeStandingZoneFaces true;
}}

snapControls
{{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
{snap_feature_block}
}}

addLayersControls
{{
}}

meshQualityControls
{{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-9;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
}}

writeFlags
(
    scalarLevels
    layerSets
    layerFields
);

mergeTolerance 1e-6;

// ************************************************************************* //
"""


def generate_controldict_text() -> str:
    # Multi-region controlDict template (OpenFOAM.org dev):
    # - no 'application' entry
    # - foamSetupCHT generates system/regionSolvers
    return """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}
// Generated by openfoam_automation.py

#includeIfPresent "regionSolvers"

startFrom       startTime;
startTime       0;

stopAt          endTime;
endTime         500;

deltaT          1;

writeControl    timeStep;
writeInterval   50;

purgeWrite      5;

writeFormat     ascii;
writePrecision  8;
writeCompression off;

timeFormat      general;
timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""


def generate_fvschemes_text() -> str:
    return """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// Generated by openfoam_automation.py

ddtSchemes
{
}

gradSchemes
{
}

divSchemes
{
}

laplacianSchemes
{
}

interpolationSchemes
{
}

snGradSchemes
{
}

fluxRequired
{
}

// ************************************************************************* //
"""


def generate_fvsolution_text() -> str:
    return """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// Generated by openfoam_automation.py

PIMPLE
{
    nOuterCorrectors 1;
}

// ************************************************************************* //
"""


def generate_materialproperties_text(
    *,
    fluid_material: str = "air",
    solid_left_material: str = "copper",
    solid_right_material: str = "copper",
) -> str:
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      materialProperties;
}}
// Generated by openfoam_automation.py

fluid
{{
    solver          fluid;
    material        {fluid_material};
}}

solid_left
{{
    solver          solid;
    material        {solid_left_material};
}}

solid_right
{{
    solver          solid;
    material        {solid_right_material};
}}

// ************************************************************************* //
"""


def generate_surfacefeatureextractdict_text(
    tri_surface_dir: Path, mesh: MeshingParams
) -> str:
    stl_files = _expected_stl_files(tri_surface_dir)

    entries: list[str] = []
    for f in stl_files:
        entry = (
            f"{f.name}\n"
            + "{\n"
            + "    extractionMethod    extractFromSurface;\n"
            + "\n"
            + "    extractFromSurfaceCoeffs\n"
            + "    {\n"
            + f"        includedAngle   {mesh.feature_extract_included_angle_deg:.1f};\n"
            + "    }\n"
            + "\n"
            + "    writeObj            yes;\n"
            + "}\n\n"
        )
        entries.append(entry)

    entries_block = "".join(entries)

    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      surfaceFeatureExtractDict;
}}
// Generated by openfoam_automation.py

{entries_block}// ************************************************************************* //
"""


def generate_surfacefeaturesdict_text(
    tri_surface_dir: Path, mesh: MeshingParams
) -> str:
    stl_files = _expected_stl_files(tri_surface_dir)

    surfaces_block = "\n".join([f'    "{f.name}"' for f in stl_files])

    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      surfaceFeaturesDict;
}}
// Generated by openfoam_automation.py

surfaces
(
{surfaces_block}
);

// Identify a feature when angle between faces < includedAngle
includedAngle   {mesh.feature_extract_included_angle_deg:.1f};

writeObj        yes;

// ************************************************************************* //
"""


def write_openfoam_case_files(
    case_dir: Path, geom: GeometryParams, mesh: MeshingParams
) -> None:
    system_dir = case_dir / "system"
    constant_dir = case_dir / "constant"
    system_dir.mkdir(parents=True, exist_ok=True)
    constant_dir.mkdir(parents=True, exist_ok=True)

    # Remove stale snappy/feature dictionaries from prior workflows.
    for dict_name in [
        "surfaceFeaturesDict",
        "surfaceFeatureExtractDict",
        "snappyHexMeshDict",
    ]:
        try:
            (system_dir / dict_name).unlink()
        except FileNotFoundError:
            pass

    # Remove stale foamSetupCHT output.
    try:
        (system_dir / "regionSolvers").unlink()
    except FileNotFoundError:
        pass

    (system_dir / "controlDict").write_text(
        generate_controldict_text(), encoding="utf-8"
    )
    (system_dir / "fvSchemes").write_text(generate_fvschemes_text(), encoding="utf-8")
    (system_dir / "fvSolution").write_text(generate_fvsolution_text(), encoding="utf-8")
    (system_dir / "blockMeshDict").write_text(
        generate_blockmeshdict_text(geom, mesh), encoding="utf-8"
    )

    (constant_dir / "materialProperties").write_text(
        generate_materialproperties_text(), encoding="utf-8"
    )

    # Helper script for blockMesh-only + multi-region splitting workflow.
    allmesh = (
        "#!/bin/sh\n"
        "set -eu\n"
        "blockMesh > log.blockMesh 2>&1\n"
        "splitMeshRegions -cellZones -overwrite -defaultRegionName fluid > log.splitMeshRegions 2>&1\n"
    )
    allmesh_path = case_dir / "Allmesh"
    allmesh_path.write_text(allmesh, encoding="utf-8")
    try:
        allmesh_path.chmod(0o755)
    except OSError:
        # Best-effort; e.g. on some filesystems chmod may be disallowed.
        pass

    # Optional helper for a full CHT run (OpenFOAM.org dev workflow).
    allrun = (
        "#!/bin/sh\n"
        "set -eu\n"
        "./Allmesh\n"
        "foamSetupCHT > log.foamSetupCHT 2>&1\n"
        "foamMultiRun > log.foamMultiRun 2>&1\n"
        "paraFoam -touchAll > /dev/null 2>&1 || true\n"
    )
    allrun_path = case_dir / "Allrun"
    allrun_path.write_text(allrun, encoding="utf-8")
    try:
        allrun_path.chmod(0o755)
    except OSError:
        pass


def run_openfoam_meshing(case_dir: Path, mesh: MeshingParams) -> None:
    if shutil.which("blockMesh") is None:
        raise FileNotFoundError(
            "blockMesh not found in PATH (did you source OpenFOAM?)"
        )
    if shutil.which("splitMeshRegions") is None:
        raise FileNotFoundError(
            "splitMeshRegions not found in PATH (did you source OpenFOAM?)"
        )

    subprocess.run(["blockMesh"], cwd=str(case_dir), check=True)
    subprocess.run(
        [
            "splitMeshRegions",
            "-cellZones",
            "-overwrite",
            "-defaultRegionName",
            "fluid",
        ],
        cwd=str(case_dir),
        check=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate OpenFOAM blockMeshDict for checkerboard fin walls and run blockMesh + splitMeshRegions."
    )
    parser.add_argument("--case-dir", default=".")
    parser.add_argument("--write-only", action="store_true")

    parser.add_argument("--gap", type=float, default=0.20)
    parser.add_argument("--T", type=float, default=0.01)
    parser.add_argument("--H-wall", type=float, default=0.10)
    parser.add_argument("--W-wall", type=float, default=0.10)
    parser.add_argument("--p", type=float, default=0.02)
    parser.add_argument("--t", type=float, default=0.005)
    parser.add_argument("--L", type=float, default=0.12)

    parser.add_argument("--background-cell-size", type=float, default=0.0005)

    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    geom = GeometryParams(
        gap=args.gap,
        T=args.T,
        H_wall=args.H_wall,
        W_wall=args.W_wall,
        p=args.p,
        t=args.t,
        L=args.L,
    )
    mesh = MeshingParams(
        background_cell_size=args.background_cell_size,
    )

    write_openfoam_case_files(case_dir, geom, mesh)
    if not args.write_only:
        run_openfoam_meshing(case_dir, mesh)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
